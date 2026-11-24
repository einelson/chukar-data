"""
Builds Chukar's NBM tiles: NOAA National Blend of Models (2.5 km, lower 48) 2 m temperature and 10 m wind, in the same
2° tiles and format as the HRRR tiles (hrrr/build_hrrr.py), published beside them as a third set, `nbm`. Older app
builds only read the `hourly` and `main` sets, so they never see these.

NBM is NOAA's statistically corrected blend: scored against 4,000 stations over 30 days (the app repo's
docs/WEATHER_ACCURACY_2026-11-NBM.md) its wind beat HRRR's everywhere and its temperature 12 hours out. The app blends
the two.

Reads one run from NOAA's open data bucket on AWS (noaa-nbm-grib2-pds), only the records needed, by byte range. NBM's
own ground height comes from URMA (noaa-urma-pds), which is on the same 2.5 km grid.

Set:
  nbm  the latest 00/06/12/18Z run: hours 1-36 every hour, then every 3 hours to 72

Output, in --out (default build/):
  nbm/<run>/<tile>.zip   hrrr.json (the metadata) + hrrr.bin (the grids) + ground.bin (NBM's ground height)
  nbm.json               the index: run, times, the tiles and where they are

The zip's files keep the HRRR names so one reader serves both. hrrr.bin, per forecast time in order: temperature, then
u, then v, each ny rows from north to south of nx int8 values in steps of `scale` (°C or m/s; u toward the east, v
toward the north; NBM's wind direction is already from true north). -128 means no data. ground.bin: ny x nx
little-endian int16 meters, rows north to south.

Usage: python build_nbm.py [--run 2026112412] [--out build] [--if-newer-than URL]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pygrib
import pyproj

BUCKET = "https://noaa-nbm-grib2-pds.s3.amazonaws.com"
URMA = "https://noaa-urma-pds.s3.amazonaws.com"
FORMAT = 1
SET = "nbm"
SCALE = 0.5
NO_DATA = -128
TILE_DEG = 2
STEP = 0.025
PER_SIDE = int(round(TILE_DEG / STEP)) + 1
# The same lower-48 box as the HRRR tiles: 30 x 13 tiles.
WEST, EAST, SOUTH, NORTH = -126, -66, 24, 50
HOURS = list(range(1, 37)) + list(range(39, 73, 3))
WANTED = ["TMP:2 m above ground", "WIND:10 m above ground", "WDIR:10 m above ground"]
USER_AGENT = "Chukar NBM builder (github.com/einelson/chukar-data)"


def get(url: str, headers: dict | None = None, tries: int = 4) -> bytes:
    for attempt in range(tries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as err:
            if err.code == 404 or attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
        time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def step_url(run: dt.datetime, hour: int) -> str:
    return f"{BUCKET}/blend.{run:%Y%m%d}/{run:%H}/core/blend.t{run:%H}z.core.f{hour:03d}.co.grib2"


def latest_run(now: dt.datetime) -> dt.datetime:
    """The newest 00/06/12/18Z run whose last hour is out."""
    hour = now.replace(minute=0, second=0, microsecond=0)
    for back in range(24):
        run = hour - dt.timedelta(hours=back)
        if run.hour % 6:
            continue
        try:
            get(step_url(run, HOURS[-1]) + ".idx")
            return run
        except urllib.error.HTTPError:
            continue
    raise SystemExit("No complete NBM run in the last 24 hours")


def records(url: str, wanted: list[str]) -> dict[str, pygrib.gribmessage]:
    """The deterministic records (not the ensemble spread or probabilities) named in `wanted`."""
    lines = get(url + ".idx").decode().splitlines()
    starts = [int(line.split(":")[1]) for line in lines]
    found = {}
    for i, line in enumerate(lines):
        fields = line.split(":")
        key = f"{fields[3]}:{fields[4]}"
        if key in wanted and key not in found and fields[6] == "":
            end = starts[i + 1] - 1 if i + 1 < len(starts) else ""
            found[key] = pygrib.fromstring(get(url, {"Range": f"bytes={starts[i]}-{end}"}))
    missing = set(wanted) - set(found)
    if missing:
        raise RuntimeError(f"{url}: missing {sorted(missing)}")
    return found


def nbm_ground(run: dt.datetime) -> pygrib.gribmessage:
    """NBM's ground height: URMA's (the same grid), from the newest analysis found."""
    for back in range(48):
        t = run - dt.timedelta(hours=back)
        url = f"{URMA}/urma2p5.{t:%Y%m%d}/urma2p5.t{t:%H}z.2dvaranl_ndfd.grb2_wexp"
        try:
            return records(url, ["HGT:surface"])["HGT:surface"]
        except (urllib.error.HTTPError, RuntimeError):
            continue
    raise SystemExit("No URMA analysis for NBM's ground height in the last 48 hours")


class Sampler:
    """Where each tile point falls on NBM's Lambert grid."""

    def __init__(self, message: pygrib.gribmessage, lons: np.ndarray, lats: np.ndarray):
        proj = pyproj.Proj(message.projparams)
        lat0, lon0 = message.latlons()
        x0, y0 = proj(lon0[0, 0], lat0[0, 0])
        x, y = proj(lons, lats)
        fx = (x - x0) / message["DxInMetres"]
        fy = (y - y0) / message["DyInMetres"]
        self.ny, self.nx = message.values.shape
        self.inside = (fx >= 0) & (fy >= 0) & (fx <= self.nx - 1) & (fy <= self.ny - 1)
        self.x0 = np.clip(np.floor(fx).astype(int), 0, self.nx - 2)
        self.y0 = np.clip(np.floor(fy).astype(int), 0, self.ny - 2)
        self.tx = fx - self.x0
        self.ty = fy - self.y0

    def sample(self, values: np.ndarray) -> np.ndarray:
        if values.shape != (self.ny, self.nx):
            raise RuntimeError(f"grid {values.shape} is not NBM's {(self.ny, self.nx)}")
        v = np.ma.filled(values.astype(np.float32), np.nan)
        x0, y0, tx, ty = self.x0, self.y0, self.tx, self.ty
        top = v[y0, x0] + (v[y0, x0 + 1] - v[y0, x0]) * tx
        bottom = v[y0 + 1, x0] + (v[y0 + 1, x0 + 1] - v[y0 + 1, x0]) * tx
        out = top + (bottom - top) * ty
        return np.where(self.inside, out, np.nan)

    def nearest(self, values: np.ndarray) -> np.ndarray:
        """For a direction, which can't be averaged across north."""
        v = np.ma.filled(values.astype(np.float32), np.nan)
        out = v[self.y0 + (self.ty > 0.5), self.x0 + (self.tx > 0.5)]
        return np.where(self.inside, out, np.nan)


def tiles() -> list[tuple[str, int, int]]:
    out = []
    for south in range(SOUTH, NORTH, TILE_DEG):
        for west in range(WEST, EAST, TILE_DEG):
            name = f"{'w' if west < 0 else 'e'}{abs(west)}{'n' if south >= 0 else 's'}{abs(south)}"
            out.append((name, west, south))
    return out


def quantize(values: np.ndarray, scale: float = SCALE) -> np.ndarray:
    q = np.round(values / scale)
    return np.where(np.isfinite(q), np.clip(q, -127, 127), NO_DATA).astype(np.int8)


def build(run: dt.datetime, out: Path) -> None:
    tile_list = tiles()
    offsets = np.arange(PER_SIDE) * STEP
    lons = np.concatenate([np.tile(west + offsets, PER_SIDE) for _, west, _ in tile_list])
    lats = np.concatenate([np.repeat(south + TILE_DEG - offsets, PER_SIDE) for _, _, south in tile_list])
    per_tile = PER_SIDE * PER_SIDE

    first = records(step_url(run, HOURS[0]), WANTED)
    sampler = Sampler(first["TMP:2 m above ground"], lons, lats)
    ground = sampler.sample(nbm_ground(run).values)

    grids = []  # per hour: temperature, u, v at every point, int8
    for hour in HOURS:
        g = first if hour == HOURS[0] else records(step_url(run, hour), WANTED)
        t = sampler.sample(g["TMP:2 m above ground"].values) - 273.15
        speed = sampler.sample(g["WIND:10 m above ground"].values)
        towards = np.radians(sampler.nearest(g["WDIR:10 m above ground"].values))
        # Blowing from `towards`: u, v point the other way.
        u, v = -speed * np.sin(towards), -speed * np.cos(towards)
        grids.append((quantize(t), quantize(u), quantize(v)))
        print(f"  f{hour:03d}", flush=True)

    stamp = run.strftime("%Y%m%d%H")
    run_iso = run.strftime("%Y-%m-%dT%H:%M:%SZ")
    times = [int((run + dt.timedelta(hours=h)).replace(tzinfo=dt.timezone.utc).timestamp() * 1000) for h in HOURS]
    folder = out / SET / stamp
    folder.mkdir(parents=True, exist_ok=True)
    index_tiles = []
    total = 0
    for k, (name, west, south) in enumerate(tile_list):
        part = slice(k * per_tile, (k + 1) * per_tile)
        if not sampler.inside[part].any():
            continue
        meta = {
            "format": FORMAT,
            "model": "NOAA NBM 2.5 km, 2 m temperature and 10 m wind",
            "set": SET,
            "run": run_iso,
            "west": west,
            "east": west + TILE_DEG,
            "south": south,
            "north": south + TILE_DEG,
            "dx": STEP,
            "dy": STEP,
            "nx": PER_SIDE,
            "ny": PER_SIDE,
            "scale": SCALE,
            "noData": NO_DATA,
            "times": times,
        }
        data = b"".join(t[part].tobytes() + u[part].tobytes() + v[part].tobytes() for t, u, v in grids)
        heights = np.where(np.isfinite(ground[part]), np.round(ground[part]), 0).astype("<i2")
        path = folder / f"{name}.zip"
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            z.writestr("hrrr.json", json.dumps(meta))
            z.writestr("hrrr.bin", data)
            z.writestr("ground.bin", heights.tobytes())
        size = path.stat().st_size
        total += size
        index_tiles.append({"id": name, "west": west, "south": south, "bytes": size})
    index = {
        "format": FORMAT,
        "model": "NOAA NBM 2.5 km",
        "set": SET,
        "run": run_iso,
        "builtAt": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "folder": f"{SET}/{stamp}",
        "tileDeg": TILE_DEG,
        "times": times,
        "tiles": index_tiles,
    }
    (out / f"{SET}.json").write_text(json.dumps(index, indent=1) + "\n")
    print(f"{SET} {run_iso}: {len(index_tiles)} tiles, {total / 1e6:.1f} MB, {len(HOURS)} hours")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", help="NBM run as YYYYMMDDHH (default: the latest complete 00/06/12/18Z one)")
    parser.add_argument("--out", default="build")
    parser.add_argument("--if-newer-than", help="URL of a published index: exit 3 without building if its run is current")
    args = parser.parse_args()
    if args.run:
        run = dt.datetime.strptime(args.run, "%Y%m%d%H")
    else:
        run = latest_run(dt.datetime.now(dt.timezone.utc).replace(tzinfo=None))
    print(f"NBM run {run:%Y-%m-%d %HZ}")
    if args.if_newer_than:
        try:
            if json.loads(get(args.if_newer_than)).get("run") == run.strftime("%Y-%m-%dT%H:%M:%SZ"):
                print("Already published")
                sys.exit(3)
        except urllib.error.HTTPError:
            pass
    build(run, Path(args.out))


if __name__ == "__main__":
    main()
