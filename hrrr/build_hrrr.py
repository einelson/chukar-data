"""
Builds Chukar's HRRR tiles: NOAA HRRR (3 km, lower 48) 2 m temperature and 10 m wind, regridded to 2° lon/lat tiles
for the app's wind and temperature maps. See the app repo's docs/HRRR_PLAN.md.

Reads one run from NOAA's open data bucket on AWS (noaa-hrrr-bdp-pds), only the records needed, by byte range.

Two sets:
  hourly  the latest run, hours 0-18
  main    the latest 00/06/12/18Z run, hours 0-48

Output, in --out (default build/):
  <set>/<run>/<tile>.zip   hrrr.json (the metadata) + hrrr.bin (the grids) + ground.bin (HRRR's ground height)
  <set>.json               the index: run, times, the tiles and where they are

hrrr.bin, per forecast time in order: temperature, then u, then v, each ny rows from north to south of nx int8
values in steps of `scale` (°C or m/s; u toward the east, v toward the north, turned to true north). -128 means no
data. ground.bin: ny x nx little-endian int16 metres, rows north to south.

sky.json + sky.bin (added 2026-11-20; older app builds don't read them): per forecast time, visibility then cloud
ceiling, each ny x nx uint8 rows north to south. Visibility in steps of `visScale` metres, `visCap` meaning that far
or more; ceiling (the lowest broken or overcast layer, metres above sea level) in steps of `ceilingScale`, 255 = none
below 6 km.

A tile is 2° square with points every 0.025° (81 x 81), its edges shared with its neighbours'. Tile names are the
south-west corner: "w112n40" is 112-110°W, 40-42°N.

Usage: python build_hrrr.py --set hourly|main [--run 2026111612] [--out build] [--if-newer-than URL]
"""

from __future__ import annotations

import argparse
import datetime as dt
import io
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

BUCKET = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
FORMAT = 1
SCALE = 0.5
NO_DATA = -128
TILE_DEG = 2
STEP = 0.025
PER_SIDE = int(round(TILE_DEG / STEP)) + 1
# The GFS lower-48 box: 30 x 13 tiles.
WEST, EAST, SOUTH, NORTH = -126, -66, 24, 50
SETS = {"hourly": 18, "main": 48}
WANTED = ["TMP:2 m above ground", "UGRD:10 m above ground", "VGRD:10 m above ground"]
# Fog and low cloud (sky.bin): visibility at the surface, and the ceiling's height above sea level.
SKY = ["VIS:surface", "HGT:cloud ceiling"]
# Only short ranges matter for seeing game, and only ceilings that can hide terrain (the lower 48's highest point is
# 4,421 m); coarse steps keep the tiles small (the noisy visibility field compresses poorly).
VIS_SCALE = 200
VIS_CAP = 50  # 10 km: "clear"
CEILING_SCALE = 100
CEILING_TOP = 6000
NO_CEILING = 255
USER_AGENT = "Chukar HRRR builder (github.com/einelson/chukar-data)"


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
    return f"{BUCKET}/hrrr.{run:%Y%m%d}/conus/hrrr.t{run:%H}z.wrfsfcf{hour:02d}.grib2"


def latest_run(now: dt.datetime, which: str) -> dt.datetime:
    """The newest run whose last hour is out (HRRR writes its hours in order)."""
    hour = now.replace(minute=0, second=0, microsecond=0)
    for back in range(12):
        run = hour - dt.timedelta(hours=back)
        if which == "main" and run.hour % 6:
            continue
        try:
            get(step_url(run, SETS[which]) + ".idx")
            return run
        except urllib.error.HTTPError:
            continue
    raise SystemExit(f"No complete HRRR {which} run in the last 12 hours")


def records(run: dt.datetime, hour: int, wanted: list[str]) -> dict[str, pygrib.gribmessage]:
    url = step_url(run, hour)
    lines = get(url + ".idx").decode().splitlines()
    starts = [int(line.split(":")[1]) for line in lines]
    found = {}
    for i, line in enumerate(lines):
        fields = line.split(":")
        key = f"{fields[3]}:{fields[4]}"
        if key in wanted and key not in found:
            end = starts[i + 1] - 1 if i + 1 < len(starts) else ""
            found[key] = pygrib.fromstring(get(url, {"Range": f"bytes={starts[i]}-{end}"}))
    missing = set(wanted) - set(found)
    if missing:
        raise RuntimeError(f"f{hour:02d}: missing {sorted(missing)}")
    return found


class Sampler:
    """Where each tile point falls on HRRR's Lambert grid, and the turn from grid to true north there."""

    def __init__(self, message: pygrib.gribmessage, lons: np.ndarray, lats: np.ndarray):
        params = message.projparams
        proj = pyproj.Proj(params)
        lat0, lon0 = message.latlons()
        x0, y0 = proj(lon0[0, 0], lat0[0, 0])
        x, y = proj(lons, lats)
        self.fx = (x - x0) / message["DxInMetres"]
        self.fy = (y - y0) / message["DyInMetres"]
        self.ny, self.nx = message.values.shape
        self.inside = (self.fx >= 0) & (self.fy >= 0) & (self.fx <= self.nx - 1) & (self.fy <= self.ny - 1)
        angle = np.sin(np.radians(params["lat_1"])) * np.radians(((lons - params["lon_0"] + 180) % 360) - 180)
        self.cos, self.sin = np.cos(angle), np.sin(angle)
        self.x0 = np.clip(np.floor(self.fx).astype(int), 0, self.nx - 2)
        self.y0 = np.clip(np.floor(self.fy).astype(int), 0, self.ny - 2)
        self.tx = self.fx - self.x0
        self.ty = self.fy - self.y0

    def sample(self, values: np.ndarray) -> np.ndarray:
        v = np.ma.filled(values.astype(np.float32), np.nan)
        x0, y0, tx, ty = self.x0, self.y0, self.tx, self.ty
        top = v[y0, x0] + (v[y0, x0 + 1] - v[y0, x0]) * tx
        bottom = v[y0 + 1, x0] + (v[y0 + 1, x0 + 1] - v[y0 + 1, x0]) * tx
        out = top + (bottom - top) * ty
        return np.where(self.inside, out, np.nan)

    def wind(self, u: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return self.cos * u + self.sin * v, -self.sin * u + self.cos * v


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


def visibility_bytes(metres: np.ndarray) -> np.ndarray:
    q = np.round(metres / VIS_SCALE)
    return np.where(np.isfinite(q), np.clip(q, 0, VIS_CAP), VIS_CAP).astype(np.uint8)


def ceiling_bytes(metres: np.ndarray) -> np.ndarray:
    q = np.round(metres / CEILING_SCALE)
    ok = np.isfinite(metres) & (metres <= CEILING_TOP)
    return np.where(ok, np.clip(q, 0, 254), NO_CEILING).astype(np.uint8)


def build(which: str, run: dt.datetime, out: Path) -> None:
    hours = list(range(SETS[which] + 1))
    tile_list = tiles()
    # Every tile's points, tile after tile, each 81 x 81 rows north to south.
    offsets = np.arange(PER_SIDE) * STEP
    lons = np.concatenate([np.tile(west + offsets, PER_SIDE) for _, west, _ in tile_list])
    lats = np.concatenate([np.repeat(south + TILE_DEG - offsets, PER_SIDE) for _, _, south in tile_list])
    per_tile = PER_SIDE * PER_SIDE

    first = records(run, 0, WANTED + SKY + ["HGT:surface"])
    sampler = Sampler(first["TMP:2 m above ground"], lons, lats)
    ground = sampler.sample(first["HGT:surface"].values)

    grids = []  # per hour: temperature, u, v at every point, int8
    sky = []  # per hour: visibility, ceiling at every point, uint8
    for hour in hours:
        g = first if hour == 0 else records(run, hour, WANTED + SKY)
        t = sampler.sample(g["TMP:2 m above ground"].values) - 273.15
        u, v = sampler.wind(sampler.sample(g["UGRD:10 m above ground"].values), sampler.sample(g["VGRD:10 m above ground"].values))
        grids.append((quantize(t), quantize(u), quantize(v)))
        vis = sampler.sample(g["VIS:surface"].values)
        ceiling = sampler.sample(g["HGT:cloud ceiling"].values)
        sky.append((visibility_bytes(vis), ceiling_bytes(ceiling)))
        print(f"  f{hour:02d}", flush=True)

    stamp = run.strftime("%Y%m%d%H")
    run_iso = run.strftime("%Y-%m-%dT%H:%M:%SZ")
    times = [int((run + dt.timedelta(hours=h)).replace(tzinfo=dt.timezone.utc).timestamp() * 1000) for h in hours]
    folder = out / which / stamp
    folder.mkdir(parents=True, exist_ok=True)
    index_tiles = []
    total = 0
    for k, (name, west, south) in enumerate(tile_list):
        part = slice(k * per_tile, (k + 1) * per_tile)
        if not sampler.inside[part].any():
            continue
        meta = {
            "format": FORMAT,
            "model": "NOAA HRRR 3 km, 2 m temperature and 10 m wind",
            "set": which,
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
            z.writestr(
                "sky.json",
                json.dumps(
                    {
                        "format": 1,
                        "fields": ["visibility", "ceiling"],
                        "visScale": VIS_SCALE,
                        "visCap": VIS_CAP,
                        "ceilingScale": CEILING_SCALE,
                        "noCeiling": NO_CEILING,
                    }
                ),
            )
            z.writestr("sky.bin", b"".join(vis[part].tobytes() + ceil[part].tobytes() for vis, ceil in sky))
        size = path.stat().st_size
        total += size
        index_tiles.append({"id": name, "west": west, "south": south, "bytes": size})
    index = {
        "format": FORMAT,
        "model": "NOAA HRRR 3 km",
        "set": which,
        "run": run_iso,
        "builtAt": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "folder": f"{which}/{stamp}",
        "tileDeg": TILE_DEG,
        "times": times,
        "tiles": index_tiles,
    }
    (out / f"{which}.json").write_text(json.dumps(index, indent=1) + "\n")
    print(f"{which} {run_iso}: {len(index_tiles)} tiles, {total / 1e6:.1f} MB, {len(hours)} hours")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--set", choices=sorted(SETS), required=True)
    parser.add_argument("--run", help="HRRR run as YYYYMMDDHH (default: the latest complete one)")
    parser.add_argument("--out", default="build")
    parser.add_argument("--if-newer-than", help="URL of a published index: exit 3 without building if its run is current")
    args = parser.parse_args()
    if args.run:
        run = dt.datetime.strptime(args.run, "%Y%m%d%H")
    else:
        run = latest_run(dt.datetime.now(dt.timezone.utc).replace(tzinfo=None), args.set)
    print(f"HRRR {args.set} run {run:%Y-%m-%d %HZ}")
    if args.if_newer_than:
        try:
            if json.loads(get(args.if_newer_than)).get("run") == run.strftime("%Y-%m-%dT%H:%M:%SZ"):
                print("Already published")
                sys.exit(3)
        except urllib.error.HTTPError:
            pass
    build(args.set, run, Path(args.out))


if __name__ == "__main__":
    main()
