"""
Builds Chukar's wind and temperature map data from NOAA's GFS forecast: 10 m wind and 2 m temperature on the 0.25°
grid, cut to each region and packed small enough for a phone to keep for offline use.

Run by .github/workflows/wind.yml every six hours. Reads the GFS run from NOAA's open data bucket on AWS
(noaa-gfs-bdp-pds; public domain), fetching only the three records needed of each forecast hour by byte range.

Output, in --out (default build/):
  wind-<region>.zip   wind.json (the metadata below) + wind.bin (the grids)
  wind.json           the index the app checks: the run, and each region's box, file and size
  temp-<region>.zip   temp.json + temp.bin, the same way
  temp.json           the temperature index

wind.bin, per forecast time in order: the u grid then the v grid, each ny rows from north to south of nx int8
values, in units of `scale` m/s (u: toward the east, v: toward the north). -128 means no data.

temp.bin, per forecast time in order: one grid, laid out the same, in units of `scale` °C.

Usage: python build_wind.py [--out build] [--run 2026101706] [--if-newer-than https://.../wind.json]
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

BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
FORMAT = 1
SCALE = 0.5  # m/s or °C per step: int8 covers ±63.5 m/s (142 mph) and ±63.5 °C (-82..146 °F)
NO_DATA = -128

# Forecast hours kept: hourly for a day, 3-hourly to 3 days, 6-hourly to 5 days (49 grids).
HOURS = list(range(0, 25)) + list(range(27, 73, 3)) + list(range(78, 121, 6))

# Boxes on the 0.25° grid, in degrees (west and east as -180..180).
REGIONS = {
    "conus": {"west": -125.0, "east": -66.0, "south": 24.0, "north": 50.0},
    "alaska": {"west": -170.0, "east": -129.0, "south": 51.0, "north": 72.0},
    "hawaii": {"west": -161.0, "east": -154.0, "south": 18.0, "north": 23.0},
}

WANTED = {("UGRD", "10 m above ground"), ("VGRD", "10 m above ground"), ("TMP", "2 m above ground")}

USER_AGENT = "Chukar wind builder (github.com/einelson/chukar-data)"


def get(url: str, headers: dict | None = None, tries: int = 4) -> bytes:
    for attempt in range(tries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise
            if attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
        time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def step_url(run: dt.datetime, hour: int) -> str:
    return f"{BUCKET}/gfs.{run:%Y%m%d}/{run:%H}/atmos/gfs.t{run:%H}z.pgrb2.0p25.f{hour:03d}"


def run_is_complete(run: dt.datetime) -> bool:
    """A run is complete once the last hour we need has its index (GFS writes hours in order)."""
    try:
        get(step_url(run, HOURS[-1]) + ".idx")
        return True
    except urllib.error.HTTPError:
        return False


def latest_run(now: dt.datetime) -> dt.datetime:
    cycle = now.replace(minute=0, second=0, microsecond=0, hour=now.hour - now.hour % 6)
    for back in range(8):
        run = cycle - dt.timedelta(hours=6 * back)
        if run_is_complete(run):
            return run
    raise SystemExit("No complete GFS run in the last two days")


def records(run: dt.datetime, hour: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The 10 m u and v (m/s) and 2 m temperature (°C) grids (721 x 1440, north to south, 0..359.75° east) of one forecast hour."""
    url = step_url(run, hour)
    lines = get(url + ".idx").decode().splitlines()
    starts = [int(line.split(":")[1]) for line in lines]
    grids = {}
    for i, line in enumerate(lines):
        fields = line.split(":")
        if (fields[3], fields[4]) in WANTED:
            end = starts[i + 1] - 1 if i + 1 < len(starts) else ""
            data = get(url, {"Range": f"bytes={starts[i]}-{end}"})
            message = pygrib.fromstring(data)
            values = np.ma.filled(message.values.astype(np.float32), np.nan)
            lats = message.latlons()[0][:, 0]
            if lats[0] < lats[-1]:  # always store north to south
                values = values[::-1]
            grids[fields[3]] = values
    if set(grids) != {variable for variable, _ in WANTED}:
        raise RuntimeError(f"f{hour:03d}: missing {sorted({v for v, _ in WANTED} - set(grids))} in the index")
    return grids["UGRD"], grids["VGRD"], grids["TMP"] - 273.15


def crop(grid: np.ndarray, box: dict) -> np.ndarray:
    """The box's part of a global grid: rows 0 = 90°N, columns 0 = 0°E, 0.25° apart."""
    row0 = int(round((90 - box["north"]) / 0.25))
    row1 = int(round((90 - box["south"]) / 0.25)) + 1
    col0 = int(round((box["west"] % 360) / 0.25))
    nx = int(round((box["east"] - box["west"]) / 0.25)) + 1
    cols = [(col0 + i) % 1440 for i in range(nx)]
    return grid[row0:row1][:, cols]


def quantize(grid: np.ndarray) -> np.ndarray:
    q = np.round(grid / SCALE)
    q = np.where(np.isfinite(q), np.clip(q, -127, 127), NO_DATA)
    return q.astype(np.int8)


def write_product(
    out: Path, name: str, model: str, run: dt.datetime, times: list[int], parts: dict[str, list[bytes]], grids_per_time: int
) -> None:
    """Writes <name>-<region>.zip for each region, then the <name>.json index that names them."""
    built = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    index = {
        "format": FORMAT,
        "model": model,
        "run": run.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "builtAt": built.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "regions": [],
    }
    for region, box in REGIONS.items():
        nx = int(round((box["east"] - box["west"]) / 0.25)) + 1
        ny = int(round((box["north"] - box["south"]) / 0.25)) + 1
        meta = {
            "format": FORMAT,
            "region": region,
            "run": index["run"],
            "builtAt": index["builtAt"],
            **box,
            "dx": 0.25,
            "dy": 0.25,
            "nx": nx,
            "ny": ny,
            "scale": SCALE,
            "noData": NO_DATA,
            "times": times,
        }
        data = b"".join(parts[region])
        assert len(data) == len(times) * grids_per_time * nx * ny
        file = f"{name}-{region}.zip"
        with zipfile.ZipFile(out / file, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            z.writestr(f"{name}.json", json.dumps(meta))
            z.writestr(f"{name}.bin", data)
        index["regions"].append({"id": region, **box, "file": file, "bytes": (out / file).stat().st_size})
        print(f"{file}: {nx}x{ny}, {len(times)} times, {(out / file).stat().st_size / 1e6:.2f} MB")
    (out / f"{name}.json").write_text(json.dumps(index, indent=2) + "\n")


def build(run: dt.datetime, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    wind = {region: [] for region in REGIONS}
    temp = {region: [] for region in REGIONS}
    times = []
    for hour in HOURS:
        u, v, t = records(run, hour)
        times.append(int((run + dt.timedelta(hours=hour)).timestamp() * 1000))
        for region, box in REGIONS.items():
            wind[region].append(quantize(crop(u, box)).tobytes() + quantize(crop(v, box)).tobytes())
            temp[region].append(quantize(crop(t, box)).tobytes())
        print(f"  f{hour:03d}", flush=True)
    write_product(out, "wind", "NOAA GFS 0.25°, 10 m wind", run, times, wind, 2)
    write_product(out, "temp", "NOAA GFS 0.25°, 2 m temperature", run, times, temp, 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="build")
    parser.add_argument("--run", help="GFS run as YYYYMMDDHH (default: the latest complete one)")
    parser.add_argument("--if-newer-than", help="URL of a published index: exit 3 without building if its run is current")
    args = parser.parse_args()

    if args.run:
        run = dt.datetime.strptime(args.run, "%Y%m%d%H").replace(tzinfo=dt.timezone.utc)
    else:
        run = latest_run(dt.datetime.now(dt.timezone.utc))
    print(f"GFS run {run:%Y-%m-%d %HZ}")

    if args.if_newer_than:
        try:
            published = json.loads(get(args.if_newer_than))
            if published.get("run") == run.strftime("%Y-%m-%dT%H:%M:%SZ"):
                print("Already published")
                sys.exit(3)
        except urllib.error.HTTPError:
            pass  # nothing published yet

    build(run, Path(args.out))


if __name__ == "__main__":
    main()
