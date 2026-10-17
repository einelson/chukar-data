# Chukar data

Prebuilt offline data packs for the Chukar maps app. The app downloads them from the rolling
[`data` release](https://github.com/einelson/chukar-data/releases/tag/data). The packs are
built and published by Chukar's data pipeline, and the release's files are replaced in place, so their links never
change. A release holds at most 1,000 files, so newer packs live in the [`data-2` release](https://github.com/einelson/chukar-data/releases/tag/data-2)
(the manifest in `data` lists them with `release: "data-2"`).

The one piece of code here is the wind map's: `wind/build_wind.py`, run every six hours by the `Wind` workflow, turns
each new NOAA GFS forecast run into `wind.json` (the index: run time and each region's file) and `wind-conus.zip`,
`wind-alaska.zip` and `wind-hawaii.zip` (10 m wind on the 0.25° grid out 5 days) in `data-2`.

`data-2` also holds `crop-fields-<st>.pmtiles`: every farm field in each of the lower 48 with its crop each year
2018–2025, as PMTiles vector tiles the app streams by HTTP range or downloads whole.

- `manifest.json` lists every pack: its file, size, version (build time), and the source of each state's data.
- `<state>-<layer>[-<part>].zip`: one GeoJSON file per map cell for public land, forest roads (MVUM), USFS trails,
  BLM roads & trails, and OpenStreetMap roads, trails and points of interest.
- Per-state packs, one GeoJSON file each: `hunt-<st>.zip` (hunting units), `walk-in-<st>.zip` (walk-in hunting
  access), `fire-history-<st>.zip` (past fire perimeters), `power-lines-<st>.zip`, `cell-coverage-<st>.zip`,
  `big-game-<st>.zip` (migrations), `meshcore-<st>.zip` (mesh-radio nodes), plus bird-migration data for the whole US.

## Sources and licenses

Most of the data is public domain, from US federal agencies:

- Public land: USGS PAD-US and BLM Surface Management Agency
- Forest roads and trails: USDA Forest Service MVUM and National Forest System trails
- BLM roads & trails: BLM Ground Transportation Linear Features
- Wind: NOAA GFS forecast, from NOAA's open data on AWS
- Crop fields: USDA NASS Crop Sequence Boundaries 2018–2025
- Past fires: National Interagency Fire Center (Interagency Fire Perimeter History, WFIGS)
- Power lines: HIFLD
- Cell coverage: FCC Broadband Data Collection; 850 MHz license records via HIFLD
- Big-game migrations: USGS, *Ungulate Migrations of the Western United States*, volumes 1–6
- State outlines: US Census Bureau

The rest is shared under its own terms:

- **Hunting units** and **walk-in hunting access**: each state wildlife agency's own published GIS data, named for each
  state in `manifest.json`. Best representation only, with no warranty. Hunting units are simplified to about 30 m,
  except Oregon's: ODFW's terms forbid altering its boundaries.
- **Roads, trails and points of interest**: © OpenStreetMap contributors, cut from Geofabrik extracts. Available under
  the [Open Database License (ODbL)](https://opendatacommons.org/licenses/odbl/).
- **Cell towers** in the cell-coverage packs: [OpenCelliD](https://opencellid.org), shared under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
- **Bird migration**: eBird Observation Dataset, Cornell Lab of Ornithology, via GBIF.org
  ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).
- **MeshCore nodes**: the MeshCore project's public node directory (map.meshcore.io).

Map data may be out of date or wrong. Public-land boundaries are for planning only: check land status and access on
the ground.
