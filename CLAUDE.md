# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-purpose R pipeline that converts FWRI side-scan sonar microgrid + digitized habitat polygons (delivered as a `.gdb` geodatabase) into per-cell proportional habitat-coverage ASCII rasters for the West Florida Shelf Ecospace model. There is no build system, package, or test suite — it's a driver script plus a function file, run interactively in R. See `README.md` for the full data-flow narrative, function reference, and a worked example.

## Running it

```r
setwd("path/to/GFISHER")   # repo root; all other paths resolve relative to it
# edit the USER INPUTS block at the top of `process GFISHER data.R`:
#   file.gdb <- "<absolute path to your GFISHER_EAST_Universe_2026.gdb>"
#   res      <- 5     # or 15
source("process GFISHER data.R")
```

`process GFISHER data.R` is the only entry point: it sets the two user inputs, sources `R/GFISHER functions.R`, picks the matching depth raster from `data/bathymetry/`, then calls `fn.make_GFISHER_habitat_maps()` followed by `fn.plot_GFISHER_habitats()`. The geodatabase read + centroid intersection + IDW fill takes minutes; once the `.asc` outputs exist they can be reused directly (re-run only `fn.plot_GFISHER_habitats()` to redraw figures).

Requires R 4.x with CRAN packages: `raster`, `terra`, `sf`, `sp`, `lwgeom`, `gstat`, `colorRamps`, `maps`. `raster` is loaded by the driver; the rest load inside `R/GFISHER functions.R`.

## Architecture notes that aren't obvious from one file

- **Two functions matter, two are dead.** `R/GFISHER functions.R` defines four functions but the driver only uses `fn.make_GFISHER_habitat_maps` and `fn.plot_GFISHER_habitats`. `fn.make_GFISHER_habitat_maps_old` is a superseded two-geodatabase variant, and `fn.make_gfisher_videodataset` is unrelated species/video-count processing that depends on globals never defined in this repo (`modgrps`, `region`, `read_excel`, `rtruncnorm`, `sizeatage`). Don't treat those two as live code paths.
- **Habitat classes** are a 2-letter code derived from `NewHabStrat`: `{A,N}` (artificial/natural) × `{L,M,H}` (low/medium/high relief). Class `AP` is dropped; `MULTISURFACE` geometries are dropped. The six classes are reordered to `AL,AM,AH,NL,NM,NH` via a hardcoded index (`newhabs[c(2,3,1,5,6,4)]`) — if the source class set changes, that index is what breaks.
- **The two habitat families fill empty cells differently.** Natural (`NL/NM/NH`) cells outside the mapped footprint are gap-filled by local IDW interpolation (`idw_fill_raster_longlat`, projected through an Albers equal-area CRS, `idp=4, nmax=8`, only for depth ≤ 500 m, clamped to [0,1]). Artificial (`AL/AM/AH`) cells are *not* extrapolated — unmapped water is set to 0, because artificial structure is point-like, not a gradient. This asymmetry is intentional; preserve it.
- **Depth cutoffs are baked in:** habitat proportions are forced to 0 where `depth > 200` m, and IDW only predicts where `depth ≤ 500` m. Land/no-depth cells are `NA` throughout.
- **`res` is overloaded.** The user-input global `res` (5 or 15) shares a name with `raster::res()`. Inside the function `res.min <- round(res(depth)[1]*60,0)` calls the *function*; the global `res` is only read by the driver to pick the depth file. Renaming one without the other will silently break.
- **`terra` is loaded mid-function** (`library('terra')` near the end of `fn.make_GFISHER_habitat_maps`), which masks several `raster` generics. Keep `terra` loading late; loading it earlier can change which `plot`/`crs`/`extent` methods dispatch.
- **The driver opens a Windows graphics device** (`windows(record=T)`) and the functions `plot()` intermediate rasters as side effects. This is Windows-only and assumes an interactive session — headless/non-Windows runs need those calls neutralized.

## Inputs, outputs, and what's gitignored

- **User-supplied input:** the `.gdb` geodatabase is large, machine-local, and never committed (`.gitignore` excludes `data/**/*.gdb/`). Expose its location only via the `file.gdb` USER INPUT — do not hardcode or commit geodatabase paths. The `data/` folder also holds several local-only `.gdb` directories and spreadsheets that are not repo inputs.
- **Ships with the repo:** depth + exclusion-mask ASCII rasters in `data/bathymetry/` (at 5min `66x78` and 15min `22x26`). Filenames embed resolution and dimensions; the driver greps for `depth <res>min` to find the grid.
- **Generated outputs** (gitignored, written under `output/`): habitat-proportion rasters and figures go to `output/maps/<res>min/` — `GFISHER_<class>_prop_<res>min_<rows>x<cols>.asc` for the six classes, `GFISHER_microgrid_<res>min_...asc` (mapping footprint area per cell), an optional `hab area vs grid area.csv` exception report, and the `Habitat prop area <res>min.png` / `mapping footprint <res>min.png` figures. MaxN heatmap rasters go to `output/maps/GFISHER/<res>min/maxn/`.
- `output/maps/GFISHER/`, `output/maps/<res>min/`, `output/affinity_selratio/`, `data/**/*.gdb/`, and `hoard/` (author scratch) are gitignored.
