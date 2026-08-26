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

`process GFISHER data.R` is the only entry point. It sets the user inputs, sources the map-stage files, picks the matching depth raster from `data/bathymetry/`, and runs four numbered stages. Stage 1 (basemaps) takes minutes — the geodatabase read dominates; once the `.asc` outputs exist, re-run only `fn.plot_habitat_basemaps()` to redraw figures.

Requires R 4.x with CRAN packages: `raster`, `terra`, `sf`, `sp`, `lwgeom`, `gstat`, `colorRamps`, `maps`, `FNN`, `xlsx`, `reshape2`, `truncnorm`, and `mgcv` only for the non-default `target='smooth'`.

## Architecture notes that aren't obvious from one file

- **The code is split by pipeline stage, one file per stage.** `R/habitat_basemaps.R` (stage 1), `R/video_dataset.R` (2), `R/maxn_maps.R` (3), and the three affinity modules (4). The former monolithic `R/GFISHER functions.R` no longer exists — it was split into those files, and its cell-area habitat maps were retired to `R/legacy/habitat_maps_cellarea.R`. Don't add new map code to the affinity modules or vice versa.
- **`R/legacy/` is superseded, not dead.** `fn.make_GFISHER_habitat_maps` still runs and is kept so older Ecospace runs stay reproducible; the driver shows the commented call. `fn.make_GFISHER_habitat_maps_old` and `fn.plot_GFISHER_habitats` are genuinely dead — the latter reads the retired `GFISHER_<CLS>_prop_*` filenames and expects a microgrid footprint raster, so it breaks on the new basemaps.
- **Habitat classes** are a 2-letter code derived from `NewHabStrat`: `{A,N}` (artificial/natural) × `{L,M,H}` (low/medium/high relief). `REEF.CLASSES` in `R/habitat_basemaps.R` fixes the order as `AL,AM,AH,NL,NM,NH`; the legacy module instead reorders via a hardcoded index (`newhabs[c(2,3,1,5,6,4)]`).
- **The basemap layers sum to exactly 1 per water cell** — six reef classes + `RCK` + `UNC` + `SGR`. Anything added has to enter the normalization in `fn.make_habitat_basemaps`, or the invariant breaks silently. The QC block prints the min/max row sum for exactly this reason.
- **Reef density divides by SCANNED area, not cell area.** GFISHER has scanned 2.95% of the domain. Unmapped ground is filled by shrinking toward a region × depth-bin stratum mean (`target='stratum'`). `target='smooth'` and `target='idw'` exist but were both evaluated and rejected — the file header records why, including that the GAM was 48× miscalibrated. Don't switch the default back without re-reading it.
- **dbSeabed is two variables, not one composition.** `GVL+SND+MUD` is a closed grain-size triangle summing to ~99; `RCK` is a separate areal fraction of rock outcrop. Read the RAW grids (`Gmf_<CLS>/gmf_<CLS>_val.asc`, 1.2 arc-min) and treat `< 0` as NA — the processed `gmf_*_prop_*.asc` layers in the Ecospace maps tree renormalize all four together and convert the `-99` NODATA flag to 0.
- **Depth cutoffs:** `depth.max.reef` (default 300 m, the deepest cell holding any observed reef) forces reef to 0 below it. The legacy module instead used 200 m, and IDW-predicted only where depth ≤ 500 m. Land/no-depth cells are `NA` throughout.
- **`res` is overloaded.** The user-input global `res` (5 or 15) shares a name with `raster::res()`. Inside the function `res.min <- round(res(depth)[1]*60,0)` calls the *function*; the global `res` is only read by the driver to pick the depth file. Renaming one without the other will silently break.
- **`terra` is loaded mid-function** in the legacy `fn.make_GFISHER_habitat_maps`, which masks several `raster` generics. `R/habitat_basemaps.R` avoids this by namespacing every call (`raster::`), and accepts a `SpatRaster` template by coercing it up front.
- **The driver opens a Windows graphics device** (`windows(record=T)`) and the functions `plot()` intermediate rasters as side effects. This is Windows-only and assumes an interactive session — headless/non-Windows runs need those calls neutralized.

## Inputs, outputs, and what's gitignored

- **User-supplied input:** the `.gdb` geodatabase is large, machine-local, and never committed (`.gitignore` excludes `data/**/*.gdb/`). Expose its location only via the `file.gdb` USER INPUT — do not hardcode or commit geodatabase paths. The `data/` folder also holds several local-only `.gdb` directories and spreadsheets that are not repo inputs.
- **Ships with the repo:** depth + exclusion-mask ASCII rasters in `data/bathymetry/` (at 5min `66x78` and 15min `22x26`). Filenames embed resolution and dimensions; the driver greps for `depth <res>min` to find the grid.
- **Generated outputs** (written under `output/`):
  - `output/basemaps/<res>min/` — the live habitat layers: `habitat_<CODE>_<res>min_<rows>x<cols>.asc` (nine of them), `habitat_basemap_QC_<res>min.csv`, and the `habitat_basemaps_<res>min.png` panel figure.
  - `output/maps/GFISHER/<res>min/maxn/<scheme>/` — MaxN heatmap rasters, tagged by grouping scheme so `ecospace` and `mice` coexist.
  - `output/affinity_*_<scheme>/` — the three affinity routes.
  - `output/maps/<res>min/` — only if you deliberately re-run the retired legacy maps.
- `output/maps/GFISHER/`, `output/maps/<res>min/`, `output/affinity_selratio*/`, `data/**/*.gdb/`, `data/April2026/`, and `hoard/` (author scratch) are gitignored.
- **`output/basemaps/5min/` is deliberately NOT gitignored** — it holds the sum-to-1 layers that feed Ecospace and is treated as a versioned deliverable, so a rebuild shows up in `git diff`. Don't add a blanket `output/` rule that would swallow it.
