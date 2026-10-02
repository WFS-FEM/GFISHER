# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An R pipeline that turns the FWRI GFISHER side-scan habitat mapping and the 3LABS video survey into West Florida Shelf Ecospace inputs: sum-to-1 habitat basemaps, per-group MaxN heatmaps, and habitat affinities. There is no build system, package, or test suite. One driver script plus one file per pipeline stage, run interactively or with `Rscript`. See `README.md` for inputs, configuration, outputs, and caveats; see `docs/issue2-review-gfisher-repo-plan.md` for the October 2026 review findings and evidence.

## Running it

```r
# open GFISHER.Rproj (or setwd() to the repo root), then
source("process GFISHER data.R")
```

`process GFISHER data.R` is the only entry point. Order of events: root-anchor check (`GFISHER.Rproj` must be in the working directory), `R/_setup.R` (package check, plot device, input manifest), source the stage files, SETTINGS block of repo-relative defaults, `source('config.local.R')` if present, resolve derived paths, download dbSEABED grids into the default folder if absent, `fn.check_inputs()` (stops before any slow work if a required input is missing), then stages 1, 2, 3, 4a, 4b, 4c. A full run is about 16 minutes; stage 1's geodatabase read dominates.

Requires R 4.x with `sf`, `sp`, `raster`, `FNN`, `colorRamps`, `reshape2`, `truncnorm`, `readxl`. Optional: `mgcv` (only for `target='smooth'`), and `terra`, `gstat`, `lwgeom`, `maps` (only for `R/legacy/`). `fn.check_packages()` prints the install line.

## Conventions that matter

- **Machine-specific paths never go in tracked files.** Defaults in the driver are repo-relative; overrides go in `config.local.R` (gitignored), documented key by key in `config.local.example.R`. Three people must be able to run the same tracked code: the author (Dave, whose layout the defaults match), Holden, and a fresh GitHub clone. If you add an input, add it to `fn.data_manifest()` in `R/_setup.R`, to `config.local.example.R`, and to README "Getting the data" with its size and source.
- **Code is split by stage, one file each:** `R/habitat_basemaps.R` (1), `R/video_dataset.R` (2), `R/maxn_maps.R` (3), `R/selection_ratio_affinities.R` (4a), `R/site_level_affinities.R` (4b), `R/substrate_affinities.R` (4c). Don't put map code in the affinity modules or vice versa. `R/_setup.R` is base R only and is sourced before anything else.
- **Every stage file namespaces its calls** (`raster::`, `sf::`), coerces a `SpatRaster` template to a `RasterLayer` on entry, and nothing needs `terra` attached. Don't add `library('terra')` to the driver; it is only used by `R/legacy/habitat_maps_cellarea.R`.
- **Nothing Windows-only or interactive-only runs unguarded.** `fn.plot_device()` opens a recording window only when `interactive()` on Windows; `plot(depth)` in the driver is behind `if(interactive())`. Figures go to files.
- **Stage 2 is seeded.** `fn.make_gfisher_videodataset(..., seed=1)` seeds the stanza length draws on entry. Changing or removing the seed changes the 11 gag / red grouper stanza maps and everything downstream; the single-stanza groups, totals, and basemaps do not depend on it. Issue #5 tracks how much that matters for the affinities.
- **Commits are split by kind** (code / docs / regenerated outputs / housekeeping), subject ends with `(issue #N)`, body says why and what evidence was checked.

## Architecture notes that aren't obvious from one file

- **`R/legacy/` is superseded, not dead.** `fn.make_GFISHER_habitat_maps` still runs (commented call in the driver) so older Ecospace runs stay reproducible; `fn.make_GFISHER_habitat_maps_old` and `fn.plot_GFISHER_habitats` are dead. The legacy module loads `terra` mid-function and reorders classes with a hardcoded index `newhabs[c(2,3,1,5,6,4)]`.
- **Habitat classes** are a 2-letter code from `NewHabStrat`: `{A,N}` x `{L,M,H}`. `REEF.CLASSES` in `R/habitat_basemaps.R` fixes the order `AL,AM,AH,NL,NM,NH`. Class `AP` is dropped.
- **The basemap layers sum to exactly 1 per water cell** (six reef + `RCK` + `UNC` + `SGR`). Anything added must enter the normalisation in `fn.make_habitat_basemaps` or the invariant breaks silently; the QC block prints min/max row sums for this reason.
- **Reef density divides by SCANNED area, not cell area**, then shrinks toward a region x depth-bin mean (`target='stratum'`). `'smooth'` and `'idw'` exist and were rejected; the file header says why. `depth.max.reef=300` forces reef to 0 deeper than that.
- **dbSEABED is two variables, not one composition.** Read the RAW grids (`Gmf_<CLS>/gmf_<CLS>_val.asc`, 1.2 arc-min) and treat `< 0` as NA; the processed 5-min `gmf_*_prop_*.asc` layers elsewhere renormalise rock and turn `-99` into 0.
- **Seagrass.** Stage 1 reads `file.seagrass` if it exists, else `SGR = 0`. The committed basemaps used the author's rasters; the EcospaceBasemap `Seagrass_Statewide` raster reproduces them in all but 45 cells.
- **`res` is overloaded.** The global `res` (5 or 15) shares a name with `raster::res()`; inside functions `res(depth)` resolves to the function because R skips non-function objects when looking up a call. Renaming one without the other breaks silently.
- **Grid template.** Effort raster, MaxN maps and basemaps must share one grid; the driver checks `raster::compareRaster(hab, depth)`. Eleven survey stations lie exactly on 5-minute latitude lines, so a template with a rounded cell size assigns them differently (plan doc R5).
- **Stage 2 output holds presence records only** (the melt drops `maxn == 0`); anything needing true zeros must reinstate them against the station list, as `R/site_level_affinities.R` does.

## Inputs, outputs, and what's gitignored

- **Cannot ship (gitignored, by request from FWRI):** `data/April2026/` with `GFISHER_EAST_Universe_2026.gdb` (275 MB) and the three `*3LABS_93to24.csv` survey files (42 MB). Location set by `dir.data`; the geodatabase is found by `fn.find_gdb()`.
- **Public, tracked:** `data/dbseabed/` (4.4 MB, four raw dbSEABED grids; provenance in its `SOURCE.md`). `fn.pull_dbseabed()` re-downloads them only if the folder is removed; the CSDMS server was down during the Oct 2026 review, which is why they ship.
- **Ships with the repo:** `data/bathymetry/` depth grids (`depth <res>min <rows>x<cols>.asc`) and `data/Master Species List.xlsx`.
- **Tracked, optional at run time:** `data/seagrass/seagrass_5min.asc` (44 KB, from EcospaceBasemap; provenance in its `SOURCE.md`). Without a seagrass raster `SGR = 0`. There is no 15-min version yet.
- **Tracked outputs (deliverables):** `output/basemaps/<res>min/` (the nine layers, QC table, panel figure), `output/maps/GFISHER/<res>min/maxn/<scheme>/`, `output/affinity_*_<scheme>/`. Stage 1 overwrites `output/basemaps/` in place so a rebuild shows in `git diff`.
- **Gitignored outputs:** `output/maps/<res>min/` (legacy maps only), `output/GFISHER_species_fg*.csv`, `Rplots.pdf`, `config.local.R`.
- Two legacy geodatabases (`East_Master_Hab_data_Dissolve_byMicro_13Sept24.gdb`, `FWRI_East_Gulf_Mapping_2023.gdb`) were untracked in October 2026 and remain only in git history; no current code reads them.
