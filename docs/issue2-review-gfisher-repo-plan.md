# GFISHER issue #2: Review GFISHER repo

Date opened: 1 Oct 2026. Branch: `2-review-gfisher-repo`. Pull request: #4 (draft).
Authors: Holden Harris, with Claude Code. Original code author: David Chagaris.

Sections 1 to 4 describe what was found. Sections 5 to 8 are the work plan. Section 9 logs
decisions as they are made. This file is updated as the work proceeds and travels with the PR.

## 1. Summary

GFISHER is Dave's R pipeline that turns FWRI side-scan habitat polygons and video-survey MaxN
counts into Ecospace habitat basemaps, per-group MaxN heatmaps, and habitat affinities. The
review has three standing goals, which apply to every WFS-FEM repo:

1. Holden can run it (Windows 11, R 4.5.1, RStudio).
2. Dave can still run it with his existing paths and habits.
3. Any new GitHub user can run it, accepting that the geodatabase and survey CSVs cannot
   live on GitHub.

The branch starts from `main` and was fast-forwarded onto Dave's unmerged
`rework-habitat-basemaps` (8 commits, Aug 2026), so this PR also brings the staged pipeline
layout into `main`. Whether the basemap stage belongs in this repo at all is tracked
separately in issue #3 and is out of scope here.

## 2. What the code does (merged branch)

`process GFISHER data.R` is the only entry point. It runs four numbered stages:

| Stage | File | What it produces |
|---|---|---|
| 1 | `R/habitat_basemaps.R` | Nine sum-to-1 habitat layers per cell (`AL..NH`, `RCK`, `UNC`, `SGR`) in `output/basemaps/<res>min/` |
| 2 | `R/video_dataset.R` | Station x model-group MaxN table from the 3LABS survey CSVs and the species list |
| 3 | `R/maxn_maps.R` | One MaxN heatmap raster per model group in `output/maps/GFISHER/<res>min/maxn/<scheme>/` |
| 4a | `R/selection_ratio_affinities.R` | Habitat affinities from use/availability selection ratios |
| 4b | `R/site_level_affinities.R` | Affinities from paired camera + habitat reads at each station |
| 4c | `R/substrate_affinities.R` | Rock/gravel/sand/mud affinities from raw dbSEABED grids |

`R/legacy/habitat_maps_cellarea.R` holds the retired cell-area habitat maps (kept so older
Ecospace runs can be reproduced). `R/estimate_habitat_affinities.R` is sourced by nothing.

## 3. Findings

Line numbers refer to the merged branch at commit `3e4dea1`.

### 3.1 Portability (blocks goals 1 and 3)

| # | Where | Problem |
|---|---|---|
| P1 | `process GFISHER data.R:26` | `dir.ecospace.maps` hardcoded to Dave's OneDrive. Required by stage 1 (seagrass layer). |
| P2 | `process GFISHER data.R:78` | `dir.dbseabed` hardcoded to `C:/dchagaris/GitHub/...`. Required by stages 1 and 4c. |
| P3 | `process GFISHER data.R:1,105` | `rm(list=ls()); rm(.SavedPlots); windows(record=T)`: wipes the user's workspace, warns when `.SavedPlots` is absent, fails off Windows or under `Rscript`. |
| P4 | `process GFISHER data.R:31-39` | Survey inputs expected in gitignored `data/April2026/`; no check that they exist, so a missing file fails deep inside stage 2 after stage 1 has run for minutes. |
| P5 | `R/video_dataset.R:14,72,90` | Species list read with `xlsx`, which needs a Java runtime. `readxl` reads the same sheets without Java. |
| P6 | `R/selection_ratio_affinities.R:414-419`, `R/estimate_habitat_affinities.R:191-196` | Standalone driver blocks with Dave's paths; a second, undocumented entry point. |
| P7 | README, CLAUDE.md | List neither `FNN`, `reshape2`, `truncnorm`, `xlsx` (README) nor the four-stage layout (README). No package check, so a fresh install fails on the first missing library. |
| P8 | repo root | No root-anchor check; the driver assumes `getwd()` is the repo root. |

### 3.2 Reproducibility (blocks goal 1, objective 2)

| # | Where | Problem |
|---|---|---|
| R1 | `R/video_dataset.R:161-169` | `rtruncnorm` and `sample.int` draw random lengths with no `set.seed`, so stage 2 and everything downstream (MaxN maps, affinities) differ on every run. The affinity modules do seed their bootstraps (`selection_ratio_affinities.R:233`, `site_level_affinities.R:258`, `substrate_affinities.R:167`). |
| R2 | `output/basemaps/` | Committed basemaps include seagrass (`SGR` mean about 0.012), which needs `seagrass_<res>min.asc` from the Ecospace maps tree. Without it the code runs with `SGR = 0` and cannot reproduce the committed layers. |
| R3 | `.gitignore` vs `git ls-files` | 30 output files and 97 geodatabase files are tracked although `.gitignore` matches them. Two legacy geodatabases (`East_Master_Hab_data_Dissolve_byMicro_13Sept24.gdb`, `FWRI_East_Gulf_Mapping_2023.gdb`, about 200 MB) are referenced by no code. `data/Video Count Data4ChagarisTake2.xlsx` (14 MB) is referenced only in a comment. |
| R4 | `output/affinity_selratio/` | Untagged output folder that no current code writes (current code writes `affinity_selratio_<scheme>/`). |

### 3.3 Bugs and fragility

| # | Where | Problem |
|---|---|---|
| B1 | `R/legacy/habitat_maps_cellarea.R` | On `main`, the driver passed a terra `SpatRaster` into this `raster`-based function, which fails (`x[is.na(depth)] <- NA`: "object of type 'S4' is not subsettable"; `rasterize`, `projectRaster`, `addLayer` all error). The rework branch retires the function, and the new stage functions coerce `SpatRaster` to `RasterLayer` on entry. Legacy module still has `library('terra')` mid-function (line 236) and non-recursive `dir.create` (lines 241, 463). |
| B2 | `R/legacy/habitat_maps_cellarea.R:22` | Self-referencing default `depth=depth` gives "promise already under evaluation" if the argument is omitted. |
| B3 | `R/video_dataset.R` | Edge cases: a species whose age falls below the first stanza gives "replacement has length zero"; a multistanza species with no length records makes `rtruncnorm` error; inconsistent `[1]` indexing on size-at-age rows; count mismatch only `message()`s instead of stopping; bare REPL calls (`names()`, `length()`) left in the function body; a lookup CSV is written into `data/`. |
| B4 | `R/maxn_maps.R` | `save.format` argument is accepted but ignored (always writes ASCII). Header says `fun='sum'`; the driver passes `mean`. |
| B5 | `R/estimate_habitat_affinities.R` | Dead file; defines a second `fn.build_effort_raster` that would silently win if sourced after `selection_ratio_affinities.R`. Sole remaining holder of hardcoded `_5min_` strings. |
| B6 | `R/habitat_basemaps.R:78-79`, legacy `:34-35` | Geodatabase layers found by `grep('Microgrid')` / `grep('Hab_Data_FINAL')` with no check for zero or multiple matches. |
| B7 | `process GFISHER data.R:7` | `library('terra')` loaded after the `raster`-based stage files. On raster 3.6.32 / terra 1.8.80 the generics are shared so masking is not the issue, but terra is only needed for `rast()` on line 47 and in the legacy module. |

### 3.4 Documentation

README describes the pre-split `R/GFISHER functions.R` layout and a two-input USER INPUTS
block. CLAUDE.md was updated on the rework branch but still says the geodatabases are never
committed and still describes the USER INPUTS block.

## 4. Evidence

- `SpatRaster` scratch file found untracked in the repo root on 1 Oct 2026 recorded the B1
  errors verbatim (`rasterize sp-class: try-error`, `projectRaster SpatRaster: try-error`,
  `addLayer: try-error`, `Error in x[is.na(depth)] <- NA : object of type 'S4' is not
  subsettable`). Deleted after recording.
- `git ls-files --cached --ignored --exclude-standard` lists the tracked-but-ignored files (R3).
- Package check on Holden's machine before this work: `truncnorm` and `FNN` missing;
  `xlsx`, `readxl`, `raster 3.6-32`, `terra 1.8-80`, `sf 1.0-21` present.
- Baseline run results: *to be added in section 4.1 after the first end-to-end run.*

### 4.1 Baseline run (merged branch, unmodified)

*Pending.* Setup on Holden's machine (1 Oct 2026): `data/April2026` is a Windows directory
junction to the OneDrive copy of Dave's folder (`New-Item -ItemType Junction`), so the files
are not duplicated and git ignores the path. The dbSEABED and seagrass inputs are taken from
the local `EcospaceBasemap` repo (section 4.2); Dave's OneDrive Ecospace maps tree is not
synced to Holden's machine and is not needed.

### 4.2 Data not in the repository

Everything the pipeline reads that is **not** under `data/` in this repo, with size and where
an external user obtains it. The two legacy geodatabases that *are* tracked
(`East_Master_Hab_data_Dissolve_byMicro_13Sept24.gdb`, 78 MB, one dissolved layer;
`FWRI_East_Gulf_Mapping_2023.gdb`, 125 MB, 2023 vintage of the two layers below) are read by
no current code and are untracked in this PR.

| Input | Size | Used by | Source | How a user gets it |
|---|---|---|---|---|
| `GFISHER_EAST_Universe_2026.gdb` (layers `East_Master_Hab_Data_FINAL_2026`, `East_Master_Microgrid_Mapped_2026`) | 275 MB | Stage 1, 4b, 4c | FWRI (Sean Keenan), GFISHER side-scan habitat mapping, April 2026 delivery | Request from FWRI or Dave; place in `dir.data` |
| `maxn3LABS_93to24.csv` | 23 MB | Stage 2 | FWRI 3LABS video survey 1993 to 2024, April 2026 delivery | Same |
| `env3LABS_93to24.csv` | 6 MB | Stage 2, 4a, 4b, 4c | Same | Same |
| `lens3LABS_93to24.csv` | 13 MB | Stage 2 | Same | Same |
| `3LABS_METADATA_93to24.xlsx` | 1.6 MB | Not read by code; documents the CSVs | Same | Same (optional) |
| dbSEABED raw grids `Gmf_{RCK,GVL,SND,MUD}/gmf_*_val.asc` (891 x 383 cells at 0.02 deg) | 8.6 MB total | Stage 1, 4c | CSDMS dbSEABED "Data for Modellers", https://csdms.colorado.edu/wiki/DBSEABED | Public download; `EcospaceBasemap` has `fn.pull_dbseabed()` that fetches the four zips. Could ship with this repo (small) |
| Seagrass raster `seagrass_<res>min.asc` on the model grid | 44 KB | Stage 1 | Derived from FWC "Seagrass Habitat in Florida" (shapefile, 329 MB), https://geodata.myfwc.com/datasets/myfwc::seagrass-habitat-in-florida ; rasterised to 5 min by `EcospaceBasemap` (`seagrass_coverage_Seagrass_Statewide_5min.asc`) | Copy from `EcospaceBasemap/output/5min/habitat/seagrass/`, or ship with this repo (small). Match to the committed SGR layer: same 212 non-zero cells, r = 0.985; exactness confirmed by the baseline run |

**To do (README, Phase 4):** add a "Getting the data" section with this table (ships / public
download / request from FWRI), the expected `data/` tree, and the `config.local.R` keys that
point at each item, so an external user can locate every input and reproduce the analysis.
Decide with Dave whether the two small public-derived inputs (dbSEABED grids, seagrass raster)
should ship in `data/` so only the FWRI files need requesting.

## 5. Fix design

Mirrors `RedTideMaps` and `EcospaceBasemap`.

### 5.1 Configuration

Repo-relative defaults in the driver, then `if(file.exists('config.local.R')) source(...)`.
`config.local.R` is gitignored; `config.local.example.R` is tracked with every override
commented out and explained. Dave's current layout becomes three uncommented lines in his
local copy, so every variable he uses today keeps its name.

| Key | Default | Purpose |
|---|---|---|
| `dir.data` | `data/April2026` | Folder with the geodatabase and the three 3LABS CSVs (cannot ship) |
| `file.gdb` | `NULL` = find `GFISHER_EAST_Universe*.gdb` in `dir.data` | The FWRI geodatabase |
| `dir.dbseabed` | `data/dbseabed` | Raw dbSEABED grids `Gmf_<CLS>/gmf_<CLS>_val.asc` |
| `file.seagrass` | `data/seagrass/seagrass_<res>min.asc` | Seagrass layer for stage 1 |
| `dir.ecospace.maps` | `NULL` | Optional: Dave's Ecospace maps tree; sets `file.seagrass` and `dir.ewemaps` the old way |
| `dir.ewemaps` | `output/maps` | Where MaxN heatmaps are written |
| `dir.bathy` | `data/bathymetry` | Depth grids (ship with the repo) |
| `res` | `5` | Grid resolution in arc-minutes (5 or 15) |
| `group.scheme` | `'mice'` | Species grouping scheme |
| `seed` | `1` | Seed for the stage 2 length draws |

### 5.2 Setup checks (`R/_setup.R`, base R only)

- `fn.check_packages()`: stops with the exact `install.packages(c(...))` line.
- `fn.plot_device()`: `if(interactive() && .Platform$OS.type=='windows') try(dev.new(record=TRUE))`
  replaces `windows(record=T)`.
- `fn.data_manifest()` / `fn.check_inputs()`: one row per input (bathymetry, species list,
  seagrass, dbSEABED, geodatabase, survey CSVs) with where it comes from; prints a status
  table and stops before stage 1 if a required input is missing.
- Root anchor: `if(!file.exists('GFISHER.Rproj')) stop(...)`.

### 5.3 Code changes

- Driver: remove workspace wipe and `windows()`; drop `library('terra')` and read depth with
  `raster::raster()`; guard `plot()` with `interactive()`; source config; run checks.
- `R/video_dataset.R`: `readxl::read_excel` instead of `xlsx::read.xlsx`; new `seed=1`
  argument, `set.seed(seed)` on entry. Existing callers remain valid.
- Delete `R/estimate_habitat_affinities.R` and the standalone block at the foot of
  `R/selection_ratio_affinities.R`.
- Legacy module: `recursive=TRUE` on `dir.create`; comment that it is the only terra user.
- B3, B4, B6 fixed individually with evidence from the baseline run.

### 5.4 Housekeeping

- `git rm --cached` the two legacy geodatabases and the unused xlsx (stay on disk and in
  history; reversible). No history rewrite in this PR.
- Remove stale `output/affinity_selratio/` (untagged) from tracking.
- Rewrite `.gitignore` so that what is tracked is exactly what is not ignored.

## 6. Implementation steps

1. Commit this plan doc. (docs)
2. Copy inputs in; write a local `config.local.R`; run the merged driver as-is; record
   failures and timings in 4.1; baseline MD5s of all tracked `.asc` outputs. (no commit)
3. Portability commit: `GFISHER.Rproj` anchor, `R/_setup.R`, `config.local.example.R`,
   driver edits, `.gitignore` adds `config.local.R`. (code)
4. Reproducibility commit: `readxl` swap, `seed`, drop terra from driver. (code)
5. Dead-code commit: remove `estimate_habitat_affinities.R` and the standalone block. (code)
6. Bug-fix commits for B3, B4, B6 as confirmed. (code)
7. Docs commit: README rewrite, CLAUDE.md sync. (docs)
8. Regenerate outputs with the seeded pipeline; commit. Basemaps are expected unchanged;
   MaxN maps and affinities change once (random to seeded) and are then stable. (outputs)
9. Housekeeping commit: untrack legacy geodatabases, unused xlsx, stale output folder;
   `.gitignore` rewrite. (housekeeping)
10. Fresh-clone verification (section 7); update the PR checklist; mark ready; request
    Dave's review.

Constraints: no new dependencies beyond `readxl` (already installed); function signatures
stay backward compatible; `tools::md5sum()` for comparisons; Windows 11 / R 4.5.1.

## 7. Acceptance criteria

- [ ] `Rscript "process GFISHER data.R"` from a fresh clone with only `config.local.R`
      added completes without error on Holden's machine.
- [ ] Run from the wrong folder stops with the anchor message; missing input stops before
      stage 1 with a table naming the file and where to get it.
- [ ] `git grep -nE "dchagaris|OneDrive"` over `*.R` matches only comments in
      `config.local.example.R`. `git grep -n "windows("` matches nothing.
- [ ] `output/basemaps/` has no diff after a full run (`git status --porcelain output/basemaps`
      empty); MD5s identical.
- [ ] Stages 2 and 3 run twice give identical MD5s for the MaxN rasters.
- [ ] `git ls-files --cached --ignored --exclude-standard` is empty.
- [ ] README and CLAUDE.md describe the four stages, inputs, outputs, and tested environment.
- [ ] Dave runs the branch with his three-line `config.local.R` and reports a clean
      `git status` on `output/basemaps/`.

## 8. GitHub workflow

- Issue #2 (this review) with sub-issue #3 (basemap redundancy, separate PR later).
- Branch `2-review-gfisher-repo` created from the issue's Development sidebar.
- Draft PR #4 against `main`, body carries `Fixes #2` and the checklist above.
- Commits are split by kind (code / docs / outputs / housekeeping), subject ends with
  `(issue #2)`, body explains why, trailer `Co-Authored-By: Claude Fable 5.1`.
- Merge with a merge commit once Dave approves; delete the branch.

## 9. Decisions log

1. **1 Oct 2026 (Holden):** Review the rework branch, merged into the review branch by
   fast-forward, rather than `main` alone.
2. **1 Oct 2026 (Holden):** Swap `xlsx` for `readxl`.
3. **1 Oct 2026 (Holden):** Untrack the legacy geodatabases and unused xlsx with
   `git rm --cached`; keep history as is.
4. **1 Oct 2026 (Holden):** Keep `data/April2026` as the default survey-input folder name.
5. **1 Oct 2026 (Holden):** Basemap-code relocation to `EcospaceBasemap` is issue #3, a
   separate PR after this review.
