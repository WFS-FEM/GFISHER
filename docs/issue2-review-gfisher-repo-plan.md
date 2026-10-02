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
| R5 | `output/affinity_selratio_mice/GFISHER_survey_effort_5min_66x78.asc` | The committed effort raster was built on a template whose header has `CELLSIZE 0.0833333333329999` (the author's older Ecospace maps tree), while every other grid in the pipeline uses the exact 1/12 degree. Seven stations sit exactly on a 5-minute latitude line and fall on opposite sides of the boundary under the two templates. The current code is self-consistent (template = basemaps = depth grid); the committed file is simply from an older template. Found by the baseline run, section 4.1. |

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

Run on 1 Oct 2026 on Holden's machine (Windows 11, R 4.5.1, raster 3.6-32, terra 1.8-80,
sf 1.0-21) with `Rscript` on a copy of the driver at commit `b544821` in which exactly two
lines were changed: `dir.dbseabed` pointed at `EcospaceBasemap/data/dbseabed` and
`file.seagrass` at `EcospaceBasemap/output/5min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_5min.asc`.
`data/April2026` is a Windows directory junction to the OneDrive copy of Dave's folder, so
the files are not duplicated and git ignores the path. Dave's OneDrive Ecospace maps tree is
not synced to Holden's machine and was not needed.

**Result: all stages ran to completion. Exit code 0. Wall time 15 min 39 s.**

| Stage | Ran | Log evidence |
|---|---|---|
| 1 basemaps | yes | 309,348 microgrids, 141,053 habitat polygons; 1,495 of 3,838 water cells mapped; QC row sum min 1 max 1; wrote 9 layers |
| 2 video dataset | yes | "Dropping 52471 record(s) with missing modnumber, maxn, or coordinates" |
| 3 MaxN maps | yes | 19 group rasters + PDF written |
| 4a selection ratios | yes | red-grouper-1 pooled with red-grouper-0 (3 cells with MaxN>0); 4 files written |
| 4b site affinities | yes | 4 files written |
| 4c substrate affinities | yes | 5 files written |

Warnings only: `rm(.SavedPlots)` object not found (P3); packages built under a newer R
patch release; one GDAL `organizePolygons()` performance message while reading the
geodatabase. `windows(record=T)` ran under `Rscript` on Windows without error.

**Comparison with the committed outputs** (MD5 of all 54 tracked `.asc`/`.csv` files
before and after; no line-ending-only differences were found):

| Output group | Identical | Changed | Cause |
|---|---|---|---|
| `output/basemaps/15min/` (10 files) | 10 | 0 | not regenerated at `res = 5` |
| `output/basemaps/5min/` (10 files) | 0 | 10 | seagrass input only, see below |
| `output/maps/.../maxn/mice/` (19 rasters) | 8 | 11 | unseeded length draws (R1): the 8 unchanged are the single-stanza groups plus red-grouper-0; the 11 changed are the gag and red grouper age stanzas |
| `output/affinity_selratio_mice/` (5 files) | 0 | 5 | downstream of the above, plus the effort raster (R5) |
| `output/affinity_site_mice/`, `affinity_substrate_mice/` (9 files) | 0 | 9 | downstream of stage 2 |
| `output/affinity_selratio/` (untagged, 5 files) | 5 | 0 | stale folder, not written by current code (R4) |

**Basemaps: the seagrass raster is the only source of difference.** The 45 water cells
where the SGR layer differs are the only cells where any layer differs; in the other 3,793
water cells all nine layers are bit-identical. In those 45 cells the reef and rock layers
scale by exactly `(1 - SGR_new) / (1 - SGR_old)` (max deviation 7e-8), which is the sum-to-1
normalisation. Largest absolute differences: UNC 0.080, SGR 0.080, NL 0.006, others below
1e-3. The committed SGR mean is 0.01218; ours 0.01236. So the geodatabase read, the dbSEABED
grids, and the whole stage 1 algorithm reproduce exactly on Holden's machine; what is
missing is Dave's `seagrass_5min.asc` (and `seagrass_15min.asc`), which differ slightly
from the EcospaceBasemap `Seagrass_Statewide` raster in 45 cells. **Ask of Dave:** add the
two files (about 44 KB each) to `data/seagrass/` so the basemaps can be reproduced exactly.

**New finding R5, survey effort raster.** `GFISHER_survey_effort_5min_66x78.asc` keeps the
same station total (14,613) but 14 cells differ by one station each, in seven pairs of
vertically adjacent cells. Seven survey stations in the domain have a latitude exactly on a
5-minute row boundary (28.50, 30.00, 26.50, 29.75, 30.00, 26.75, 28.25 deg N); the env file
has 11 such stations and 25 on column boundaries overall. **Cause (settled 2 Oct 2026):**
the committed file's header reads `CELLSIZE 0.0833333333329999`, while the basemaps, the
MaxN maps and the depth grid all carry the exact `0.0833333333333333`. Dave's effort raster
was therefore built on a template from his older Ecospace maps tree (the `input_ascii_sum1`
layers), whose rounded cell size places those seven stations just above a row boundary;
with the exact template they sit on the boundary and go to the row below. The same seven
stations carry MaxN records in the single-stanza groups, and those maps are identical
between Dave's run and ours, which confirms that the current code assigns boundary points
consistently when the template is the same. No code change is needed; the committed file
will simply be regenerated from the current template. A guard that the effort template
matches the depth grid is cheap insurance and is included in Phase 3.

A copy of the baseline outputs is kept outside the repo at
`%LOCALAPPDATA%\Temp\gfisher_baseline\outputs_baseline_run\` with the MD5 tables
(`baseline_md5_committed.csv`, `baseline_md5_compare.csv`) and the full log.

### 4.1b Stochasticity check: is the seed the only thing moving the stanza maps?

Asked by Holden on 2 Oct 2026 after the portability run reproduced the 8 single-stanza MaxN
maps exactly but not the 11 gag / red grouper age-stanza maps. Stages 2 and 3 were rerun on
their own (no geodatabase needed) three times: seed 1, seed 1 again, and seed 2.

| Test | Result |
|---|---|
| Same seed twice | Stage 2 tables `identical()`; all 19 rasters byte-identical |
| Grand total MaxN, seed 1 vs seed 2 | 1,043,543 in both |
| Per-species total MaxN | identical for every species |
| Per-station x species total MaxN (60,374 rows) | 0 rows differ |
| Union of occupied cells per species | gag 458 cells, red grouper 747 cells, identical sets |
| 7 single-stanza maps | identical across seeds |
| Rows with a missing group number | 17,215 rows (MaxN 139,311) across 160 taxa, identical across seeds; these are taxa outside the model groups (e.g. `baitfish unk`, `pagrus pagrus`), not multistanza fish. No gag or red grouper record lacks a group |
| Function's own "counts did not sum back" check | silent for both seeds |

So the seed changes exactly one thing: which age stanza each individual gag or red grouper is
assigned to, through the random length draw (`rtruncnorm`) or the random pairing of observed
lengths with individuals (`sample.int`). Nothing is gained or lost.

**The per-cell differences are not small, and that is inherent to the method, not a bug.**
Cell values are the mean MaxN per record, so moving one fish between stanzas changes a
sparse cell by about 1. Seed 1 vs seed 2 gives, per stanza map, 13 to 229 cells changed,
mean |difference| 0.3 to 1.0, max 1 to 4, and correlations from 0.14 (red grouper 1, 10
occupied cells) to 0.98 (red grouper 5+, 695 cells). Dave's committed maps vs our baseline
run show the same pattern (e.g. gag 1: 114 cells, mean 0.97, r = 0.79), confirming they are
two draws from the same process. Red grouper 0 is identical across all runs: its assignment
does not depend on the draw.

Implication for the fix: `set.seed()` makes the stanza maps reproducible but does not make
them less noisy. A single draw is one realisation of the stanza split. If the stanza maps
matter downstream (they feed stage 4), a deterministic alternative is to average the maps
over many seeds (the expected stanza composition) or to assign the expected fraction of each
record to each stanza instead of drawing. **This is a methods decision for Dave, tracked in issue #5 (seed-sensitivity experiment
across 10 seeds, compared with the bootstrap intervals); not changed in this PR.**

### 4.1c Issue #5 experiment: ten seeds through stages 2, 3 and 4a (2 Oct 2026)

Setup: seeds 1 to 10; stage 2 and 3 as in the driver (`fun=mean`, `background=0`); stage 4a
against the committed `output/basemaps/5min` with 1,000 bootstrap draws per group, no prior
constraints; outputs written outside the repo. 1.7 min per seed. Notes log and per-seed
tables in `%LOCALAPPDATA%\Temp\gfisher_baseline\exp_issue5\`; summary in
`summary_across_seeds.csv` (96 rows = 11 stanza groups x 8 layers).

**Control.** The 7 single-stanza groups and red grouper 0 have exactly the same selection
ratio and affinity in all ten seeds (range 0). Whatever moves below is the stanza draw alone.

**How to read the numbers.** For each stanza group and habitat layer, "seed range" is the
spread of the selection ratio across the ten seeds, and "bootstrap width" is the 95% interval
the method already reports for sampling uncertainty. Ratio < 1 means the seed adds less
uncertainty than the method already acknowledges. A is the 0 to 1 affinity Ecospace uses;
"best layers" counts how many different layers came out on top (A = 1) across the ten seeds.

| Group | MaxN>0 cells | seed range / bootstrap width, median (max) | A range, median (max) | distinct best layers in 10 seeds | rank agreement with seed 1, min Spearman |
|---|---|---|---|---|---|
| gag 0 | 39 | 0.41 (0.95) | 0.35 (0.95) | 3 | 0.36 |
| gag 1 | 145 | 0.56 (0.92) | 0.18 (0.41) | 4 | 0.76 |
| gag 2 | 190 | 0.49 (0.63) | 0.11 (0.26) | 4 | 0.86 |
| gag 3 | 192 | 0.49 (0.95) | 0.15 (0.26) | 4 | 0.76 |
| gag 4 | 176 | 0.35 (0.72) | 0.14 (0.34) | 2 | 0.81 |
| gag 5+ | 342 | 0.31 (0.51) | 0.08 (0.17) | 1 | 0.98 |
| red grouper 0 | 108 | 0 (0) | 0 (0) | 1 | 1.00 |
| red grouper 1 | 10 | 2.23 (152) | 0.95 (0.97) | 4 | -0.25 |
| red grouper 2 | 41 | 0.79 (11.4) | 0.41 (0.96) | 3 | 0.24 |
| red grouper 3 | 162 | 0.40 (0.67) | 0.18 (0.69) | 2 | 0.76 |
| red grouper 4 | 262 | 0.41 (0.56) | 0.18 (0.33) | 4 | 0.81 |
| red grouper 5+ | 695 | 0.14 (0.21) | 0.03 (0.08) | 3 | 0.95 |

Across all 96 group x layer rows, 88 have a seed range smaller than the bootstrap width
(median ratio 0.43). The 8 exceptions all belong to red grouper 1 and red grouper 2.

**Plain-language reading.**
- For the well-sampled stanzas (gag 5+, red grouper 5+, red grouper 3 and 4, gag 1 to 4),
  the seed moves the affinities by less than half of the uncertainty the method already
  reports, the ranking of habitats barely changes, and the top habitat is stable or nearly so.
  These outputs are robust to the draw.
- For the two sparse young stanzas (red grouper 1 with 10 occupied cells, red grouper 2 with
  41), the seed dominates: the affinity of a layer can swing from 0 to 1 and the "best"
  habitat changes with the seed (red grouper 1: NL, NH, AH or AL depending on the draw).
  These two groups' affinities are not robust as currently estimated. The code already pools
  red grouper 1 with red grouper 0 for being sparse; red grouper 2 (41 cells) is not pooled.
- gag 0 (39 cells) sits in between: the ratios are within the bootstrap width, but the top
  habitat flips between NH and NL.

**What this suggests for the decision in issue #5.** A fixed seed is sufficient for the
well-sampled stanzas. For gag 0, red grouper 1 and red grouper 2 the answer depends on the
draw, so either (a) average the stanza maps over many seeds before estimating affinities (the
expected stanza composition), (b) raise the pooling threshold so stanzas with fewer than
about 50 occupied cells borrow from their neighbour, or (c) assign each record's expected
stanza fractions instead of drawing. Each gives a deterministic answer; which is appropriate
is a modelling judgement for Dave.

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
