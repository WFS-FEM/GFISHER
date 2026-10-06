# GFISHER issue #3: Remove redundant basemap code from GFISHER

Date opened: 2 Oct 2026 (work started 5 Oct 2026). Branch: `3-remove-redunant-basemap-code`.
Companion: `WFS-FEM/EcospaceBasemap#2`, branch `2-changes-made-with-gfisher`, draft PR
`WFS-FEM/EcospaceBasemap#3`.
Authors: Holden Harris, with Claude Code. Original code author of both pipelines: David Chagaris.

Sections 1 to 4 describe what was found. Sections 5 to 8 are the work plan. Section 9 logs
decisions as they are made. This file is updated as the work proceeds and travels with the PR.

## 1. Summary

Issue #3 asks that all code building the WFS Ecospace habitat basemaps live in
`EcospaceBasemap`, and that GFISHER keep only the survey-specific stages (video dataset, MaxN
maps, affinities) and consume basemap layers rather than build them. The same three goals as
the issue #2 review apply to both repositories: Holden can run it, Dave can run it with his
existing paths, any new GitHub user can run it.

The comparison of the two repositories (section 3) found that the overlap is not a copy and a
port. Dave wrote two different sum-to-1 habitat products on consecutive days: GFISHER stage 1
(`R/habitat_basemaps.R`, commit `c6de956`, 26 Aug 2026) and the EcospaceBasemap habitat
pipeline (`R/GFISHER functions.R` + `R/sum_to_1_functions.R`, commit `00cec74`, 27 Aug 2026).
They disagree on rock, on how unmapped cells are filled, on seagrass, on the depth cutoff, on
the geodatabase layer read, and on the layer set (nine layers versus ten). Their natural
low-relief reef layers correlate at r = 0.10 at 5 min. Removing GFISHER stage 1 outright
would change what the affinity stage reads, which is a method change and the author's call.

The work is therefore split:

- **Unconditional (this PR):** delete the retired legacy module, remove the copied dbSEABED
  download code, record provenance and MD5s for every input that EcospaceBasemap produces or
  ships (dbSEABED grids, seagrass rasters, depth grids), ship the 15-minute seagrass raster
  and regenerate the 15-minute basemaps, and let `config.local.R` point GFISHER at a sibling
  EcospaceBasemap clone instead of the shipped copies. The stage 1 method and the 5-minute
  basemaps do not change.
- **Gated on Dave (follow-up issue):** whether stage 1 is ported into EcospaceBasemap as a
  second product, dropped in favour of EcospaceBasemap's sum1 layers, or kept in GFISHER as
  affinity covariates (section 5.5).

Ecospace reads EcospaceBasemap's `output/5min/habitat/sum1/habitat_*_5min.asc` (decision 1),
so those ten files are the byte-identity target for the companion PR.

## 2. What the code does (main at `62f7a08`)

| Stage | File | Reads basemaps? |
|---|---|---|
| 1 | `R/habitat_basemaps.R` | builds nine layers in `output/basemaps/<res>min/` |
| 2 | `R/video_dataset.R` | no |
| 3 | `R/maxn_maps.R` | no (template is the depth grid) |
| 4a | `R/selection_ratio_affinities.R` | yes: reads eight of the nine layers from disk (`BASEMAP.SPEC`, SGR excluded) |
| 4b | `R/site_level_affinities.R` | no (geodatabase attributes) |
| 4c | `R/substrate_affinities.R` | no (raw dbSEABED grids at native resolution) |

Nothing uses the in-memory object stage 1 returns. Stage 4c needs the raw dbSEABED grids
whatever happens to stage 1. `R/legacy/habitat_maps_cellarea.R` (525 lines) is sourced by
nothing; a commented call in the driver pointed at it.

## 3. Findings: inventory of overlap

Line numbers refer to GFISHER `main` at `62f7a08` and EcospaceBasemap `main` at `3d97024`.

### 3.1 Basemap-building code

| Item | GFISHER | EcospaceBasemap | Verdict |
|---|---|---|---|
| Geodatabase read, reef per microgrid | `R/habitat_basemaps.R:76-110` `fn.microgrid_reef`: layer `Hab_Data_FINAL`, reef area / scanned area | `R/GFISHER functions.R:373-399`: layer `Hab_Data_Dissolve_Site_Selection`, reef area / scanned area | Different polygon layer, same denominator |
| Unmapped-cell fill | `fn.shrink_reef` `:260-302`: empirical-Bayes shrinkage toward the region x depth-bin mean (`k='auto'`, `target='stratum'`); reef forced to 0 below 300 m | `fn.fill_habitat_gaps` `:52-291`, driver `make_WFS_basemaps.R:174-179`: IDW idp 4 nmax 8, `anchor.zero='both'`, fill limited to 200 m, no cutoff on the final layers | Different method |
| dbSEABED to the grid | `fn.dbseabed_layers` `:316-355`: raw grids, `< 0` as NA, mean of native cells whose centres fall in each model cell, NA filled from the stratum mean | `fn.rasterize_dbseabed` `R/dbSEABED_functions.R:236-292`: raw grids, `-99` as NA, block `aggregate` then nearest `resample`; NA filled by a 3x3 focal mean in the sum step | Both keep NODATA; different aggregation and fill |
| Rock and sediment in the sum | RCK kept as its own layer, `UNC = 1 - rock`, every layer divided by the row sum (header point 4: reef and rock are uncorrelated, so neither is subtracted) | rock dissolved into NL/NM/NH by stratum relief shares (`fn.combine_habitats_sum1` `R/sum_to_1_functions.R:195-355`); Sand, Mud, Gravel renormalised among themselves to the remainder | Different method |
| Seagrass | FWC statewide raster, area-averaged, diluted by the row normalisation (`:393-400`) | `max(GulfwideSAV, FWC)` combined raster, takes `min(seagrass, 1 - reef)` | Different source and rule |
| Artificial reef database | not used | added to AL/AM/AH (`R/artificial_reef_functions.R`) | EcospaceBasemap only |
| Layer set | AL AM AH NL NM NH RCK UNC SGR (9) | AL AM AH NL NM NH seagrass Sand Mud Gravel (10) + `RCK_raw` outside the sum | Different products |

Measured at 5 min on the committed outputs of both repositories (both share the same 78 x 66
grid, water mask and depth values to 5e-12):

| Layer | GFISHER total (cell fractions) | EcospaceBasemap total | Correlation |
|---|---|---|---|
| NL | 50.6 | 328.0 | 0.10 |
| NM | 0.9 | 16.9 | 0.31 |
| NH | 10.3 | 55.9 | 0.36 |
| AL / AM / AH | 0.12 / 0.29 / 0.03 | 0.14 / 0.35 / 0.09 | 0.85 / 0.62 / 0.31 |
| seagrass | 47.5 (SGR) | 191.9 | 0.79 |

EcospaceBasemap has reef in 289 cells deeper than 300 m (rock-derived); GFISHER has none.
GFISHER's README (line 123-124) says its nine layers are "the versioned deliverable that
feeds Ecospace", which contradicts decision 1; the wording is corrected in this PR and the
question is put to the author.

### 3.2 Legacy cell-area functions

`R/legacy/habitat_maps_cellarea.R` holds `fn.make_GFISHER_habitat_maps` (live but called only
from a comment), `fn.make_GFISHER_habitat_maps_old` (dead) and `fn.plot_GFISHER_habitats`
(dead for current outputs). Its header says reef is divided by cell area; the code (line 219)
divides by the summed microgrid area, the same scanned-area denominator as everything else.
The real differences from stage 1 are the IDW fill (idp 4, nmax 8), the 200 m cutoff, and
the absence of rock, UNC, SGR and the sum-to-1. EcospaceBasemap's `R/GFISHER functions.R` is
the terra port of this file; `anchor.zero='both'` reproduces its behaviour, and its README
section 2.3 documents the port. The legacy module also carried the only uses of `terra`,
`gstat`, `lwgeom` and `maps` in GFISHER (listed as optional packages) and the only reason for
the `output/maps/5min/` and `output/maps/15min/` ignore rules.

### 3.3 Data shipped or pulled in both

| Item | GFISHER | EcospaceBasemap | Same? |
|---|---|---|---|
| `fn.pull_dbseabed` | `R/_setup.R:65-81` (copy) + driver download block `:79-87` | `R/dbSEABED_functions.R:45-72` (original) | Same URLs and layout |
| Raw dbSEABED grids | `data/dbseabed/` tracked, `SOURCE.md` | gitignored, downloaded | Byte-identical (MD5 GVL `5ce2c0b6`, MUD `a173130e`, RCK `46c0bd22`, SND `bf01968f`) |
| Seagrass raster | `data/seagrass/seagrass_5min.asc` | `output/5min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_5min.asc` | Byte-identical (MD5 `dff2ff98`); no 15-min version on either side |
| Bathymetry | `data/bathymetry/depth 5min 66x78.asc` (written by `raster`, `CELLSIZE 0.0833333333333333`, NODATA -9999), `depth 15min 22x26.asc`; two `excl layer *.asc` no code reads | `output/<res>min/depth/depth_<res>min.asc` (written by `terra`, `cellsize 0.083333333333`, `NODATA_value nan`) | Values identical to 5e-12; headers differ |
| Input manifest and check | `R/_setup.R` `fn.data_manifest`, `fn.check_inputs` | `R/data_setup_functions.R` | Parallel helpers for different inputs; keep both |

**Header precision (plan doc issue #2, R5).** Seven survey stations sit exactly on 5-minute
row boundaries. A template whose header cell size is rounded to 12 decimals assigns them to
a different row than the exact `0.0833333333333333` does, which moves the effort raster and
the MaxN maps. terra writes the ASCII header through GDAL's AAIGrid driver, which prints the
cell size with 12 decimals and offers no option to change it (tested 5 Oct 2026, terra
1.8.80, GDAL 3.11.4; `SIGNIFICANT_DIGITS` has no effect on the header). GFISHER therefore
keeps its own raster-written depth file as the template and reads EcospaceBasemap-written
grids only for their values.

### 3.4 Other observations

- `LEGACY.SPEC` in `R/selection_ratio_affinities.R:75-87` describes "the retired
  `input_ascii_sum1/` layers", but its patterns (`AL_prop.*\.asc$`, `gmf_RCK_val.*\.asc$`)
  match EcospaceBasemap's intermediate products `output/<res>min/habitat/gfisher/` and
  `output/<res>min/habitat/dbseabed/`. The comment is stale; the spec is kept and re-documented.
- `output/basemaps/` holds three figures no current code writes (`depth_cutoff_comparison.png`,
  `target_comparison_NL.png`, `target_comparison_NH.png`) and a duplicate
  `habitat_basemaps_5min.png`; they came in with `c6de956` and look like the author's
  method-comparison figures. Moved to `docs/img/` only if the author confirms.
- GFISHER tracks 15-minute basemaps built by the author from a `seagrass_15min.asc` that is
  not in either repository; EcospaceBasemap has never been run at 15 min. Its regions section
  reads a hand-edited 5-minute grid, so a 15-minute run there covers sections 1 to 2.5 only.

## 4. Evidence

### 4.1 Comparison of committed outputs (5 Oct 2026)

Base R comparison of `output/basemaps/5min/*.asc` (GFISHER) against
`output/5min/habitat/sum1/*.asc` (EcospaceBasemap), both on disk at the commits named in
section 3: identical NA masks (1,310 land cells), depth max abs diff 5.0e-12, row sums 1.0 in
both (EcospaceBasemap within 2e-6), reef class totals and correlations as tabulated in 3.1.
MD5s of the dbSEABED grids and the seagrass raster agree across repositories.

### 4.2 Verification runs

**EcospaceBasemap, 5 Oct 2026 (PR #3).** Full driver at 5 min after the code changes: of 141
pre-existing `.asc/.prj/.csv` files, 139 byte-identical, the two depth grids differ in header
and NoData token only (values identical, max diff 0); the ten sum1 MD5s unchanged. Two 15-min
runs (sections 1 to 2.5; section 5 stops on the hand-edited 5-min grid) agree on all 65
tracked files; `md5sum -c CHECKSUMS.md5` passes at both resolutions; the 15-min depth grid
matches GFISHER's `depth 15min 22x26.asc` to 5e-12 with the same land mask. Line endings:
with `core.autocrlf=true` a Windows checkout turned stored-LF grids into CRLF on disk
(GFISHER's `data/seagrass/seagrass_5min.asc` read MD5 45144e86 on disk against the stored
dff2ff98), so both repositories pin `.asc/.prj/.csv/.md5` to LF; `write.csv()` writes CRLF on
Windows, so EcospaceBasemap hashes the LF form of each file.

**GFISHER, 6 Oct 2026.** Stage 1 at `res = 15` (driver lines 1 to the stage 1 plot) with the
shipped `data/seagrass/seagrass_15min.asc`: against the author's committed 15-min layers, 10 of
428 water cells changed, exactly the cells where the seagrass fraction differs (max 0.028);
the other layers moved only through the row normalisation in those cells (NL up to 0.0021, UNC
up to 0.027, the rest below 3e-4); shrinkage constants identical; row sums exactly 1. Stored
dbSEABED grid MD5s (LF form, equal in both repositories): RCK e1332f6d, GVL 1242a335, SND
8147b226, MUD 16e40884; the CSDMS zips hold CRLF files whose hashes differ, which is why the
first draft of both `SOURCE.md` files had to be corrected. The 5-min full runs (shipped copies,
then `dir.ecospace.basemap`) are recorded below when complete.

**Fresh clone, 6 Oct 2026** (branch at G7, cloned into a temporary folder, no `config.local.R`):
every checked-out grid is LF on disk (28 files, `git ls-files --eol`); the eight on-disk MD5s
(seagrass 5 and 15 min, depth 5 and 15 min, four dbSEABED grids) equal the `SOURCE.md` tables,
and the two seagrass values equal EcospaceBasemap's `CHECKSUMS.md5` entries; `Rscript` on the
driver prints the input check with bathymetry, species list, dbSEABED and seagrass OK, lists
only the geodatabase and the three survey CSVs as missing, and stops before stage 1 (exit 1).

**5-min full run on the shipped copies, 6 Oct 2026:** stages 1 to 4a completed and 4b was
running when the one-hour job limit stopped it (two R jobs shared the machine, so the run was
far slower than the 16 minutes a lone run takes); up to that point no tracked `.asc` or `.csv`
under `output/` changed, only the three PDFs, which embed timestamps.

**5-min full run with `dir.ecospace.basemap` set, 6 Oct 2026:** all six stages, exit 0 after
61 minutes (alone on the machine for most of it); the input check printed the clone as the
basemap source and resolved dbSEABED and seagrass inside it; afterwards no tracked `.asc` or
`.csv` under `output/` differed from the branch, so the clone path and the shipped copies
produce the same basemaps, MaxN maps and affinities. The four PDFs regenerate with new embedded
dates and are otherwise the same size; they are not part of the comparison.

**Hygiene, 6 Oct 2026:** `git ls-files --cached --ignored --exclude-standard` empty; the only
machine path in a tracked `.R` file was in `docs/issue5_seed_experiment.R` (a temp folder and
a `setwd()` from PR #4), fixed in G8 to `tempdir()` and the repo-root check.

## 5. Fix design

### 5.1 Owner and interface

| Item | Owner | EcospaceBasemap produces | GFISHER reads | Override in `config.local.R` |
|---|---|---|---|---|
| Raw dbSEABED grids | EcospaceBasemap owns the download path; both repos track a copy | `data/dbseabed/Gmf_<CLS>/gmf_<CLS>_val.asc` | the same path under `data/` | `dir.dbseabed`, or `dir.ecospace.basemap` which resolves to `<clone>/data/dbseabed` |
| Seagrass raster | EcospaceBasemap | `output/<res>min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_<res>min.asc` | `data/seagrass/seagrass_<res>min.asc` (copy) | `file.seagrass`, or `dir.ecospace.basemap` |
| Depth grid | EcospaceBasemap makes the values; GFISHER keeps its raster-written copy as the template | `output/<res>min/depth/depth_<res>min.asc` | `data/bathymetry/depth <res>min <rows>x<cols>.asc`, unchanged | none, deliberately (R5) |
| sum1 layers | EcospaceBasemap | `output/<res>min/habitat/sum1/habitat_*_<res>min.asc` + `CHECKSUMS.md5` | nothing today (gated option B would read copies under `data/basemaps/`) | gated |
| Nine-layer stage 1 set | GFISHER, until the author decides | | `output/basemaps/<res>min/`, written and read in place | |

MD5s: EcospaceBasemap writes `output/<res>min/CHECKSUMS.md5` from its driver; GFISHER records
an MD5 table in each `data/*/SOURCE.md` (file, MD5, EcospaceBasemap path, EcospaceBasemap
commit) and the README's reproducibility check compares them.

### 5.2 GFISHER changes (this PR)

1. Delete `R/legacy/habitat_maps_cellarea.R`; drop the commented call, the four optional
   packages it alone needed, the two ignore rules, and every documentation mention; README
   caveat 5 records the last commit holding the file and points at EcospaceBasemap's port.
2. Delete `fn.pull_dbseabed` and the driver's download block. The tracked copy is the input;
   `git checkout -- data/dbseabed` restores it, `dir.dbseabed` or `dir.ecospace.basemap`
   redirect it, and EcospaceBasemap's `fn.pull_dbseabed()` is the download path.
3. `data/dbseabed/SOURCE.md`, `data/seagrass/SOURCE.md` and a new `data/bathymetry/SOURCE.md`
   gain MD5 tables and the EcospaceBasemap commit they were taken from; the bathymetry file
   explains why the template is kept (R5).
4. Ship `data/seagrass/seagrass_15min.asc` (from the companion PR) and regenerate
   `output/basemaps/15min/` with it. The 5-minute basemaps must not change.
5. New setting `dir.ecospace.basemap` (default `NULL`): when set, `dir.dbseabed` and
   `file.seagrass` resolve inside that clone unless given explicitly. `dir.ecospace.maps`
   keeps only its `dir.ewemaps` role. `fn.check_inputs` prints which source is in use.
6. Documentation: README, CLAUDE.md, `config.local.example.R`, `R/_setup.R` header, the
   `LEGACY.SPEC` comment, this plan.

### 5.3 EcospaceBasemap changes (companion PR, `WFS-FEM/EcospaceBasemap#2`)

Track `output/<res>min/depth/` and `output/<res>min/habitat/**` (`.asc/.prj/.csv`), write
`CHECKSUMS.md5` per resolution from the driver, allow the dbSEABED writer to overwrite, write
depth with `NAflag=-9999`, track `data/dbseabed/*.asc` with a `SOURCE.md`, run and commit the
15-minute grids for sections 1 to 2.5, document the output tree and the downstream contract.
The ten 5-minute sum1 layers must regenerate byte-identically before any output is committed.

### 5.4 Verification

- EcospaceBasemap: the ten sum1 MD5s before and after the companion PR are identical; the
  15-minute outputs reproduce across two runs; `md5sum -c CHECKSUMS.md5` passes.
- GFISHER: a full 5-minute run on this branch leaves `output/basemaps/5min/`, `output/maps/`
  and `output/affinity_*_mice/` byte-identical to `main`; a second run with
  `dir.ecospace.basemap` set gives the same result; every `SOURCE.md` MD5 matches the file on
  disk and EcospaceBasemap's `CHECKSUMS.md5`.
- Fresh clones of both repositories: GFISHER with no `config.local.R` stops at the input
  check naming only the FWRI files; with a three-line `config.local.R` it runs to completion
  under `Rscript`. EcospaceBasemap's `fn.check_inputs()` reports only the geodatabase and the
  artificial-reef structures missing.
- Hygiene: `git ls-files --cached --ignored --exclude-standard` empty in both; no machine
  paths in tracked `.R` files.

### 5.5 Gated on the author: the fate of stage 1

| | A. Port stage 1 into EcospaceBasemap as a second product | B. Drop stage 1; 4a reads sum1 | C. Keep stage 1 in GFISHER as affinity covariates |
|---|---|---|---|
| EcospaceBasemap | `R/habitat_basemaps.R` moved in with namespaced `raster::`/`sp::`/`FNN::` calls, a driver section writing `output/<res>min/habitat/sum1_gfisher/` | nothing further | README notes the nine-layer set is GFISHER's |
| GFISHER | remove `R/habitat_basemaps.R` and the stage 1 block; track copies under `data/basemaps/<res>min/`; `dir.hab` resolves copies or the clone; drop `data/seagrass/` and `FNN` | as A, plus a `SUM1.SPEC` in `R/selection_ratio_affinities.R` (reef family AL..NH, sediment family Sand/Mud/Gravel, seagrass excluded, no RCK) | documentation only; optional rename of `output/basemaps/` |
| Guard (A and B) | effort raster built on the depth template rather than `hab[[1]]`; `fn.load_layer_stack` checks `compareRaster(rowcol=TRUE, extent=FALSE)` and re-stamps the extent | same | not needed |
| Verification | layers byte-identical if written with `raster`; otherwise values within 5e-7 with header-only differences; stages 2 to 4 byte-identical | stages 2, 3, 4b, 4c byte-identical; 4a becomes a new baseline | nothing moves |
| Risk | two sum-to-1 products in one repository | affinity results change; seagrass exclusion must be re-checked on the combined raster | identical layer names in two repositories |

Questions put to the author on issue #3: (1) which product feeds Ecospace today; (2) which
seagrass raster the basemap should use; (3) whether the nine-layer set survives and where;
(4) whether the three comparison figures are his and can move to `docs/img/`; (5) whether the
`excl layer *.asc` files are still needed.

## 6. Implementation steps

EcospaceBasemap first (GFISHER copies its outputs and cites its commit), then GFISHER.
Commits are one logical change each, subject ending `(issue #N)`.

EcospaceBasemap (done 5 Oct 2026, PR #3): E1 ignore rules; E2 checksum helpers + driver call;
E3 dbSEABED writer overwrite; E4 depth NoData flag; 5-minute regeneration and MD5 comparison;
E5 commit the 5-minute grids; three fixes found on the way (LF pin via `.gitattributes`,
`dir.basemaps` resolved after `config.local.R`, checksums taken on the LF form of each file);
E6 ship the dbSEABED grids; two 15-minute runs; E7 commit the 15-minute grids; E8 documentation.

GFISHER (6 Oct 2026): G0 this plan; G1 retire the legacy module; G2 stop downloading dbSEABED;
G3 dbSEABED provenance; G3b pin LF line endings (`.gitattributes`, same reason as in
EcospaceBasemap); G4 ship the 15-minute seagrass raster; G5 bathymetry provenance; G6
sibling-clone override; G7 regenerate the 15-minute basemaps; 5-minute full runs both ways;
G8 documentation; G9 move the comparison figures (only if the author confirms); G10 record
verification.

## 7. Acceptance criteria

- [x] `git diff main -- output/basemaps/5min` is empty; `output/maps/` and
      `output/affinity_*_mice/` are byte-identical to `main` after a full run (4.2).
- [x] Every `data/*/SOURCE.md` MD5 matches the file on disk and EcospaceBasemap's
      `CHECKSUMS.md5` (fresh clone, 4.2).
- [x] A run with `dir.ecospace.basemap` set produces the same outputs as the shipped copies (4.2).
- [x] A fresh clone without `config.local.R` stops at the input check naming the FWRI files
      (4.2); the full run under `Rscript` is the override run above, in the working clone.
- [x] `git ls-files --cached --ignored --exclude-standard` is empty;
      `git grep -nE "dchagaris|OneDrive" -- '*.R'` matches only `config.local.example.R`.
- [ ] The author's answers to the five questions are recorded in section 9 and the gated
      option is tracked in a follow-up issue.
- [ ] Dave runs the branch with his `config.local.R` and confirms the outputs.

## 8. GitHub workflow

- Issue #3 (sub-issue of #2); companion issue `WFS-FEM/EcospaceBasemap#2`.
- Branches `3-remove-redunant-basemap-code` (GFISHER) and `2-changes-made-with-gfisher`
  (EcospaceBasemap), both created from the issues' Development sidebars.
- Draft PRs against `main` in each repository, bodies carry `Fixes #N`, the checklist above
  and a link to the other PR. EcospaceBasemap merges first.
- Commits split by kind (code / docs / outputs / data / housekeeping), subject ends with
  `(issue #N)`, body explains why, trailer `Co-Authored-By: Claude Fable 5.1`.

## 9. Decisions log

1. **2 Oct 2026 (Holden):** The basemaps Ecospace reads are EcospaceBasemap's
   `output/5min/habitat/sum1/` layers; they are the byte-identity target.
2. **2 Oct 2026 (Holden):** EcospaceBasemap PR #1 (`feature/portable-data-setup`) merged into
   `main` (`3d97024`); the companion branch starts from `main`.
3. **5 Oct 2026 (Holden):** Split issue #3: unconditional removals in this PR; the fate of
   stage 1 (port, drop, or keep) waits for the author's answer and gets a follow-up issue.
4. **5 Oct 2026 (Holden):** Delete `R/legacy/habitat_maps_cellarea.R`; README records the
   last commit holding it.
5. **5 Oct 2026 (Holden):** The 15-minute grid is in scope for both repositories.
6. **5 Oct 2026 (Holden):** Interface is "ship copies plus sibling override": EcospaceBasemap
   tracks its small deliverables; GFISHER tracks copies with MD5s and can point at a clone.
7. **5 Oct 2026 (Holden):** Holden creates the EcospaceBasemap issue from drafted text and
   posts the comment to the author on issue #3.
