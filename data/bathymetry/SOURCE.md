# Depth grids (tracked copies, kept as the template)

| File | Grid | MD5 (LF form) | Read by |
|---|---|---|---|
| `depth 5min 66x78.asc` | 78 x 66 cells, 5 arc-min, lon -87.5 to -81, lat 25 to 30.5, NODATA -9999 | `0418033b87d27b82ff79d0227ff5d0c5` | the driver, as `depth`: template for every stage |
| `depth 15min 22x26.asc` | 26 x 22 cells, 15 arc-min, same extent | `166a46ed6a46887070991ca2b6b1c386` | the same at `res = 15` |
| `excl layer 5min 66x78.asc` | 1 where depth > 500 m, else NODATA | `532230b93723576abe823561ac948325` | nothing in this repo |
| `excl layer 15min 22x26.asc` | the same at 15 arc-min | `28f8227d288643fe28721b174d9fe863` | nothing in this repo |

Depth is positive metres below the surface; land is NODATA. The values are the ETOPO 2022
bathymetry, fetched by `marmap::getNOAA.bathy()` and processed (land to NA, sign flipped,
exact zeros bumped to the shallowest positive depth) by the `WFS-FEM/EcospaceBasemap`
repository, section 1 of `make_WFS_basemaps.R`. These files were written by the original
author with `raster::writeRaster()`; EcospaceBasemap now tracks its own `terra`-written
versions, `output/<res>min/depth/depth_<res>min.asc` (MD5 `67040dbe0884f6ca344918a000359c0d`
at 5 min, `b1ed354cf5c0cce2153ca5126675fd08` at 15 min, EcospaceBasemap PR #3), whose values
agree with these to 5e-12 at both resolutions with the same land mask (GFISHER issue #3
plan document, section 4).

**Why GFISHER keeps its own copy instead of reading EcospaceBasemap's.** The two files differ
only in their header: `raster` writes the exact cell size `0.0833333333333333`, while `terra`
goes through GDAL's ASCII-grid driver, which prints twelve decimals, `0.083333333333`, and
offers no option to print more (tested October 2026, terra 1.8.80, GDAL 3.11.4). Seven video
stations sit exactly on 5-minute row boundaries, and the rounded cell size puts them in the
row above where the exact one puts them, which moves the survey-effort raster and the MaxN
maps (issue #2 plan document, finding R5). The depth grid is the template every other grid
is built on, so it stays the raster-written file. EcospaceBasemap-written grids (seagrass,
and in future the basemap layers) are read for their values only and never used as the
template.

`dir.bathy` in `config.local.R` points the driver at another folder holding files of the same
names. MD5s are of the LF form, which is how git stores and checks out `.asc` files here.
