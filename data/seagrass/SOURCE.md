# Seagrass rasters on the model grid (tracked copies)

Fraction of each Ecospace cell covered by seagrass, read by stage 1 (`fn.make_habitat_basemaps`,
argument `file.sgr`) and area-averaged into the `SGR` layer of the sum-to-1 habitat basemaps.
Optional: without the file for the current `res`, `SGR = 0`.

| File | Grid | MD5 (LF form) | EcospaceBasemap source file |
|---|---|---|---|
| `seagrass_5min.asc` | 78 x 66 cells, lon -87.5 to -81, lat 25 to 30.5, NODATA -9999 | `dff2ff98fea2cfdd4e5166a1aa145561` | `output/5min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_5min.asc` |
| `seagrass_15min.asc` | 26 x 22 cells, same extent | `fb9370fa52f7807eed115343f99bda42` | `output/15min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_15min.asc` |

Source: FWC "Seagrass Habitat in Florida" polygons,
https://geodata.myfwc.com/datasets/myfwc::seagrass-habitat-in-florida (public), rasterised to
the model grid (`terra::rasterize(cover = TRUE)`) by the `WFS-FEM/EcospaceBasemap` repository,
section 2.1 of `make_WFS_basemaps.R`. Both files are byte-identical copies of EcospaceBasemap's
tracked outputs at its commit `e887a1b` (pull request #3), and the same MD5s appear in its
`output/<res>min/CHECKSUMS.md5`. The 5-minute copy was taken on 2 Oct 2026, the 15-minute
copy on 6 Oct 2026. Set `dir.ecospace.basemap` in `config.local.R` to read them from a clone
of EcospaceBasemap instead of these copies.

Note that EcospaceBasemap's own sum-to-1 basemap uses a different seagrass layer, the cell-wise
maximum of this FWC layer and the NOAA GulfwideSAV layer
(`seagrass_coverage_combined_<res>min.asc`). Which layer the basemap should use is one of the
questions put to the author under issue #3; stage 1 here keeps the FWC-only layer the author
used.

Why the copies are here: the basemaps the original author committed were built from his own
`seagrass_5min.asc` (an earlier rasterisation of the same FWC layer, not available in this
repo). The 5-minute copy reproduces those basemaps in all but 45 of 3,838 water cells, where
the seagrass fraction differs by up to 0.08. The author's `seagrass_15min.asc` was never
available either, so the 15-minute basemaps committed before issue #3 were his originals; from
issue #3 on they are rebuilt from the 15-minute copy and differ from his in the seagrass cells.
