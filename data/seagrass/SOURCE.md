# Seagrass raster on the 5-minute model grid (tracked copy)

`seagrass_5min.asc`: fraction of each 5-arc-minute Ecospace cell covered by seagrass
(66 x 78 cells, lon -87.5 to -81, lat 25 to 30.5, NODATA -9999). Read by stage 1
(`fn.make_habitat_basemaps`, argument `file.sgr`) and area-averaged into the `SGR` layer of
the sum-to-1 habitat basemaps.

Source: FWC "Seagrass Habitat in Florida" polygons,
https://geodata.myfwc.com/datasets/myfwc::seagrass-habitat-in-florida (public), rasterised to
the model grid by the `WFS-FEM/EcospaceBasemap` repository as
`output/5min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_5min.asc`. This file is a
byte-identical copy of that output, taken 2 Oct 2026.

Why it is here: the basemaps the original author committed were built from his own
`seagrass_5min.asc` (an earlier rasterisation of the same FWC layer, not available in this
repo). This copy reproduces those basemaps in all but 45 of 3,838 water cells, where the
seagrass fraction differs by up to 0.08. Shipping it lets a fresh clone build stage 1 without
any external path. There is no 15-minute version yet, so the committed 15-minute basemaps
remain those built by the author. EcospaceBasemap is the authoritative producer of this layer;
the duplication is tracked in issue #3.
