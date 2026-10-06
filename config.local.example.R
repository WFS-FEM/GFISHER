# Local path overrides for `process GFISHER data.R`.
#
# Copy this file to config.local.R in the repo root and edit it. config.local.R is
# gitignored, so machine-specific paths never reach the repository. This is the only
# file you should need to touch to run the pipeline on another machine.
#
# The driver sources it AFTER setting its repo-relative defaults and BEFORE reading any
# input, so uncomment only the lines you want to change. `dir.gfisher` (the repo root)
# and `res` are already defined at that point and can be used below.
#
# Inputs that cannot ship with the repo (and how to get them) are listed in README.md
# under "Getting the data".


# Grid resolution in arc-minutes: 5 or 15 (depth grids for both ship in data/bathymetry/).
# Default: 5
# res <- 15


# Species grouping scheme: 'mice' or 'ecospace' (see the SPECIES GROUPING SCHEME block).
# Default: 'mice'
# group.scheme <- 'ecospace'


# Seed for the random length draws in stage 2 (fn.make_gfisher_videodataset). Any fixed
# integer makes stages 2-4 reproducible run to run; NULL restores the unseeded behaviour.
# Default: 1
# seed <- 1


# Folder holding the inputs that cannot ship: the FWRI geodatabase and the three 3LABS
# survey CSVs (maxn/env/lens 3LABS_93to24.csv). Point this at a shared drive or a OneDrive
# sync to avoid keeping a second 320 MB copy. (The repo default matches the author's layout.)
# Default: file.path(dir.gfisher, 'data', 'April2026')
# dir.data <- file.path('C:/Users/<you>/University of Florida',
#                       'Chagaris, David - WFS Fisheries Ecosystem Modeling',
#                       'data/GFISHER/April2026')


# The geodatabase, if it is not inside dir.data.
# Default: NULL, which finds the single GFISHER_EAST_Universe*.gdb inside dir.data.
# file.gdb <- 'C:/path/to/GFISHER_EAST_Universe_2026.gdb'


# RAW dbSEABED grids (Gmf_<CLS>/gmf_<CLS>_val.asc for RCK, GVL, SND, MUD) used by stages 1
# and 4c. They ship with the repo (data/dbseabed/, 4.4 MB, MD5s in SOURCE.md); override only
# to use another copy. If the folder is ever emptied, `git checkout -- data/dbseabed` restores
# it; EcospaceBasemap's fn.pull_dbseabed() is the download path.
# Default: file.path(dir.gfisher, 'data', 'dbseabed')
# dir.dbseabed <- 'C:/Repos/WFS-FEM/EcospaceBasemap/data/dbseabed'                      # a sibling repo
# dir.dbseabed <- 'C:/dchagaris/GitHub/WFS-FEM/EnvironmentalDrivers2EwE/data/dbSEABED'  # the author's layout


# Seagrass raster on the model grid, used by stage 1. Optional: without it the SGR layer is
# zero and the basemaps will not match the committed ones.
# Default: NULL, which resolves to data/seagrass/seagrass_<res>min.asc, or, when
# dir.ecospace.maps is set, to <dir.ecospace.maps>/input_ascii_sum1/<res>min/seagrass_<res>min.asc
# file.seagrass <- 'C:/Repos/WFS-FEM/EcospaceBasemap/output/5min/habitat/seagrass/seagrass_coverage_Seagrass_Statewide_5min.asc'


# The author's external Ecospace maps tree. Setting it reproduces the former lookups:
# seagrass from <tree>/input_ascii_sum1/<res>min/ and MaxN heatmaps written to <tree>/GFISHER/.
# Default: NULL (everything stays inside the repo)
# dir.ecospace.maps <- 'C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps'


# Where the stage 3 MaxN heatmaps are written.
# Default: file.path(dir.gfisher, 'output', 'maps'), or dir.ecospace.maps when that is set
# dir.ewemaps <- 'D:/Ecospace/maps'


# Depth grids. Default: file.path(dir.gfisher, 'data', 'bathymetry')  (ships with the repo)
# dir.bathy <- file.path(dir.ecospace.maps, 'bathymetry')
