# process GFISHER data.R -- driver for the GFISHER habitat / MaxN / affinity pipeline.
#
# Every path below is repo-relative, so a fresh clone runs without editing this file once the
# inputs that cannot ship (geodatabase, survey CSVs) are in place. To read them from, or write
# outputs to, somewhere else, put a config.local.R in the repo root: it is gitignored, and
# config.local.example.R shows what goes in it. README.md > Getting the data lists every input.
#
# Written to be stepped through interactively: each stage leaves its result in the workspace.
rm(list=ls()); graphics.off(); gc()

# Everything downstream is resolved from the working directory, so check it first.
if(!file.exists('GFISHER.Rproj'))
  stop('Set the working directory to the repo root - open GFISHER.Rproj, or setwd() there. ',
       'Currently: ', getwd())

source(file.path('R','_setup.R'))             # package check, plot device, input manifest
fn.check_packages()                           # stops with an install.packages() line if any is missing
fn.plot_device()                              # recording plot window when interactive on Windows

# Map-building code is split by stage. The former 'R/GFISHER functions.R' was divided into
# video_dataset.R / maxn_maps.R, with its cell-area habitat maps retired to R/legacy/.
# Every call in these files is namespaced (raster::, sf::), so nothing needs terra attached;
# the depth template is read with raster::raster() below. (terra is only used by R/legacy/.)
source(file.path('R','video_dataset.R'))      # STAGE 2: station x group MaxN table
source(file.path('R','maxn_maps.R'))          # STAGE 3: per-group MaxN heatmaps
source(file.path('R','habitat_basemaps.R'))   # STAGE 1: sum-to-1 habitat basemaps

#=========================== SETTINGS (repo-relative defaults) ======================================
# Change these in config.local.R, not here. The config is sourced right after this block.
dir.gfisher  <- getwd()
res          <- 5          # map resolution in arc-minutes; 5 and 15 ship in data/bathymetry/
group.scheme <- 'mice'     # species grouping scheme; see the SPECIES GROUPING SCHEME block
seed         <- 1          # seed for the stage 2 length draws; NULL = unseeded (not reproducible)

# Inputs that cannot ship with the repo (too large / FWRI data): the GFISHER East Universe
# geodatabase and the three 3LABS survey CSVs. dir.data is the folder that holds them.
dir.data     <- file.path(dir.gfisher,'data','April2026')
file.gdb     <- NULL       # NULL finds the single GFISHER_EAST_Universe*.gdb inside dir.data

# Inputs that ship with the repo.
dir.bathy    <- file.path(dir.gfisher,'data','bathymetry')
file.spplist <- file.path(dir.gfisher,'data','Master Species List.xlsx')

# RAW dbSEABED grids (Gmf_<CLS>/gmf_<CLS>_val.asc at their native 1.2 arc-min), read by stages 1
# and 4c. NOT the processed 5-min gmf_*_prop_*.asc layers, which renormalize rock against the
# grain-size triangle and convert the -99 NODATA flag to zero. Downloaded from CSDMS on first
# run if the default folder is empty.
dir.dbseabed <- file.path(dir.gfisher,'data','dbseabed')

# Seagrass raster on the model grid, used by stage 1 (optional: without it SGR = 0).
# NULL resolves to data/seagrass/seagrass_<res>min.asc, or to the author's Ecospace maps tree
# when dir.ecospace.maps is set.
file.seagrass     <- NULL
dir.ecospace.maps <- NULL  # the author's external Ecospace maps tree; sets file.seagrass + dir.ewemaps

# Where outputs go. dir.ewemaps holds the stage 3 MaxN heatmaps (NULL = dir.maps).
dir.maps     <- file.path(dir.gfisher,'output','maps')
dir.ewemaps  <- NULL

#--------------------------- local overrides --------------------------------------------------------
if(file.exists('config.local.R')){
  source('config.local.R')
  message('Applied local overrides from config.local.R')
}

#--------------------------- resolve derived paths and check inputs ---------------------------------
if(is.null(file.gdb))      file.gdb <- fn.find_gdb(dir.data)
if(is.null(file.seagrass)) file.seagrass <- if(!is.null(dir.ecospace.maps))
  file.path(dir.ecospace.maps,'input_ascii_sum1',paste0(res,'min'),paste0('seagrass_',res,'min.asc')) else
  file.path(dir.gfisher,'data','seagrass',paste0('seagrass_',res,'min.asc'))
if(is.null(dir.ewemaps))   dir.ewemaps <- if(!is.null(dir.ecospace.maps)) dir.ecospace.maps else dir.maps
file.maxn  <- file.path(dir.data,'maxn3LABS_93to24.csv')
file.env   <- file.path(dir.data,'env3LABS_93to24.csv')
file.len   <- file.path(dir.data,'lens3LABS_93to24.csv')
file.depth <- list.files(dir.bathy, pattern=paste0('^depth ',res,'min.*\\.asc$'), full.names=TRUE)
if(length(file.depth)!=1) stop("Expected one 'depth ",res,"min*.asc' raster in ",dir.bathy,
                               ", found ",length(file.depth))

# Public dbSEABED grids download themselves into the default folder only; a custom
# dir.dbseabed that is missing is reported by fn.check_inputs instead.
if(dir.dbseabed == file.path(dir.gfisher,'data','dbseabed') &&
   !all(file.exists(fn.dbseabed_files(dir.dbseabed)))) fn.pull_dbseabed(dir.dbseabed)

inputs <- fn.check_inputs(fn.data_manifest(dir.data, file.gdb, file.spplist, file.depth,
                                           dir.dbseabed, file.seagrass, res))

if(!dir.exists(dir.maps)) dir.create(dir.maps, recursive=TRUE, showWarnings=FALSE)
depth <- raster::raster(file.depth)
if(interactive()) raster::plot(depth)

##geographic domain----
region <- 'WFS'
if(region=='GOM') bbox <- c(latN=30.5, latS=25, lonW=-97.5, lonE=-81)
if(region=='WFS') bbox <- c(latN=30.5, latS=25, lonW=-87.5, lonE=-81)

## SPECIES GROUPING SCHEME --------------------------------------------------------------------------
# Which columns of the spplist sheet in Master Species List.xlsx define the model groups that the MaxN
# heatmaps and habitat-affinity analysis are built for. This is the one knob to change groupings; add
# a new scheme by adding a row to group.schemes (spplist group cols + matching spp_stanzas_sizes fg col).
#   'ecospace' -> full Ecospace groups (modnum / modname / fg)
#   'mice'     -> MICE model groups    (modnum_mice / modname_mice / fg_mice)
# group.scheme itself is set in the SETTINGS block above (overridable in config.local.R).
group.schemes <- list(
  ecospace = c(modnum='modnum',      modname='modname',      fg='fg'),
  mice     = c(modnum='modnum_mice', modname='modname_mice', fg='fg_mice'))
if(!group.scheme %in% names(group.schemes)) stop("unknown group.scheme: ", group.scheme)
group.cols <- group.schemes[[group.scheme]]

#====================================================================================================
#STAGE 1 -- HABITAT BASEMAPS------------------------------------------------------------------------
# Nine layers summing to exactly 1 in every water cell: the six GFISHER reef classes, rock and
# unconsolidated bottom from raw dbSeabed, and seagrass. See the header of R/habitat_basemaps.R
# for what this changed relative to the retired cell-area maps (now in R/legacy/).
#
# dir.dbseabed and file.seagrass are set in the SETTINGS block above (overridable in config.local.R).
dir.basemaps <- file.path(dir.gfisher,'output','basemaps',paste0(res,'min'))
basemaps <- fn.make_habitat_basemaps(depth=depth, file.gdb=file.gdb, dir.raw=dir.dbseabed,
              dir.out=dir.basemaps,
              file.sgr=if(file.exists(file.seagrass)) file.seagrass else NULL,
              k='auto',            # per-class shrinkage constant from the variance structure
              target='stratum',    # borrow from the region x depth-bin mean; see header
              depth.max.reef=300)  # deepest cell in which GFISHER observed any reef class
fn.plot_habitat_basemaps(dir.basemaps)

# The retired cell-area maps, if you need to reproduce an older Ecospace run:
#   source(file.path('R','legacy','habitat_maps_cellarea.R'))
#   fn.make_GFISHER_habitat_maps(depth=depth, file.gdb=file.gdb, dir.maps=dir.maps)
#   fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,paste0(res,'min')))

#STAGE 2 -- PREPARE VIDEO DATASET-------------------------------------------------------------------
# Grouping columns come from the SPECIES GROUPING SCHEME block above, so the scheme is set in one place.
maxn <- fn.make_gfisher_videodataset(file.maxn, file.env, file.len, bbox, file.spplist,
                                     col.modnum=group.cols[['modnum']],
                                     col.modname=group.cols[['modname']],
                                     col.fg=group.cols[['fg']], seed=seed)

#STAGE 3 -- FISH MAXN HEATMAPS----------------------------------------------------------------------
# Outputs are scheme-tagged (.../maxn/<scheme>/) so different groupings coexist without clobbering.
graphics.off(); fn.plot_device()
dir.maxn <- file.path(dir.ewemaps,'GFISHER',paste0(res,'min'),'maxn',group.scheme)
maxn.stack <- fn.make_GFISHER_maxn_maps(maxn, depth, plot=T, fun=mean, background=0,
                                        dir.out=dir.maxn,
                                        save.format='all')        # one layer per model group

#STAGE 4a -- HABITAT AFFINITIES FROM SELECTION RATIOS (raster route)---------------------------------
# Sourcing only defines the functions (its own driver block is guarded), so we call the batch
# directly on the scheme's MaxN maps + the STAGE 1 basemaps built above.
source("R/selection_ratio_affinities.R")
dir.hab <- dir.basemaps                      # the sum-to-1 layers from stage 1
dir.aff <- file.path(dir.gfisher,'output',paste0('affinity_selratio_',group.scheme))
if(!dir.exists(dir.aff)) dir.create(dir.aff, recursive=TRUE)
# A is normalized WITHIN each layer family -- 6 reef (AL..NH) and 3 bottom (RCK/UNC/SGR) -- so
# A=1 marks the best reef layer AND, separately, the best bottom layer. The two are on separate
# scales and must not be compared across families: reef layers average 3e-5..1.3e-2 of a cell
# against UNC's 0.86, so a single global max would hand A=1 to a rare reef layer every time.
# Stanzas with too few MaxN>0 cells to fit (e.g. red-grouper-1) are pooled with their nearest
# same-species stanza; the long table's `pooled_with` column flags those rows.
#
# The MUD/SAND prior constraints that used to sit here are INERT against the basemaps, which
# carry unconsolidated bottom as one UNC layer rather than a GVL/SND/MUD split. The constraint
# machinery still exists (fn.constrain_w silently skips codes that are absent) -- it applies to
# LEGACY.SPEC runs. For the rock/gravel/sand/mud split, use STAGE 4c below, which reads raw
# dbSeabed at native resolution where that contrast actually exists.
affinity.constraints <- list(apply=FALSE)
hab <- fn.load_layer_stack(dir.hab)           # spec defaults to BASEMAP.SPEC
eff <- fn.build_effort_raster(file.env, hab[[1]],
         save.as=file.path(dir.aff, paste0('GFISHER_survey_effort_',res,'min_',nrow(hab),'x',ncol(hab),'.asc')))
aff <- fn.batch_selection_ratios(hab, dir.emp=dir.maxn, dir.out=dir.aff, effort=eff, n.boot=1000,
                                 con=affinity.constraints)
print(aff$wide)

#STAGE 4b -- HABITAT AFFINITIES FROM SITE-LEVEL PAIRED CAMERA + HABITAT DATA-------------------------
# Independent of the raster route above: works from the video stations themselves rather than from
# 5-min cells, in two families read AT the site --
#   reef      : HAB_STRAT, the side-scan habitat class, collapsed to the six Ecospace classes.
#   substrate : ARTI_PER/ROCK_PER/SED_PER percent cover scored from the video, as ART/RCK/SED.
#               The finer in-situ reads (shell-gravel, silt-sand-clay) exist on too few stations
#               to separate GVL/SND/MUD, so that split still has to come from dbSeabed.
# Because the survey targeted geoforms, availability is set by the sampler, so this estimates
# relative DENSITY (mean MaxN per station, zeros retained) per class rather than a use/availability
# ratio, normalized within each family. Stanzas seen at fewer than min.pos stations are pooled with
# their nearest same-species stanza (the long table's `pooled_with` column flags those rows).
# Habitat class is confounded with depth (mean station depth 75.8 m in NM vs 26.8 m in AM), so two
# adjustments are produced side by side; classes ranking the same way under both are the ones the
# habitat signal actually supports. See the header of R/site_level_affinities.R.
source("R/site_level_affinities.R")
dir.site <- file.path(dir.gfisher,'output',paste0('affinity_site_',group.scheme))
site <- fn.batch_site_affinities(maxn, file.env, bbox, dir.out=dir.site,
                                 file.gdb=file.gdb,     # stratum weights = mapped microgrid shares
                                 n.boot=300, tag=group.scheme)
print(site$wide.strat)   # A, stratified by SPACE_STRAT
print(site$wide.depth)   # A, depth as a continuous covariate





#STAGE 4c -- SUBSTRATE AFFINITIES FROM RAW dbSEABED-------------------------------------------------
# The rock/gravel/sand/mud split, estimated independently of both routes above: nothing conditions
# on reef type, and A is normalized only within the substrate families. Reads the RAW dbSeabed
# grids at their native 1.2 arc-min (dir.dbseabed, set in STAGE 1), because the processed 5-min
# layers renormalize rock against the grain-size triangle and treat NODATA as zero -- at 5 min only
# 13 stations in ONE cell sit on >=50% rock, versus 914 stations in 162 cells at native resolution.
#
# Two families: GVL/SND/MUD as a closed grain-size triangle, and RCK scored against its complement,
# because dbSeabed's rock is a separate areal fraction rather than a fourth member of the triangle.
# Requires R/site_level_affinities.R (sourced in STAGE 4b) for the shared station frame.
source("R/substrate_affinities.R")
dir.sub <- file.path(dir.gfisher,'output',paste0('affinity_substrate_',group.scheme))
sub <- fn.batch_substrate_affinities(maxn, file.env, bbox, dir.raw=dir.dbseabed, dir.out=dir.sub,
                                     file.gdb=file.gdb, n.boot=1000, tag=group.scheme)
print(sub$wide.strat)   # A, stratified by SPACE_STRAT
print(sub$wide.depth)   # A, standardized over depth bins
# NOTE: check sub$coverage before using the MUD column -- the survey avoided mud (targeting ~0.20),
# so it is flagged 'poor' and its ratio is not identifiable.
