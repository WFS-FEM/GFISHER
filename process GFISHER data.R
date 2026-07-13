rm(list=ls());rm(.SavedPlots);graphics.off();gc();windows(record=T)
source(file.path('R','GFISHER functions.R'))
library('terra')

#=========================== USER INPUTS ============================================================
#setup----------------------------------------------------------------------------------------------
# Set the working directory to the GFISHER repo root before running this script (setwd("path/to/GFISHER")).
# All other paths below resolve relative to it.
dir.gfisher <- getwd()
dir.maps <- file.path(dir.gfisher,'output','maps')
if(!dir.exists(dir.maps)) dir.create(dir.maps, recursive=TRUE, showWarnings=FALSE)

# dir.ewemaps: where the MaxN heatmap outputs are written. Defaults to the in-repo output/maps/ so
# the script runs out of the box. To write them into an external Ecospace maps tree instead, set
# dir.ewemaps.ext to its path; it is used only when it exists, otherwise the in-repo output/maps/ is used.
dir.ewemaps <- dir.maps
dir.ewemaps.ext <- ""   # e.g. "C:/Users/<you>/OneDrive .../WFS EwE/Ecospace/maps"
if(nzchar(dir.ewemaps.ext) && dir.exists(dir.ewemaps.ext)) dir.ewemaps <- dir.ewemaps.ext

# dir.ecospace.maps: external Ecospace maps tree that holds the 10 habitat/sediment layers the
# habitat-affinity step reads (input_ascii_sum1/). Read-only; not produced by this repo.
dir.ecospace.maps <- "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps"

# file.gdb: ABSOLUTE path to your local copy of the GFISHER East Universe geodatabase.
# This file is NOT shipped with the repo (too large). Users running this code are expected to
# have their own copy of GFISHER_EAST_Universe_2026.gdb (or equivalent) and to point this at it.
dir.data <- file.path(dir.gfisher,'data','April2026')
dir.scripts <- file.path(dir.gfisher,'R')
file.gdb <- file.path(dir.data,"GFISHER_EAST_Universe_2026.gdb")
file.spplist <- file.path(dirname(dir.data),"Master Species List.xlsx")
file.sizeatage <- file.path(dirname(dir.data),'size_at_age.csv')
file.maxn = file.path(dir.data,'maxn3LABS_93to24.csv')
file.env = file.path(dir.data,'env3LABS_93to24.csv')
file.len = file.path(dir.data,'lens3LABS_93to24.csv')
dir.bathy <- file.path(dir.gfisher,'data','bathymetry')
dir.bathy.ext <- ""   # e.g. "C:/Users/<you>/OneDrive .../WFS EwE/Ecospace/maps/bathymetry"
if(nzchar(dir.bathy.ext) && dir.exists(dir.bathy.ext)) dir.bathy <- dir.bathy.ext

# res: map resolution in arc-minutes. 5 and 15 ship with the repo (see data/bathymetry/).
res <- 5
file.depth <- list.files(dir.bathy,pattern=paste0('depth ',res,'min'),full.names=TRUE)
if(length(file.depth)==0) stop(paste0("No 'depth ",res,"min' raster found in ",dir.bathy))
depth <- rast(file.depth)
plot(depth)

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
group.scheme <- 'mice'
group.schemes <- list(
  ecospace = c(modnum='modnum',      modname='modname',      fg='fg'),
  mice     = c(modnum='modnum_mice', modname='modname_mice', fg='fg_mice'))
if(!group.scheme %in% names(group.schemes)) stop("unknown group.scheme: ", group.scheme)
group.cols <- group.schemes[[group.scheme]]

#====================================================================================================
#MAKE HABITAT MAPS---------------------------------------------------------------------------------------
# file.gdb and res come from the USER INPUTS block at the top of this script.
fn.make_GFISHER_habitat_maps(depth=depth, file.gdb=file.gdb, dir.maps=dir.maps)
fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,paste0(res,'min')))

#PREPARE VIDEO DATASET------------------------------------------------------------------------------
# Grouping columns come from the SPECIES GROUPING SCHEME block above, so the scheme is set in one place.
maxn <- fn.make_gfisher_videodataset(file.maxn, file.env, file.len, bbox, file.spplist,
                                     col.modnum=group.cols[['modnum']],
                                     col.modname=group.cols[['modname']],
                                     col.fg=group.cols[['fg']])

#FISH MAXN HEATMAPS-----------------------------------------------------------------------------------
# Outputs are scheme-tagged (.../maxn/<scheme>/) so different groupings coexist without clobbering.
class(maxn)
graphics.off();rm(.SavedPlots);windows(record=T)
dir.maxn <- file.path(dir.ewemaps,'GFISHER',paste0(res,'min'),'maxn',group.scheme)
maxn.stack <- fn.make_GFISHER_maxn_maps(maxn, depth, plot=T, fun=mean, background=0,
                                        dir.out=dir.maxn,
                                        save.format='all')        # one layer per model group

#HABITAT AFFINITIES FROM SELECTION RATIOS-----------------------------------------------------------
# Sourcing only defines the functions (its own driver block is guarded), so we call the batch directly
# on the scheme's MaxN maps + the 10 habitat/sediment layers from the external Ecospace maps tree.
source("R/selection_ratio_affinities.R")
dir.hab <- file.path(dir.ecospace.maps,'input_ascii_sum1',paste0(res,'min'))
dir.aff <- file.path(dir.gfisher,'output',paste0('affinity_selratio_',group.scheme))
if(!dir.exists(dir.aff)) dir.create(dir.aff, recursive=TRUE)
# Prior habitat constraints imposed on the empirical selection ratios before forming affinities:
#   - MUD forced to 0 (not identifiable from the reef-targeted video survey)
#   - SAND forced strictly below both ROCK and GRAVEL (hard substrate preferred), 5% margin.
# Set apply=FALSE for the raw, unconstrained affinities.
affinity.constraints <- list(apply=TRUE, zero.codes='MUD', sand.code='SND',
                             hard.codes=c('RCK','GVL'), sand.margin=0.05)
hab <- fn.load_layer_stack(dir.hab)
eff <- fn.build_effort_raster(file.env, hab[[1]],
         save.as=file.path(dir.aff, paste0('GFISHER_survey_effort_',res,'min_',nrow(hab),'x',ncol(hab),'.asc')))
aff <- fn.batch_selection_ratios(hab, dir.emp=dir.maxn, dir.out=dir.aff, effort=eff, n.boot=1000,
                                 con=affinity.constraints)
print(aff$wide)




