rm(list=ls());rm(.SavedPlots);graphics.off();gc();windows(record=T)
source(file.path('R','GFISHER functions.R'))
library('terra')

#=========================== USER INPUTS ============================================================
#setup----------------------------------------------------------------------------------------------
# Set the working directory to the GFISHER repo root before running this script (setwd("path/to/GFISHER")).
# All other paths below resolve relative to it.
dir.gfisher <- getwd()
dir.maps <- file.path(dir.gfisher,'maps','GFISHER')
if(!dir.exists(dir.maps)) dir.create(dir.maps, recursive=TRUE, showWarnings=FALSE)

# file.gdb: ABSOLUTE path to your local copy of the GFISHER East Universe geodatabase.
# This file is NOT shipped with the repo (too large). Users running this code are expected to
# have their own copy of GFISHER_EAST_Universe_2026.gdb (or equivalent) and to point this at it.
dir.data <- file.path(dir.gfisher,'data')
dir.scripts <- file.path(dir.gfisher,'R')
file.gdb <- "C:/path/to/your/GFISHER_EAST_Universe_2026.gdb"
file.spplist <- file.path(dir.data,"Master Species List.xlsx")
file.sizeatage <- file.path(dir.data,'size_at_age.csv')
file.maxn = file.path(dir.data,'April2026','maxn3LABS_93to24.csv')
file.env = file.path(dir.data,'April2026','env3LABS_93to24.csv')
file.len = file.path(dir.data,'April2026','lens3LABS_93to24.csv')
dir.bathy <- file.path(dir.gfisher,'maps','bathymetry')
dir.bathy.ext <- ""   # e.g. "C:/Users/<you>/OneDrive .../WFS EwE/Ecospace/maps/bathymetry"
if(nzchar(dir.bathy.ext) && dir.exists(dir.bathy.ext)) dir.bathy <- dir.bathy.ext

# res: map resolution in arc-minutes. 5 and 15 ship with the repo (see maps/bathymetry/).
res <- 5
file.depth <- list.files(dir.bathy,pattern=paste0('depth ',res,'min'),full.names=TRUE)
if(length(file.depth)==0) stop(paste0("No 'depth ",res,"min' raster found in ",dir.bathy))
depth <- rast(file.depth)
plot(depth)

##geographic domain----
region <- 'WFS'
if(region=='GOM') bbox <- c(latN=30.5, latS=25, lonW=-97.5, lonE=-81)
if(region=='WFS') bbox <- c(latN=30.5, latS=25, lonW=-87.5, lonE=-81)

#====================================================================================================
#MAKE HABITAT MAPS---------------------------------------------------------------------------------------
# file.gdb and res come from the USER INPUTS block at the top of this script.
fn.make_GFISHER_habitat_maps(depth=depth, file.gdb=file.gdb, dir.maps=dir.maps)
fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,paste0(res,'min')))

#PREPARE VIDEO DATASET------------------------------------------------------------------------------
maxn <- fn.make_gfisher_videodataset(file.maxn, file.env, file.len, bbox,spplist)

#FISH MAXN HEATMAPS-----------------------------------------------------------------------------------
class(maxn)
graphics.off();rm(.SavedPlots);windows(record=T)
maxn.stack <- fn.make_GFISHER_maxn_maps(maxn, depth, plot=T, fun=mean, background=0, dir.out="maps/GFISHER/maxn")        # one layer per model group
names(maxn.stack)
par(mfrow=c(3,3))
for(i in 1:nlayers(maxn.stack)){
  plot(maxn.stack[[i]],colNA='black', main=names(maxn.stack)[i])
}
