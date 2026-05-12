rm(list=ls());rm(.SavedPlots);graphics.off();gc();windows(record=T)
library('raster')

#=========================== USER INPUTS ============================================================
# Set the working directory to the GFISHER repo root before running this script (setwd("path/to/GFISHER")).
# All other paths below resolve relative to it.

# file.gdb: ABSOLUTE path to your local copy of the GFISHER East Universe geodatabase.
# This file is NOT shipped with the repo (too large). Users running this code are expected to
# have their own copy of GFISHER_EAST_Universe_2026.gdb (or equivalent) and to point this at it.
file.gdb <- "C:/path/to/your/GFISHER_EAST_Universe_2026.gdb"

# res: map resolution in arc-minutes. 5 and 15 ship with the repo (see maps/bathymetry/).
res <- 5
#====================================================================================================

source(file.path('R','GFISHER functions.R'))

#setup----------------------------------------------------------------------------------------------
dir.gfisher <- getwd()
dir.maps <- file.path(dir.gfisher,'maps','GFISHER')
if(!dir.exists(dir.maps)) dir.create(dir.maps, recursive=TRUE, showWarnings=FALSE)

dir.data <- file.path(dir.gfisher,'data')
dir.scripts <- file.path(dir.gfisher,'R')
file.spplist <- file.path(dir.data,"Master Species List.xlsx")
file.sizeatage <- file.path(dir.data,'size_at_age.csv')

##geographic domain----
region <- 'WFS'
if(region=='GOM') bbox <- c(latN=30.5, latS=25, lonW=-97.5, lonE=-81)
if(region=='WFS') bbox <- c(latN=30.5, latS=25, lonW=-87.5, lonE=-81)

#HABITAT MAPS---------------------------------------------------------------------------------------
# file.gdb and res come from the USER INPUTS block at the top of this script.
file.depth <- list.files(file.path(dir.gfisher,'maps','bathymetry'), pattern=paste0(" ",res,"min"), full.names=T)
file.depth <- file.depth[grep("depth",basename(file.depth))]
file.depth <- file.depth[grep(".asc",basename(file.depth))]
depth <- raster(file.depth)

fn.make_GFISHER_habitat_maps(depth=depth, file.gdb=file.gdb, dir.maps=dir.maps)
fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,paste0(res,'min')))





