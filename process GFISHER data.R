rm(list=ls());rm(.SavedPlots);graphics.off();gc();windows(record=T)
library('raster')
source("C:\\Users\\dchagaris\\Github\\WFS-FEM\\GFISHER\\R\\GFISHER habitat maps.R")

#setup----------------------------------------------------------------------------------------------
dir.gfisher <- getwd()
dir.maps <- "C:\\Users\\dchagaris\\OneDrive - University of Florida\\WFS Fisheries Ecosystem Modeling\\WFS EwE\\Ecospace\\maps\\GFISHER"
#dir.out <- "./maps/GFISHER"
file.gdb <- paste0(dir.gfisher,"/data/FWRI_East_Gulf_Mapping_2023.gdb")
file.gdb.micro <- paste0(dir.gfisher,"/data/East_Master_Hab_data_Dissolve_byMicro_13Sept24.gdb")

depth.1min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 1min 330x390.asc"))
depth.4min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 4min 82x97.asc"))
depth.6min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 6min 55x65.asc"))
depth.10min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 10min 33x39.asc"))
depth.5min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 5min 66x78.asc"))



fn.make_GFISHER_habitat_maps(depth=depth.5min, file.gdb=file.gdb, file.gdb.micro=file.gdb.micro, dir.maps=dir.maps)
fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,'5min'))

fn.make_GFISHER_habitat_maps(depth=depth.4min, file.gdb=file.gdb, file.gdb.micro=file.gdb.micro, dir.maps=dir.maps)
fn.plot_GFISHER_habitats(dir.maps=file.path(dir.maps,'4min'))

