rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
.libPaths("C:\\R\\win-library")
library('sf')
library('sp')
library('raster')
library('rgdal')
library('colorRamps')
library('maps')

#setup----------------------------------------------------------------------------------------------
dir.gfisher <- getwd()
dir.maps <- "C:\\Users\\dchagaris\\OneDrive - University of Florida\\WFS Fisheries Ecosystem Modeling\\WFS EwE\\Ecospace\\maps\\GFISHER"
#dir.out <- "./maps/GFISHER"
file.gdb <- paste0(dir.gfisher,"/FWRI_East_Gulf_Mapping_2023.gdb")
file.gdb.micro <- paste0(dir.gfisher,"/East_Master_Hab_data_Dissolve_byMicro_13Sept24.gdb")

depth.1min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 1min 330x390.asc"))
depth.6min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 6min 55x65.asc"))
depth.10min <- raster(paste0(dirname(dir.maps),"/bathymetry/depth 10min 33x39.asc"))

#import and prepare geodatabase file------------------------------------------------------------------------------
# The directory contains a zip file with a geodatabase and supporting metadata to describe the datasets. Two shapefiles are in the geodatabase:
# 1.	East_Master_microgrid_mapped – this shapefile consists of microgrids (0.1 x 0.1 nm) which we’ve mapped with side scan sonar.  Each grid represents one of our primary sampling units that has been scanned (mapping footprint).
# 2.	FWRI_East_Master_Hab_Data – this shapefile consists of our digitized polygons for our habitat classes (i.e. Geoforms).  The habitat classes are described in the metadata file and consist of both natural and artificial reefs.
st_layers(file.gdb)
st_layers(file.gdb.micro)

#extract shapefile from gdb
print("Reading Geodatabase Files")
microgrid <- st_read(file.gdb, layer="East_Master_Microgrid_Mapped_2023")   
#habitat <- st_read(file.gdb, layer="East_Master_Hab_Data_FINAL_2023")
habitat2 <- st_read(file.gdb.micro, layer="East_Master_Hab_Data_Dissolve_Site_Selection_2023")

#one record is a MULTISURFACE class, instead of MULTIPOLYGON.  Don't know how to deal so remove it.
habitat2.shapes = as.data.frame(matrix(unlist(lapply(habitat2$Shape,FUN=function(x)class(x))),ncol=3,byrow=T))
habitat2 = habitat2[-which(habitat2.shapes$V2=='MULTISURFACE'),]

#convert to sp, spatialpolygonsdataframe, there are many many polygons so everything is slow until it
#gets into a raster format
print("Converting to sp...")
microgrid.sp <- as_Spatial(microgrid)
#habitat.sp <- as_Spatial(habitat)
habitat2.sp <- as_Spatial(habitat2)

#get depth and habitat on same CRS
print("Transforming CRS to match depth grid...")
crs(habitat.sp); crs(habitat2.sp); crs(depth)
#habitat.sp <- spTransform(habitat.sp,crs(depth))
habitat2.sp <- spTransform(habitat2.sp,crs(depth))
microgrid.sp <- spTransform(microgrid.sp, crs(depth))
#habitat.sp <- habitat.sp[habitat.sp$MicroGrid!=" ",]


#convert polygons to points-------------------------------------------------------------------------
#footprint polys to points - since the polygons are small, we can take the centroid of each poly and
#convert to spatial points dataframe, then average centroids for each grid cell
micro.centroids <- coordinates(microgrid.sp)
microgrid.pts <- microgrid.sp@data
microgrid.pts$cent.long <- micro.centroids[,1]
microgrid.pts$cent.lat <- micro.centroids[,2]
coordinates(microgrid.pts) <- ~cent.long+cent.lat

# micro2.centroids <- coordinates(microgrid2.sp)
# microgrid2.pts <- microgrid2.sp@data
# microgrid2.pts$cent.long <- micro2.centroids[,1]
# microgrid2.pts$cent.lat <- micro2.centroids[,2]
# coordinates(microgrid2.pts) <- ~cent.long+cent.lat

#habitat polys to points
# tapply(habitat.sp$Shape_Area,habitat.sp$NewHab,sum)
# hab.centroids <- coordinates(habitat.sp)
# habitat.pts <- habitat.sp@data
# habitat.pts$cent.long <- hab.centroids[,1]
# habitat.pts$cent.lat <- hab.centroids[,2]
# coordinates(habitat.pts) <- ~cent.long+cent.lat

tapply(habitat2.sp$Shape_Area,habitat2.sp$NewHabStrat,sum)
habitat2.sp$NewHab <- substr(habitat2.sp$NewHabStrat,1,2)  #paste0(substr(habitat2.sp$NewHabStrat,1,1),substr(habitat2.sp$NewHabStrat,3,3))
tapply(habitat2.sp$Shape_Area,habitat2.sp$NewHab,sum)
hab2.centroids <- coordinates(habitat2.sp)
habitat2.pts <- habitat2.sp@data
habitat2.pts$cent.long <- hab2.centroids[,1]
habitat2.pts$cent.lat <- hab2.centroids[,2]
coordinates(habitat2.pts) <- ~cent.long+cent.lat

#get microgrid for reach hab centroid
habitat2.pts.sf <- st_as_sf(habitat2.pts)
microgrid.sf <- st_as_sf(microgrid.sp)
st_crs(habitat2.pts.sf) <- st_crs(microgrid.sf)
habitat2.int <- st_intersects(habitat2.pts.sf$geometry,microgrid.sf$geometry)
npolys <- unlist(lapply(habitat2.int,length))
table(npolys)
microgrid.id <- unlist(habitat2.int)
habitat2$MicroGrid <- microgrid$MicroGrid[microgrid.id]


#checks---------------------------------------------------------------------------------------------
hab <- habitat2
mcg <- microgrid
hab$Shape<-NULL;mcg$Shape<-NULL
names(hab)[which(names(hab)=='Shape_Area')] <- 'Hab_Area'
names(mcg)[which(names(mcg)=='Shape_Area')] <- 'Grid_Area'

#do all the records match?
length(unique(hab$MicroGrid))
length(unique(mcg$MicroGrid))
length(which(unique(hab$MicroGrid) %in% mcg$MicroGrid))
length(which(!hab$MicroGrid %in% mcg$MicroGrid))
print(paste("N habitat records with matching microgrids:",length(which(hab$MicroGrid %in% mcg$MicroGrid))))
print(paste("N habitat records with non-matching microgrid:",nrow(hab)-length(which(hab$MicroGrid %in% mcg$MicroGrid))-length(which(hab$MicroGrid==" "))))
print(paste("N habitat records with missing microgrid:",length(which(hab$MicroGrid==" "))))
print(paste("N microgrids without any habitat records:",length(which(!mcg$MicroGrid %in% hab$MicroGrid))))

#does habitat area exceed mapped area for any microgrids?
hab.sum <- aggregate(Hab_Area~MicroGrid, data=hab, sum)
chk1 <- merge(hab.sum,mcg[,c('MicroGrid','Grid_Area')],by='MicroGrid',all=T)
chk1$habpct <- round(chk1$Hab_Area/chk1$Grid_Area,4)
#print(paste("N microgrids with Hab_Area>Grid_Area:",length(which(round(chk1$Hab_Area,0)>round(chk1$Grid_Area,0)))))
print(paste("N microgrids with Hab_Area>Grid_Area:",length(which(chk1$habpct>1))))
chk2 <- chk1[which(chk1$habpct>1),]
if(nrow(chk2)>0) write.csv(chk2,'hab area vs grid area.csv',row.names=F)

#rasterize and make maps for ecospace---------------------------------------------------------------
#rasterize mapping footprint
tiff(filename=paste0(dir.gfisher,"/Mapping footprint microgrids.tiff"),height=7,width=7,units='in',res=600,compression='lzw')
plot(microgrid$Shape,main="Mapping Footprint (microgrids)"); map(database='state',fill=T,add=T,col='lightgray')
dev.off()

microgrid.ras1 = rasterize(microgrid.pts, depth.1min, field='Shape_Area', fun='sum', background=0)
microgrid.ras6 <- rasterize(microgrid.pts, depth.6min, field='Shape_Area', fun='sum', background=0)#resample(microgrid.ras, depth10, method='bilinear')
microgrid.ras10 <- rasterize(microgrid.pts, depth.10min, field='Shape_Area', fun='sum', background=0)#resample(microgrid.ras, depth10, method='bilinear')
microgrid.ras1n = rasterize(microgrid.pts, depth.1min, field='Shape_Area', fun='count', background=0)
microgrid.ras6n <- rasterize(microgrid.pts, depth.6min, field='Shape_Area', fun='count', background=0)#resample(microgrid.ras, depth10, method='bilinear')
microgrid.ras10n <- rasterize(microgrid.pts, depth.10min, field='Shape_Area', fun='count', background=0)#resample(microgrid.ras, depth10, method='bilinear')

microgrid.ras1[is.na(depth.1min)] <- NA
microgrid.ras6[is.na(depth.6min)] <- NA
microgrid.ras10[is.na(depth.10min)] <- NA

tiff(filename=paste0(dir.gfisher,"/Mapping footprint plots.tiff"),height=7,width=7,units='in',res=600,compression='lzw')
par(mfcol=c(3,2),mar=c(1,2,4,4))
plot(microgrid.ras1,main="Area mapped (m2)\nper 1-min grid cell"); map(database='state',fill=T,add=T,col='lightgray')
plot(microgrid.ras1n,main="N microgrids\nper 1-min grid cell (m2)"); map(database='state',fill=T,add=T,col='lightgray')
plot(microgrid.ras6,main="Area mapped (m2)\nper 6-min grid cell"); map(database='state',fill=T,add=T,col='lightgray')
plot(microgrid.ras6n,main="N microgrids\nper 6-min grid cell"); map(database='state',fill=T,add=T,col='lightgray')
plot(microgrid.ras10,main="Area mapped (m2)\nper 10-min grid cell"); map(database='state',fill=T,add=T,col='lightgray')
plot(microgrid.ras10n,main="N microgrids\nper 10-min grid cell"); map(database='state',fill=T,add=T,col='lightgray')
dev.off()

#rasterize habitats
# newhabs = sort(unique(habitat.pts$NewHab))
# newhabs = newhabs[c(2,3,1,5,6,4)]
# habpct.stack.1min <- habpct.stack.6min <- stack()
# for(i in 1:length(newhabs)){
# #i=1
# #1-min res
# hab.ras = rasterize(habitat.pts[habitat.pts$NewHab==newhabs[i],],depth,field='Shape_Area',fun='sum', background=0)
# hab.ras[is.na(depth)] <- NA
# plot(hab.ras)
# range(getValues(hab.ras),na.rm=T)
# 
# #6-min res
# hab.ras2 = rasterize(habitat.pts[habitat.pts$NewHab==newhabs[i],],depth10,field='Shape_Area',fun='sum', background=0)
# hab.ras2[is.na(depth10)] <- NA
# plot(hab.ras2)
# range(getValues(hab.ras2),na.rm=T)
# 
# #convert to percent of area mapped
# habpct.ras <- 100*hab.ras/microgrid.ras
# habpct.ras2 <- 100*hab.ras2/microgrid.ras2
# habpct.ras[is.na(habpct.ras)] <- 0
# habpct.ras2[is.na(habpct.ras2)] <- 0
# habpct.ras[is.na(depth)] <- NA
# habpct.ras2[is.na(depth10)] <- NA
# #plot(habpct.ras)
# #plot(habpct.ras2)
# 
# #cleanup
# habpct.stack.1min <- addLayer(habpct.stack.1min, habpct.ras)
# habpct.stack.6min <- addLayer(habpct.stack.6min,habpct.ras2)
# rm(hab.ras, hab.ras2, habpct.ras, habpct.ras2); gc()
# }
# names(habpct.stack.1min) <- names(habpct.stack.6min) <- newhabs
# which(habpct.stack.1min[[4]]>100)

#for new dataset
newhabs = sort(unique(habitat2.pts$NewHab))
newhabs = newhabs[c(2,3,1,5,6,4)]
habpct2.stack1 <- habpct2.stack6 <- habpct2.stack10 <- stack()
for(i in 1:length(newhabs)){
  #i=1
  
  hab.ras1 = rasterize(habitat2.pts[habitat2.pts$NewHab==newhabs[i],],depth.1min,field='Shape_Area',fun='sum', background=0)
  hab.ras6 = rasterize(habitat2.pts[habitat2.pts$NewHab==newhabs[i],],depth.6min,field='Shape_Area',fun='sum', background=0)
  hab.ras10 = rasterize(habitat2.pts[habitat2.pts$NewHab==newhabs[i],],depth.10min,field='Shape_Area',fun='sum', background=0)
  
  hab.ras1[is.na(depth.1min)] <- NA
  hab.ras6[is.na(depth.6min)] <- NA
  hab.ras10[is.na(depth.10min)] <- NA

  #convert to percent of area mapped
  habpct.ras1 <- 100*hab.ras1/microgrid.ras1
  habpct.ras6 <- 100*hab.ras6/microgrid.ras6
  habpct.ras10 <- 100*hab.ras10/microgrid.ras10
  
  habpct.ras1[is.na(habpct.ras1)] <- 0
  habpct.ras6[is.na(habpct.ras6)] <- 0
  habpct.ras10[is.na(habpct.ras10)] <- 0
  
  #stack
  habpct2.stack1 <- addLayer(habpct2.stack1, habpct.ras1)
  habpct2.stack6 <- addLayer(habpct2.stack6, habpct.ras6)
  habpct2.stack10 <- addLayer(habpct2.stack10, habpct.ras10)
  
  #cleanup
  rm(hab.ras1, hab.ras6, hab.ras10, habpct.ras1, habpct.ras6, habpct.ras10); gc()
}
names(habpct2.stack1) <- names(habpct2.stack6) <- names(habpct2.stack10) <- newhabs




tiff(filename=paste0(dir.gfisher,"/Habitat pct area 1min.tiff"),height=7,width=7,units='in',res=600,compression='lzw')
par(mfcol=c(3,2),mar=c(3,3,3,5))
for(i in 1:nlayers(habpct2.stack1)){
  #i=1
  hab.i = habpct2.stack1[[i]]
  #brks = unique(quantile(hab.ras2,c(0,seq(.8,1,.01))))
  #cols = matlab.like(n=length(brks)-1)
  #plot(hab.ras2, colNA='lightgray', main=paste(newhabs[i], 'sum'))#, col=cols, breaks=brks)
  #brks = unique(round(c(quantile(habpct.ras2,c(0,seq(.8,1,.01))),max(habpct.ras2[habpct.ras2<Inf],na.rm=T)),4))

  brks = pretty(getValues(hab.i),n=50)
  if(length(which(brks==Inf))>0) brks = brks[-which(brks==Inf)]
  colv = c(colorRamps::matlab.like2(n=length(brks)))
  col.bias=5
  #if(i==1) col.bias=3
  funpal  = colorRampPalette(colv,bias=col.bias,interpolate='spline')
  cols   = funpal(length(brks)-1)
  plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'), col=cols, breaks=brks,legend=T)
  #plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'))
  map(database='state',add=T, fill=T, col='lightgray')
}
dev.off()

#----------------------------------------------------
tiff(filename=paste0(dir.gfisher,"/Habitat pct area 6min.tiff"),height=7,width=7,units='in',res=600,compression='lzw')
par(mfcol=c(3,2),mar=c(3,3,3,5))
for(i in 1:nlayers(habpct2.stack6)){
  #i=1
  hab.i = habpct2.stack6[[i]]
  #brks = unique(quantile(hab.ras2,c(0,seq(.8,1,.01))))
  #cols = matlab.like(n=length(brks)-1)
  #plot(hab.ras2, colNA='lightgray', main=paste(newhabs[i], 'sum'))#, col=cols, breaks=brks)
  #brks = unique(round(c(quantile(habpct.ras2,c(0,seq(.8,1,.01))),max(habpct.ras2[habpct.ras2<Inf],na.rm=T)),4))
  
  brks = pretty(getValues(hab.i),n=50)
  if(length(which(brks==Inf))>0) brks = brks[-which(brks==Inf)]
  colv = c(colorRamps::matlab.like2(n=length(brks)))
  col.bias=5
  #if(i==1) col.bias=3
  funpal  = colorRampPalette(colv,bias=col.bias,interpolate='spline')
  cols   = funpal(length(brks)-1)
  plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'), col=cols, breaks=brks,legend=T)
  #plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'))
  map(database='state',add=T, fill=T, col='lightgray')
}
dev.off()

#----------------------------------------------------
tiff(filename=paste0(dir.gfisher,"/Habitat pct area 10min.tiff"),height=7,width=7,units='in',res=600,compression='lzw')
par(mfcol=c(3,2),mar=c(3,3,3,5))
for(i in 1:nlayers(habpct2.stack10)){
  #i=1
  hab.i = habpct2.stack10[[i]]
  #brks = unique(quantile(hab.ras2,c(0,seq(.8,1,.01))))
  #cols = matlab.like(n=length(brks)-1)
  #plot(hab.ras2, colNA='lightgray', main=paste(newhabs[i], 'sum'))#, col=cols, breaks=brks)
  #brks = unique(round(c(quantile(habpct.ras2,c(0,seq(.8,1,.01))),max(habpct.ras2[habpct.ras2<Inf],na.rm=T)),4))
  
  brks = pretty(getValues(hab.i),n=50)
  if(length(which(brks==Inf))>0) brks = brks[-which(brks==Inf)]
  colv = c(colorRamps::matlab.like2(n=length(brks)))
  col.bias=5
  #if(i==1) col.bias=3
  funpal  = colorRampPalette(colv,bias=col.bias,interpolate='spline')
  cols   = funpal(length(brks)-1)
  plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'), col=cols, breaks=brks,legend=T)
  #plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'))
  map(database='state',add=T, fill=T, col='lightgray')
}
dev.off()



#output---------------------------------------------------------------------------------------------
writeRaster(habpct2.stack10,filename=paste0(dir.maps,"/GFISHER"),bylayer=T,
            suffix=paste0(names(habpct2.stack10),"_pct_10min_",dim(habpct2.stack10)[1],"x",dim(habpct2.stack10)[2],"_.asc"))
#writeRaster(habitat.stack.sum,filename=paste0(dir.maps,"/GFISHER"),bylayer=T,
#            suffix=paste0(names(habitat.stack.sum),"_sum_",dim(habitat.stack.sum)[1],"x",dim(habitat.stack.sum)[2],"_.asc"))


