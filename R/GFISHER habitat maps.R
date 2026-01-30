library('sf')
library('sp')
#library('rgdal')
library('colorRamps')
library('maps')
library('lwgeom')
library('gstat')

fn.make_GFISHER_habitat_maps <- function(file.gdb, file.gdb.micro, dir.maps, depth=depth.10min){

#import and prepare geodatabase file--------------------------------------------
# The directory contains a zip file with a geodatabase and supporting metadata to describe the datasets. Two shapefiles are in the geodatabase:
# 1.	East_Master_microgrid_mapped – this shapefile consists of microgrids (0.1 x 0.1 nm) which we’ve mapped with side scan sonar.  Each grid represents one of our primary sampling units that has been scanned (mapping footprint).
# 2.	FWRI_East_Master_Hab_Data – this shapefile consists of our digitized polygons for our habitat classes (i.e. Geoforms).  The habitat classes are described in the metadata file and consist of both natural and artificial reefs.
res.min = round(res(depth)[1]*60,0)
lyrs.gdb = st_layers(file.gdb)
lyrs.gdb.micro = st_layers(file.gdb.micro)

#extract shapefile from gdb
cat("Reading Geodatabase Files...\n")
microgrid <- st_read(file.gdb, layer=lyrs.gdb$name[grep("Microgrid",lyrs.gdb$name)])   
habitat <- st_read(file.gdb.micro, layer=lyrs.gdb.micro$name[grep("Data",lyrs.gdb.micro$name)])
#habitat <- st_make_valid(habitat)

#one record is a MULTISURFACE class, instead of MULTIPOLYGON.  Don't know how to deal so remove it.
habitat.shapes <- as.data.frame(matrix(unlist(lapply(habitat$Shape,FUN=function(x)class(x))),ncol=3,byrow=T))
if(length(which(habitat.shapes$V2=='MULTISURFACE'))>0) habitat <- habitat[-which(habitat.shapes$V2=='MULTISURFACE'),]

#the next line defines the habitat type
habitat$NewHab <- substr(habitat$NewHabStrat,1,2)  #paste0(substr(habitat2.sp$NewHabStrat,1,1),substr(habitat2.sp$NewHabStrat,3,3))

#convert to sp, spatialpolygonsdataframe, there are many many polygons so everything is slow until it
#gets into a raster format
cat("Converting to spatialpolygonsdataframe...\n")
microgrid.sp <- as_Spatial(microgrid)
habitat.sp <- as_Spatial(habitat)

#get depth and habitat on same CRS
cat("Transforming CRS to match depth grid...\n")
microgrid.sp <- spTransform(microgrid.sp, crs(depth))
habitat.sp <- spTransform(habitat.sp,crs(depth))

#convert polygons to points-----------------------------------------------------
#footprint polys to points - since the polygons are small, we can take the centroid of each poly and
#convert to spatial points dataframe, then average or sum centroids for each grid cell
#microgrid centroids are already stored in the sp object
microgrid.pts <- microgrid.sp@data           #make a new points object from the data
coordinates(microgrid.pts) <- ~X+Y
# micro.centroids <- coordinates(microgrid.sp)  #get centroids of polygons
# microgrid.pts$cent.long <- micro.centroids[,1] #bring in centroids long
# microgrid.pts$cent.lat <- micro.centroids[,2] #bring in centroids lat

#habitat polys to points  currently this takes the center of all habitat polygons in multipolygon feature.
#alternatively, could cast them to individual parts, get centroid, and calculate new area
hab.centroids <- coordinates(habitat.sp)
habitat.pts <- habitat.sp@data
habitat.pts$cent.long <- hab.centroids[,1]
habitat.pts$cent.lat <- hab.centroids[,2]
coordinates(habitat.pts) <- ~cent.long+cent.lat

#get microgrid for each hab centroid, must convert to sf objects
# habitat.pts.sf <- st_as_sf(habitat.pts)
# microgrid.sf <- st_as_sf(microgrid.sp)
# st_crs(habitat.pts.sf) <- st_crs(microgrid)
# habitat.int <- st_intersects(habitat.pts.sf$geometry,microgrid.sf$geometry)
# npolys <- unlist(lapply(habitat.int,length))
# microgrid.id <- unlist(habitat.int)
# habitat$MicroGrid <- microgrid$MicroGrid[microgrid.id]

hab.centroids.sf <- st_centroid(habitat)
habitat.int <- st_intersects(hab.centroids.sf$Shape, microgrid$Shape)
npolys <- unlist(lapply(habitat.int,length))
microgrid.id <- unlist(habitat.int)
habitat$MicroGrid <- microgrid$MicroGrid[microgrid.id]

#checks-------------------------------------------------------------------------
cat("Check\n")
hab <- habitat
mcg <- microgrid
hab$Shape<-NULL;mcg$Shape<-NULL
names(hab)[which(names(hab)=='Shape_Area')] <- 'Hab_Area'
names(mcg)[which(names(mcg)=='Shape_Area')] <- 'Grid_Area'

#do all the records match?
length(unique(hab$MicroGrid))
length(unique(mcg$MicroGrid))
length(which(unique(hab$MicroGrid) %in% mcg$MicroGrid))
length(which(!hab$MicroGrid %in% mcg$MicroGrid))
cat(paste("N habitat records with matching microgrids:",length(which(hab$MicroGrid %in% mcg$MicroGrid)),"\n"))
cat(paste("N habitat records with non-matching microgrid:",nrow(hab)-length(which(hab$MicroGrid %in% mcg$MicroGrid))-length(which(hab$MicroGrid==" ")),"\n"))
cat(paste("N habitat records with missing microgrid:",length(which(hab$MicroGrid==" ")),"\n"))
cat(paste("N microgrids without any habitat records:",length(which(!mcg$MicroGrid %in% hab$MicroGrid)),"\n"))

#does habitat area exceed mapped area for any microgrids?
hab.sum <- aggregate(Hab_Area~MicroGrid, data=hab, sum)
chk1 <- merge(hab.sum,mcg[,c('MicroGrid','Grid_Area')],by='MicroGrid',all=T)
chk1$habpct <- round(chk1$Hab_Area/chk1$Grid_Area,4)
cat(paste("N microgrids with Hab_Area>Grid_Area:",length(which(chk1$habpct>1)),"\n"))
chk2 <- chk1[which(chk1$habpct>1),]
if(nrow(chk2)>0) write.csv(chk2,file.path(dir.maps,'hab area vs grid area.csv'),row.names=F)

#rasterize----------------------------------------------------------------------
##helper functions----
###function to extrapolate data to missing cells
idw_fill_raster_longlat <- function(r, idp = 2, nmax = 12,
                                    ea_crs = "+proj=aea +lat_1=24 +lat_2=31.5 +lat_0=23 +lon_0=-84 +datum=WGS84 +units=m +no_defs",
                                    depth = depth,
                                    mask = NULL) {
  #r=habpct.ras
  if (is.null(crs(r))) stop("Input raster must have a valid CRS.")
  
  # 1) Project raster to equal-area (meters)
  r_ea <- projectRaster(r, crs = ea_crs, method = "bilinear")
  depth_ea <- projectRaster(depth, crs = ea_crs, method = "bilinear") 

  # 2) Build eligibility mask on EA grid (where we will predict)
  #    If no mask supplied, eligible = NA cells in r_ea
  eligible_ea <- r_ea
  eligible_ea[] <- NA
  
  # Map original NA indices to EA grid via projected raster
  na_idx_ea <- which(is.na(r_ea[]) & !is.na(depth_ea[]) & depth_ea[]<=500)
  eligible_ea[na_idx_ea] <- 1

  # If a mask (on original grid) was provided, project it and combine
  if (!is.null(mask)) {
    mask_ea <- projectRaster(mask, r_ea, method = "ngb")  # categorical mask -> nearest neighbor
    eligible_ea[is.na(mask_ea[])] <- NA  # block outside mask
  }
  
  # If nothing to fill, return original
  if (length(na_idx_ea) == 0) {
    message("No NA cells to fill.")
    return(r)
  }
  
  # 3) Known points = all non-NA cells in r_ea (include zeros)
  known_sp <- rasterToPoints(r_ea, spatial = TRUE)
  names(known_sp) <- "z"
  known_sp <- known_sp[!is.na(known_sp$z), ]  # keep only known
  
  # 4) Prediction points = eligible NA cells
  pred_sp <- rasterToPoints(eligible_ea, spatial = TRUE)
  # keep only marked '1' (eligible)
  pred_sp <- pred_sp[!is.na(pred_sp$layer) & pred_sp$layer == 1, ]
  pred_sp$layer <- NULL
  
  if (nrow(pred_sp) == 0) {
    warning("No eligible cells for IDW. Returning original raster.")
    return(r)
  }
  
  # 5) IDW prediction (local): idp = power; nmax = nearest neighbors
  idw_pred <- gstat::idw(z ~ 1, locations = known_sp, newdata = pred_sp,
                         idp = idp, nmax = nmax)
  
  # 6) Rasterize predictions back to EA grid
  pred_ras_ea <- rasterize(pred_sp, r_ea, field = idw_pred$var1.pred, fun = "last", background = NA)
  
  # 7) Merge predictions into r_ea only at NA cells
  filled_ea <- r_ea
  filled_ea[na_idx_ea] <- pred_ras_ea[na_idx_ea]
  
  # Clamp to [0,1] for proportions
  filled_ea[filled_ea[] < 0] <- 0
  filled_ea[filled_ea[] > 1] <- 1
  
  # 8) Project the filled raster back to the original grid *exactly*
  filled_ll <- projectRaster(from = filled_ea, to = r, method = "bilinear")
  
  # 9) Replace only NA cells in the original raster; keep known values intact
  out <- r
  na_idx_ll <- which(is.na(out[]))
  out[na_idx_ll] <- filled_ll[na_idx_ll]
  out[is.na(depth)] <- NA
  #plot(out, colNA='black')
  return(out)
}


#rasterize mapping footprint----------------------------------------------------
cat("Rasterize...")
microgrid.ras = rasterize(microgrid.pts, depth, field='Shape_Area', fun='sum', background=NA, na.rm=T)
microgrid.ras[is.na(depth)] <- NA
plot(microgrid.ras, colNA='black')
# Build a mask of mapped cells for later use
mapped_mask <- !is.na(microgrid.ras)  # TRUE where mapped, FALSE where unmapped


#rasterize habitat and calculate proportion area
newhabs = sort(unique(habitat.pts$NewHab))
newhabs = newhabs[c(2,3,1,5,6,4)]
habpct.stack <- stack()
for(i in 1:length(newhabs)){
  #i=6
  hab_class <- newhabs[i]
  hab.ras = rasterize(habitat.pts[habitat.pts$NewHab==hab_class,],depth,field='Shape_Area',fun='sum', background=0, na.rm=T)
  hab.ras[is.na(depth)] <- NA
  hab.ras[!mapped_mask] <- NA
  #plot(hab.ras,colNA='black',main=hab_class)
  habpct.ras <- hab.ras/microgrid.ras
  #plot(habpct.ras,colNA='black',main=hab_class)

  if (hab_class %in% c('NL', 'NM', 'NH')) {
    habpct.ras[is.na(depth)] <- 0
    habpct.ras <- idw_fill_raster_longlat(r=habpct.ras, idp=4, nmax=8, depth=depth, mask = NULL,
                                          ea_crs = "+proj=aea +lat_1=24 +lat_2=31.5 +lat_0=23 +lon_0=-84 +datum=WGS84 +units=m +no_defs") 
  }

  if (hab_class %in% c("AL", "AM", "AH")) {
    # Set water cells to 0 even where unmapped - i.e. do not extrapolate artificial habitat as it likely isn't a gradient
    habpct.ras[is.na(habpct.ras) & !is.na(depth)] <- 0
  }
  habpct.stack <- addLayer(habpct.stack, habpct.ras)
  
  #cleanup
  rm(hab.ras, habpct.ras, mcrgrd); gc()
}
names(habpct.stack) <- newhabs

#output-------------------------------------------------------------------------
dir.out = file.path(dir.maps,paste0(res.min,"min"))
if(!dir.exists(dir.out)) dir.create(dir.out)
writeRaster(habpct.stack,filename=paste0(dir.out,"/GFISHER"),bylayer=T,format='ascii', overwrite=T,
            suffix=paste0(names(habpct.stack),"_prop_",res.min,"min_",dim(habpct.stack)[1],"x",dim(habpct.stack)[2]))
writeRaster(microgrid.ras,filename=paste0(dir.out,"/GFISHER_microgrid_",res.min,"min_",dim(microgrid.ras)[1],"x",dim(microgrid.ras)[2],".asc"), overwrite=T)
cat(paste0("Done - Ecospace ascii files saved to ",dir.out,"\n"))
}

fn.plot_GFISHER_habitats <- function(dir.maps){
  #dir.maps=file.path(dir.maps,"10min")
  files.asc = list.files(dir.maps, pattern='.asc', full.names=T)
  files.hab = files.asc[-grep("microgrid",basename(files.asc))]
  files.grid = files.asc[grep("microgrid",basename(files.asc))]
  
  habpct.stack = stack(files.hab)
  hababbr = substring(names(habpct.stack),9,10)
  names(habpct.stack) = hababbr
  habpct.stack = habpct.stack[[c('AL','AM','AH','NL','NM','NH')]]
  hablabels <- c('Artificial, low relief','Artificial, medium relief','Artificial, high relief','Natural, low relief','Natural, medium relief','Natural, high relief')
  res.min = round(res(habpct.stack)[1]*60,0)
  
  png(filename=paste0(dir.maps,"/Habitat prop area ",res.min,"min.png"),height=7,width=7,units='in',res=600)#,compression='lzw')
  par(mfcol=c(3,2),mar=c(3,3,3,5))
  for(i in 1:nlayers(habpct.stack)){
    #i=1
    hab.i = habpct.stack[[i]]
    vals.i = sort(unique(getValues(hab.i)))
    #brks = c(0,sort(unique(getValues(hab.i)))[2],pretty(sort(unique(getValues(hab.i)))[-c(1:2)],n=50))
    #brks = c(0,pretty(sort(unique(getValues(hab.i)))[-1],n=50))
    brks = c(0,min(vals.i[vals.i>0],na.rm=T),pretty(vals.i,n=50)[-1])
    if(length(which(brks==Inf))>0) brks = brks[-which(brks==Inf)]
    colv = c(colorRamps::matlab.like2(n=length(brks)))
    col.bias=2
    #if(i==1) col.bias=3
    funpal  = colorRampPalette(colv,bias=col.bias,interpolate='spline')
    cols   = funpal(length(brks)-1)
    plot(hab.i,colNA='lightgray', main=paste0(hablabels[i]), col=cols, breaks=brks,legend=T)
    #plot(hab.i,colNA='lightgray', main=paste(newhabs[i],'pct'))
    map(database='state',add=T, fill=T, col='lightgray')
  }
  dev.off()
  
  mcrgrd = raster(files.grid)
  mcrgrd = mcrgrd/1e6
  brks = pretty(sort(unique(getValues(mcrgrd))),n=50)
  if(length(which(brks==Inf))>0) brks = brks[-which(brks==Inf)]
  colv = c(colorRamps::matlab.like2(n=length(brks)))
  col.bias=2
  funpal  = colorRampPalette(colv,bias=col.bias,interpolate='spline')
  cols   = funpal(length(brks)-1)
  
  png(filename=paste0(dir.maps,"/mapping footprint ",res.min,"min.png"),height=7,width=7,units='in',res=600)#,compression='lzw')
  par(mfcol=c(1,1),mar=c(3,3,3,6))
  plot(mcrgrd,colNA='lightgray', main="mapping footprint (sq km)\n per grid cell", col=cols, breaks=brks,legend=T)
  map(database='state',add=T, fill=T, col='lightgray')
  dev.off()
  
  cat(paste0("Map figures saved in: ",dir.maps))
}



