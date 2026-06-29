library('sf')
library('sp')
#library('rgdal')
library('colorRamps')
library('maps')
library('lwgeom')
library('gstat')
library('raster')
library('xlsx')
library('reshape2')
library('truncnorm')

fn.make_gfisher_videodataset <- function(file.maxn, file.env, file.len, bbox,spplist){
  
  #fxn arguments/inputs
  # file.gfsh = "C:\\Users\\dchagaris\\Github\\WFS-FEM\\GFISHER\\data\\Video Count Data4ChagarisTake2.xlsx"
  # bbox <- bbox
  # spplist <- spplist
  
  #--------------------------------import and prepare data------------------------------------------
  # dat1 <- read_excel(path=file.gfsh, sheet='Count_Data')
  # dat.lf <- read_excel(path=file.gfsh, sheet='Length_Data')
  dat.maxn <- read.csv(file.maxn, header=T)
  names(dat.maxn) <- tolower(names(dat.maxn))
  dat.lf <- read.csv(file.len, header=T)
  names(dat.lf) <- tolower(names(dat.lf))
  dat.env <- read.csv(file.env, header=T)
  names(dat.env) <- tolower(names(dat.env))
  
  
  #dat1$taxon <- tolower(gsub("_"," ",dat1$taxon))
  #dat.lf$taxon <- tolower(gsub("_"," ",dat.lf$sciname))
  
  #---------------------------------------make dataset----------------------------------------------
  #filter for stations to use in analysis, trawls conducted with geographic box
  #use which() so rows with NA lat/lon are dropped rather than returned as phantom all-NA rows
  keep.env = which(dat.env$lat_dd>=bbox[2] & dat.env$lat_dd<=bbox[1] &
                   dat.env$lon_dd>=bbox[3] & dat.env$lon_dd<=bbox[4])
  dat.env2 = unique(dat.env[keep.env,])
  names(dat.env2)
  length(unique(dat.env$reference))
  length(unique(dat.env2$reference))
  nrow(dat.env2)
  sum(duplicated(dat.env2$reference))

  #build one station record per reference (coordinates + env), preferring complete rows.
  #NOTE: keep the 'dup'/'complete' helper columns - they are dropped by name in the merge below.
  stations = unique(dat.env2)
  stations$dup = stations$reference %in% stations$reference[duplicated(stations$reference)]
  ncheck = min(8, ncol(stations))
  stations$complete = complete.cases(stations[, seq_len(ncheck)])
  #drop duplicate-reference rows that are incomplete, then any remaining duplicate references.
  #guard each removal with length()>0 - `df[-integer(0),]` would otherwise wipe all rows.
  drop1 = which(stations$dup & !stations$complete)
  if(length(drop1)>0) stations = stations[-drop1,]
  dup2 = which(duplicated(stations$reference))
  if(length(dup2)>0) stations = stations[-dup2,]
  
  #--------------------------create species to model groupings key-------------------------------------------
  #read species-group assignments and list of model groups
  modspp  = read.xlsx(file.spplist,sheetName="spplist",stringsAsFactors=F,colIndex=1:10)
  modgrps = read.xlsx(file.spplist,sheetName="model groups",stringsAsFactors=F,colIndex=1:4)
  names(modspp)[1:2] = c('modnumber','modname')
  names(modgrps)[1:2] = c('modnumber','modname')
  spplist <- melt(modspp, id.vars=c('modnumber','modname'), measure.vars=c('species','query','og_name','class','order','family','genus'), variable.name='var',value.name='taxon')
  spplist <- spplist[,-3]
  spplist <- data.frame(lapply(spplist, tolower), stringsAsFactors = FALSE)
  spplist$taxon <- ifelse(spplist$taxon=="",NA,spplist$taxon)
  spplist <- unique(spplist)
  spplist <- spplist[spplist$modnumber!='99',]
  #spplist <- spplist[complete.cases(spplist),]
  
  ##multistanza size at age-------------------------------------------------------------------------
  sizeatage <- read.xlsx(file.spplist,sheetName="size_at_age",stringsAsFactors=F)
  
  #keep species included in the model
  keeptaxa = tolower(sort(unique(spplist$taxon)))
  names(dat.maxn)
  spp.gfsh = data.frame(taxon=sort(unique(names(dat.maxn)[-c(1:2)])),
                        taxon2 = sub("_sp$","",sort(unique(names(dat.maxn)[-c(1:2)]))),
                        taxon3 = sapply(strsplit(sort(unique(names(dat.maxn)[-c(1:2)])), "_"), `[`, 1))
  spp.gfsh$match = ifelse(spp.gfsh$taxon %in% keeptaxa | spp.gfsh$taxon2 %in% keeptaxa | spp.gfsh$taxon3 %in% keeptaxa,TRUE,FALSE)
  spp.gfsh$modnumber = ifelse(spp.gfsh$taxon %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon,spplist$taxon)],
                              ifelse(spp.gfsh$taxon2 %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon2,spplist$taxon)],
                                     ifelse(spp.gfsh$taxon3 %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon3,spplist$taxon)],NA)))
  
  spp.gfsh$modnumber = ifelse(spp.gfsh$taxon=='lutjanidae_sp',19,
                              ifelse(spp.gfsh$taxon %in% c('epinephelus_sp','epinephelus_striatus'),36,spp.gfsh$modnumber))
  spp.gfsh$modname = modgrps$modname[match(spp.gfsh$modnumber, modgrps$modnumber)]
  
  write.csv(spp.gfsh,file.path(dirname(dirname(file.maxn)),'GFISHER_species_fg.csv'),row.names=F)
  
  #melt video data----------------------------------------------------------------------------------
  dat.maxn.long <- melt(dat.maxn[,-1],id.vars='reference', variable.name='sciname', value.name='maxn')
  dat.maxn.long <- dat.maxn.long[dat.maxn.long$maxn>0,]
  dat.maxn.long$sciname <- gsub("_"," ",as.character(dat.maxn.long$sciname))
  
  #get size data for multistanza species------------------------------------------------------------
  multistanza.fg <- as.numeric(unlist(strsplit(sizeatage$fg[sizeatage$stanzas!="0"],"-")))
  multistanza.spp <- tolower(sizeatage$sciname[sizeatage$stanzas!='0'])[c(1,3,4)]
  dat.lf$sciname <- tolower(gsub("_"," ",dat.lf$sciname))
  
  lf2 = dat.lf[dat.lf$sciname %in% multistanza.spp,c('reference','sciname','length_mm')]
  lf.spp.mean <- aggregate(length_mm~sciname, lf2, mean)
  lf.spp.sd <- aggregate(length_mm~sciname, lf2, sd)
  lf.spp.min <- aggregate(length_mm~sciname, lf2, min)
  lf.spp.max <- aggregate(length_mm~sciname, lf2, max)
  lf.spp.cnt <- aggregate(length_mm~sciname, lf2, length)
  lf.spp.sd$length_mm[is.na(lf.spp.sd$length_mm)] <- lf.spp.mean$length_mm[is.na(lf.spp.sd$length_mm)]* mean(lf.spp.sd[,2]/lf.spp.mean[,2],na.rm=T)
  
  names(lf2)
  names(dat.maxn.long)
  
  lf3 <- merge(lf2, dat.maxn.long)
  
  dat3 <- dat.maxn.long[,c('reference','sciname','maxn')]
  dat3$sciname <- gsub("_"," ",dat3$sciname)
  sort(unique(dat3$sciname))
  dat3 <- dat3[dat3$sciname %in% multistanza.spp,]
  
  dat4 <- data.frame()
  for(i in 1:nrow(dat3)){
    #i=1
    #i=which(dat3$reference=='2024_NCO-131')
    dat.i <- dat3[rep(i,dat3$maxn[i]),]
    ref.i <- dat3$reference[i]
    spp.i <- dat3$sciname[i]
    mean.i <- lf.spp.mean[lf.spp.mean$sciname==spp.i,2]
    sd.i <- lf.spp.sd[lf.spp.mean$sciname==spp.i,2]
    min.i <- lf.spp.min[lf.spp.mean$sciname==spp.i,2]
    max.i <- lf.spp.max[lf.spp.mean$sciname==spp.i,2]
    
    obslen.i <- round(lf3$length_mm[lf3$sciname==spp.i & lf3$reference==ref.i])
    if(length(obslen.i)==0){
      #len.i = round(rnorm(nrow(dat.i),mean=lf.spp.mean$length_mm[lf.spp.mean$taxon==spp.i], sd=lf.spp.sd$length_mm[lf.spp.sd$taxon==spp.i]))
      len.i = round(rtruncnorm(nrow(dat.i),mean=mean.i, sd=sd.i, a=min.i, b=max.i))
      lentype = rep('rand',length(len.i))
    } else if(length(obslen.i)>=nrow(dat.i)){
      len.i = obslen.i[sample.int(nrow(dat.i))]
      lentype = rep('obs',length(len.i))
    } else{
      #obs.i = obslen.i[sample.int(length(obslen.i))]  #sample(obslen.i,length(obslen.i),replace=F)
      if(length(obslen.i)>=3) rand.i = round(rtruncnorm(nrow(dat.i)-length(obslen.i),mean=mean(obslen.i), sd=sd(obslen.i), a=min.i, b=max.i))
      if(length(obslen.i)<3) rand.i =  round(rtruncnorm(nrow(dat.i)-length(obslen.i),mean=mean.i, sd=sd.i, a=min.i, b=max.i))
      #rand.i =  round(rtruncnorm(nrow(dat.i)-length(obslen.i),mean=mean.i, sd=sd.i, a=min.i, b=max.i))
      len.i = c(obslen.i,rand.i)
      lentype = c(rep('obs',length(obslen.i)),rep('rand',length(rand.i)))
    }
    dat.i$len_mm <- round(len.i)
    dat.i$lentype <- lentype
    dat4 <- rbind(dat4,dat.i)
  }
  
  # par(mfrow=c(2,2))
  # for(i in 1:length(multistanza.spp)){
  #   dat.i = dat4[dat4$taxon==multistanza.spp[i],]
  #   hist(dat.i$len_mm, main=multistanza.spp[i], breaks=40)
  # }
  
  #--------------------------assign multistanza species to fg -------------------------------------
  for(i in 1:nrow(dat4)){
    #i=1
    spp.i = tolower(dat4$sciname[i])
    laa.i = as.numeric(sizeatage[tolower(sizeatage$sciname)==spp.i,which(substr(names(sizeatage),1,3)=='age')])
    stanzas.i = as.numeric(unlist(strsplit(sizeatage$stanzas[tolower(sizeatage$sciname)==spp.i],"-")))
    groups.i = as.numeric(unlist(strsplit(sizeatage$fg[tolower(sizeatage$sciname)==spp.i],"-")))
    names(dat4)
    size.i = dat4$len_mm[i]/10
    age.i = which.min(abs(laa.i-size.i))-1
    stz.i = tail(which(age.i-stanzas.i>=0),1)
    grp.i = groups.i[stz.i]
    grpname.i = modgrps$modname[grp.i]
    
    dat4$modnumber[i] <- grp.i
    dat4$modname[i] <- grpname.i
    
  }
  names(dat4)
  dat4$n_at_length = 1
  dat4 <- aggregate(n_at_length~reference+sciname+modnumber+modname, data=dat4, sum)
  dat4$maxn <- dat4$n_at_length
  dat4$n_at_length <- NULL
  
  #--------------------------assign NON-multistanza species to fg -------------------------------------
  #resume here...
  dat5 <-  dat.maxn.long[!dat.maxn.long$sciname %in% multistanza.spp,c('reference','sciname','maxn')]
  names(dat5)
  names(spp.gfsh)
  spp.key <- spp.gfsh[,c('taxon','modnumber','modname')]
  spp.key$taxon <- gsub("_"," ", spp.key$taxon)          # match dat5$sciname format (spaces, not underscores)
  dat5 <- merge(dat5, spp.key, by.x='sciname',by.y='taxon', all.x=T)
  
  
  #put it back together
  names(dat4)
  names(dat5)
  dat.full <- rbind(dat4, dat5)
  
  #check counts
  sum1 <- aggregate(maxn~sciname, data=dat.full, sum)
  names(sum1)[2] <- 'final_maxn'
  sum2 <- aggregate(maxn~sciname, data=dat.maxn.long, sum)
  names(sum2)[2] <- 'raw_maxn'
  chk <- merge(sum1,sum2)
  err <- which(chk$final_maxn != chk$raw_maxn)
  if(length(err)>0){
    message('MaxN counts did not sum back up after processing')
  }
  
  #merge back with env data
  which(duplicated(stations$reference))
  dat.full <- merge(dat.full, stations[,-which(names(stations) %in% c('dup','complete'))], all.x=T, by='reference')
  return(dat.full)
} #eof

fn.make_GFISHER_maxn_maps <- function(maxn, depth, lon.col='lon_dd', lat.col='lat_dd', fun='sum', background=NA, dir.out=NULL, save.format='ascii', plot=FALSE){
  # Rasterize MaxN counts into per-model-group heatmaps on the depth grid.
  #
  # maxn    : data.frame returned by fn.make_gfisher_videodataset(); must contain
  #           'modnumber', 'maxn', and station coordinates (lon.col/lat.col, decimal degrees WGS84).
  # depth   : template raster defining output dimensions, extent, and CRS.
  # lon.col,
  # lat.col : names of the longitude/latitude columns in `maxn`.
  # fun     : aggregation applied to records sharing a cell. 'sum' (default) gives a total
  #           MaxN-count heatmap; 'mean' gives mean count per observation.
  # background : value for water cells with no observation. NA (default) leaves them empty;
  #           use 0 to treat unsampled water as zero count. Land/no-depth cells are always NA.
  # dir.out : if not NULL, write the per-group .asc rasters to this directory. Created if needed.
  # save.format : 'ascii' -> one .asc per group (Ecospace grid format), named by modnumber+modname
  #                          (default). The on-screen/PDF figure is controlled separately by `plot`.
  # plot    : if TRUE, render the per-group heatmaps to a multipage PDF (3x3 panels per page,
  #           one panel per model group) written to dir.out (or the working dir if dir.out is NULL).
  #
  # Returns a RasterStack with one layer per modnumber (named 'mod<modnumber>'), masked to the
  # depth grid. Empty (unsampled) water cells are NA; land/no-depth cells are NA. A modnumber->
  # modname lookup is attached as attr(<stack>, 'modlabels') for labelling.

  #checks-------------------------------------------------------------------------
  req <- c('modnumber','maxn',lon.col,lat.col)
  miss <- req[!req %in% names(maxn)]
  if(length(miss)>0) stop(paste('maxn is missing required column(s):', paste(miss, collapse=', ')))

  # depth may arrive as a terra SpatRaster (if terra is loaded); coerce to a raster::RasterLayer
  # so the sp-points rasterize() path below dispatches to the raster method, not terra's.
  if(inherits(depth,'SpatRaster')) depth <- raster::raster(depth)

  #drop records with no group assignment, no maxn, or no location
  keep <- !is.na(maxn$modnumber) & !is.na(maxn$maxn) &
          !is.na(maxn[[lon.col]]) & !is.na(maxn[[lat.col]])
  if(sum(!keep)>0) message(paste('Dropping',sum(!keep),'record(s) with missing modnumber, maxn, or coordinates.'))
  dat <- maxn[keep,]
  if(nrow(dat)==0) stop('No records with non-missing modnumber, maxn, and coordinates.')

  #build spatial points and match the depth CRS-----------------------------------
  dat.sp <- dat
  coordinates(dat.sp) <- as.formula(paste0('~',lon.col,'+',lat.col))
  proj4string(dat.sp) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  dat.sp <- spTransform(dat.sp, crs(depth))

  #rasterize one layer per model group--------------------------------------------
  mods <- sort(unique(dat.sp$modnumber))
  maxn.stack <- stack()
  for(i in seq_along(mods)){
    mod.i <- mods[i]
    pts.i <- dat.sp[dat.sp$modnumber==mod.i,]
    ras.i <- raster::rasterize(pts.i, depth, field='maxn', fun=fun, background=background, na.rm=T)
    ras.i[is.na(depth)] <- NA            # mask land / no-depth cells
    maxn.stack <- addLayer(maxn.stack, ras.i)
  }
  #attach modnumber -> modname lookup for labelling
  modlabels <- NULL
  if('modname' %in% names(dat)){
    modlabels <- unique(dat[,c('modnumber','modname')])
    modlabels <- modlabels[match(mods, modlabels$modnumber),]
    attr(maxn.stack,'modlabels') <- modlabels
  }

  #name layers by modname (raster sanitizes spaces to '.'); fall back to modnumber
  labs <- if(!is.null(modlabels)) modlabels$modname else paste0('mod', mods)
  names(maxn.stack) <- labs

  #naming bits shared by the save and plot blocks
  res.min <- round(res(depth)[1]*60,0)
  dims    <- paste0(dim(maxn.stack)[1],'x',dim(maxn.stack)[2])

  #save---------------------------------------------------------------------------
  if(!is.null(dir.out)){
    if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
    # one .asc per group; encode modnumber + sanitized modname in the filename since ascii drops names
    # base 'GFISHER_maxn' + suffix -> GFISHER_maxn_mod<n>_<modname>_<res>min_<dims>.asc
    suff <- paste0('mod', mods, '_', gsub('[^A-Za-z0-9]+','-', labs), '_', res.min, 'min_', dims)
    raster::writeRaster(maxn.stack, filename=file.path(dir.out,'GFISHER_maxn'),
                        bylayer=TRUE, suffix=suff, format='ascii', overwrite=TRUE)
    message(paste0('Saved maxn maps (ascii) to ', dir.out))
  }

  #plot to a multipage PDF (3x3 panels per page)----------------------------------
  # A direct plot() side-effect on a multi-layer stack is unreliable in scripts (recording
  # devices, terra masking raster's plot generic), so render straight to a PDF device instead.
  if(plot){
    pdf.dir <- if(!is.null(dir.out)) dir.out else getwd()
    if(!dir.exists(pdf.dir)) dir.create(pdf.dir, recursive=TRUE)
    pdf.file <- file.path(pdf.dir, paste0('GFISHER_maxn_heatmaps_',res.min,'min_',dims,'.pdf'))
    pdf(pdf.file, onefile=TRUE, width=10, height=10)
    op <- par(mfrow=c(3,3), mar=c(3,3,3,5))
    for(i in 1:nlayers(maxn.stack)){
      plot(maxn.stack[[i]], colNA='black', main=labs[i])
    }
    par(op)
    dev.off()
    message(paste0('Saved maxn heatmap pdf to ', pdf.file))
  }
  return(maxn.stack)
} #eof

fn.make_GFISHER_habitat_maps <- function(file.gdb, dir.maps, depth=depth){
  
  #import and prepare geodatabase file--------------------------------------------
  # The directory contains a zip file with a geodatabase and supporting metadata to describe the datasets. Two shapefiles are in the geodatabase:
  # 1.	East_Master_microgrid_mapped – this shapefile consists of microgrids (0.1 x 0.1 nm) which we’ve mapped with side scan sonar.  Each grid represents one of our primary sampling units that has been scanned (mapping footprint).
  # 2.	FWRI_East_Master_Hab_Data – this shapefile consists of our digitized polygons for our habitat classes (i.e. Geoforms).  The habitat classes are described in the metadata file and consist of both natural and artificial reefs.
  res.min = round(res(depth)[1]*60,0)
  lyrs.gdb = st_layers(file.gdb)
  #lyrs.gdb.micro = st_layers(file.gdb.micro)
  
  #extract shapefile from gdb
  cat("Reading Geodatabase Files...\n")
  microgrid <- st_read(file.gdb, layer=lyrs.gdb$name[grep("Microgrid",lyrs.gdb$name)])   
  habitat <- st_read(file.gdb, layer=lyrs.gdb$name[grep("Dissolve",lyrs.gdb$name)])
  habitat <- habitat[habitat$NewHabStrat!='AP',]
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
  message("Check\n")
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
  message(paste("N habitat records with matching microgrids:",length(which(hab$MicroGrid %in% mcg$MicroGrid)),"\n"))
  message(paste("N habitat records with non-matching microgrid:",nrow(hab)-length(which(hab$MicroGrid %in% mcg$MicroGrid))-length(which(hab$MicroGrid==" ")),"\n"))
  message(paste("N habitat records with missing microgrid:",length(which(hab$MicroGrid==" ")),"\n"))
  message(paste("N microgrids without any habitat records:",length(which(!mcg$MicroGrid %in% hab$MicroGrid)),"\n"))
  
  #does habitat area exceed mapped area for any microgrids?
  hab.sum <- aggregate(Hab_Area~MicroGrid, data=hab, sum)
  chk1 <- merge(hab.sum,mcg[,c('MicroGrid','Grid_Area')],by='MicroGrid',all=T)
  chk1$habpct <- round(chk1$Hab_Area/chk1$Grid_Area,4)
  message(paste("N microgrids with Hab_Area>Grid_Area:",length(which(chk1$habpct>1)),"\n"))
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
  message("Rasterize...")
  microgrid.ras = rasterize(microgrid.pts, depth, field='Shape_Area', fun='sum', background=NA, na.rm=T)
  microgrid.ras[is.na(depth)] <- NA
  plot(microgrid.ras, colNA='black')
  # Build a mask of mapped cells for later use
  mapped_mask <- !is.na(microgrid.ras)  # TRUE where mapped, FALSE where unmapped
  
  
  #rasterize habitat and calculate proportion area
  newhabs = sort(unique(habitat.pts$NewHab))
  newhabs = newhabs[c(2,3,1,5,6,4)]
  newhabs
  habpct.stack <- stack()
  for(i in 1:length(newhabs)){
    #i=6
    hab_class <- newhabs[i]
    hab.ras = rasterize(habitat.pts[habitat.pts$NewHab==hab_class,],depth,field='Shape_Area',fun='sum', background=0, na.rm=T)
    hab.ras[is.na(depth)] <- NA
    hab.ras[!mapped_mask] <- NA
    #plot(hab.ras,colNA='black',main=hab_class)
    habpct.ras <- hab.ras/microgrid.ras
    #habpct.ras[is.na(habpct.ras)] <- 0
    habpct.ras[is.na(depth) | depth>200] <- 0

    if (hab_class %in% c('NL', 'NM', 'NH')) {
      #habpct.ras[is.na(depth)] <- 0
      habpct.ras <- idw_fill_raster_longlat(r=habpct.ras, idp=4, nmax=8, depth=depth, mask = NULL,
                                            ea_crs = "+proj=aea +lat_1=24 +lat_2=31.5 +lat_0=23 +lon_0=-84 +datum=WGS84 +units=m +no_defs") 
      #plot(habpct.ras,colNA='black',main=hab_class)
    }
    
    if (hab_class %in% c("AL", "AM", "AH")) {
      # Set water cells to 0 even where unmapped - i.e. do not extrapolate artificial habitat as it likely isn't a gradient
      habpct.ras[is.na(habpct.ras) & !is.na(depth)] <- 0
    }
    habpct.stack <- addLayer(habpct.stack, habpct.ras)
    
    #cleanup
    rm(hab.ras, habpct.ras); gc()
  }
  names(habpct.stack) <- newhabs
  library('terra')
  plot(habpct.stack)
  
  #output-------------------------------------------------------------------------
  dir.out = file.path(dir.maps,paste0(res.min,"min"))
  if(!dir.exists(dir.out)) dir.create(dir.out)
  writeRaster(habpct.stack,filename=paste0(dir.out,"/GFISHER"),bylayer=T,format='ascii', overwrite=T,
              suffix=paste0(names(habpct.stack),"_prop_",res.min,"min_",dim(habpct.stack)[1],"x",dim(habpct.stack)[2]))
  writeRaster(microgrid.ras,filename=paste0(dir.out,"/GFISHER_microgrid_",res.min,"min_",dim(microgrid.ras)[1],"x",dim(microgrid.ras)[2],".asc"), overwrite=T)
  message(paste0("Done - Ecospace ascii files saved to ",dir.out,"\n"))
} #eof


fn.make_GFISHER_habitat_maps_old <- function(file.gdb, file.gdb.micro, dir.maps, depth=depth.10min){

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
} #eof

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
} #eof





