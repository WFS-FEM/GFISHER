# maxn_maps.R ---------------------------------------------------------------------
# Rasterize the video-survey MaxN records into one heatmap per model group, on the model
# grid. Input is the table from R/video_dataset.R; output feeds the affinity modules.
#
# Extracted verbatim from the former 'R/GFISHER functions.R'; the estimator is unchanged.

suppressPackageStartupMessages({
  library('sp'); library('raster'); library('colorRamps')
})

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
