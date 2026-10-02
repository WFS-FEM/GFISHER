# make_habitat_basemaps.R ---------------------------------------------------------
# Build a complete, sum-to-1 Ecospace habitat basemap from the GFISHER side-scan mapping
# plus dbSeabed, replacing the six standalone reef-proportion rasters written by
# fn.make_GFISHER_habitat_maps().
#
# LAYERS (each a proportion of cell area; together they sum to 1 over water):
#   AL AM AH NL NM NH   GFISHER reef classes, artificial/natural x low/medium/high relief
#   RCK                 rock outcrop, from raw dbSeabed
#   UNC                 unconsolidated bottom (optionally split into GVL / SND / MUD)
#   SGR                 seagrass, from the existing seagrass raster
#
# ---------------------------------------------------------------------------------
# FOUR THINGS THIS DOES DIFFERENTLY FROM THE OLD PIPELINE, AND WHY
# ---------------------------------------------------------------------------------
# 1. REEF IS SCALED TO THE SCANNED AREA, NOT THE CELL.
#    GFISHER has scanned 2.95% of the domain: 39% of water cells contain any microgrid at
#    all, and the median such cell is only 5.5% scanned. Dividing reef area by CELL area --
#    what the old pipeline does -- reports reef found in 5% of a cell as though the whole
#    cell had been searched, understating reef by roughly 11x on average. The scanned
#    microgrids are the sample and the cell is the population, so the density estimator is
#    reef_area / scanned_area.
#
# 2. UNMAPPED GROUND IS SHRUNK TOWARD A STRATUM MEAN.
#    Each cell's density is an empirical-Bayes blend of its own scanned density and a
#    borrowing target, weighted by how much of it was scanned:
#        reef = (M*d_cell + k*target) / (M + k)
#    Heavily-scanned cells keep their own value; unscanned cells fall back to the target.
#    The default target is the pooled density of the cell's (region x depth-bin) stratum.
#
#    Two alternatives were built and rejected, and are kept only for comparison:
#      target='smooth'  a GAM, d ~ s(depth) + s(lon,lat). Badly miscalibrated: its s(lon,lat)
#                       term is unanchored in unscanned water, and on the logit scale it ran
#                       away to predict NH at 0.176 of area domain-wide against an observed
#                       0.0037 -- 48x too high. Only the depth mask hid it. DO NOT USE.
#      target='idw'     local inverse-distance mean with a hard search radius. Well calibrated,
#                       but produces a bleeding/halo effect around mapped ground that reads as
#                       an artifact of where the survey went rather than as habitat.
#
#    The stratum mean is blocky by construction -- roughly 61% of water cells end up holding
#    one of ~10 distinct values, so plateaus with hard edges on the stratum boundaries are
#    visible in the maps. That is a known and accepted cost: the assumption is simple, states
#    plainly that unscanned ground is being assigned its stratum's average, and does not
#    pretend to spatial detail the survey cannot support.
#
#    CAVEAT, and it is not resolvable from the data: FWRI scanned where they expected
#    structure, so stratum means are computed over reef-oriented ground. Extending them to
#    unmapped area probably overestimates. `damp.borrowed` scales the borrowed component so
#    that assumption is a visible dial rather than a hidden one. It defaults to 1 (undamped).
#
# 3. ROCK COMES FROM THE RAW dbSEABED GRIDS, NOT THE PROCESSED 5-MIN LAYERS.
#    dbSeabed is NOT a four-part composition. GVL+SND+MUD is a closed grain-size triangle
#    (sums to 98-100 on its own, in every cell that has any rock), and RCK is a separate
#    areal fraction of rock outcrop. The processed gmf_*_prop_* layers renormalize all four
#    to sum to 1, which costs 42.6% of the rock signal; and 04_make_ascii_dbSEABED.R maps
#    the -99 NODATA flag to 0, fabricating "no rock" across 69% of native cells. Reading the
#    raw grids at their native 1.2 arc-min and preserving NODATA raises mean rock from 0.039
#    to 0.106 and takes cells at >=50% rock from 29 to 200.
#
# 4. GFISHER REEF AND dbSEABED ROCK ARE NORMALIZED INDEPENDENTLY, NOT SUBTRACTED.
#    They are uncorrelated (Spearman 0.004 at the dbSeabed cell scale, aggregating ~36
#    microgrids per cell, so this is not a scale artifact). dbSeabed rock carries no
#    information about where GFISHER reef is -- it is regional hardbottom interpolated from
#    sparse samples, while GFISHER digitizes discrete geoforms with relief. Subtracting one
#    from the other would clip to zero wherever GFISHER mapped reef that dbSeabed missed, so
#    the layers are built as raw fractions and divided through by their sum. Some
#    double-counting of hardbottom is the accepted cost.

suppressPackageStartupMessages({library('sf'); library('raster'); library('sp')})

REEF.CLASSES <- c('AL','AM','AH','NL','NM','NH')
TRIANGLE     <- c('GVL','SND','MUD')

#--- per-microgrid reef area by class -------------------------------------------------------
# Returns one row per mapped microgrid: id, centroid, scanned area, region/depth stratum, and
# the reef area of each of the six classes.
fn.microgrid_reef <- function(file.gdb){
  ly <- sf::st_layers(file.gdb)
  l.mg <- ly$name[grep('Microgrid', ly$name, ignore.case=TRUE)]
  l.hb <- ly$name[grep('Hab_Data_FINAL', ly$name, ignore.case=TRUE)]
  if(length(l.mg) != 1) stop('could not find a unique Microgrid layer in ', file.gdb)
  if(length(l.hb) != 1) stop('could not find a unique Hab_Data_FINAL layer in ', file.gdb)
  mg <- sf::st_drop_geometry(sf::st_read(file.gdb, layer=l.mg, quiet=TRUE))
  hb <- sf::st_drop_geometry(sf::st_read(file.gdb, layer=l.hb, quiet=TRUE))
  message('  ', nrow(mg), ' mapped microgrids, ', nrow(hb), ' habitat polygons')

  cls.col <- grep('^NewHab$', names(hb), value=TRUE)
  if(!length(cls.col)){
    strat.col <- grep('^NewHabStra', names(hb), value=TRUE)[1]
    hb$NewHab <- substr(hb[[strat.col]], 1, 2)
  }
  hb$NewHab <- substr(hb$NewHab, 1, 2)
  hb <- hb[hb$NewHab %in% REEF.CLASSES & !is.na(hb$MicroGrid), ]
  hb$area <- hb[[grep('^Shape_Area$', names(hb), value=TRUE)[1]]]

  a <- tapply(hb$area, list(hb$MicroGrid, hb$NewHab), sum)
  A <- matrix(0, nrow(mg), length(REEF.CLASSES),
              dimnames=list(NULL, REEF.CLASSES))
  i <- match(rownames(a), mg$MicroGrid)
  for(k in intersect(colnames(a), REEF.CLASSES)){
    v <- a[, k]; ok <- !is.na(i) & !is.na(v)
    A[i[ok], k] <- v[ok]
  }
  dep <- if('Depth_bin' %in% names(mg)) trimws(mg$Depth_bin) else NA_character_
  reg <- if('Region'    %in% names(mg)) trimws(mg$Region)    else substr(trimws(mg$SpaceStrat),1,2)
  out <- data.frame(MicroGrid = mg$MicroGrid, X = mg$X, Y = mg$Y,
                    scanned   = mg[[grep('^Shape_Area$', names(mg), value=TRUE)[1]]],
                    region = reg, depth_bin = dep,
                    stringsAsFactors = FALSE)
  cbind(out, as.data.frame(A))
}

#--- aggregate microgrids onto the model grid ------------------------------------------------
# Sums reef area by class and total scanned area into each template cell, and records the
# modal region of the microgrids falling in it.
fn.reef_to_grid <- function(mgr, template){
  cell <- raster::cellFromXY(template, cbind(mgr$X, mgr$Y))
  ok <- !is.na(cell)
  if(!any(ok)) stop('no microgrid centroids fall inside the template grid')
  mgr <- mgr[ok, ]; cell <- cell[ok]
  scanned <- tapply(mgr$scanned, cell, sum)
  A <- sapply(REEF.CLASSES, function(k) tapply(mgr[[k]], cell, sum))
  A[is.na(A)] <- 0
  modal <- function(x){ t <- table(x[nzchar(x) & !is.na(x)]); if(!length(t)) NA_character_ else names(t)[which.max(t)] }
  reg <- tapply(mgr$region, cell, modal)
  data.frame(cell = as.integer(names(scanned)), scanned = as.numeric(scanned),
             region = as.character(reg), A, stringsAsFactors = FALSE)
}

#--- strata for every water cell -------------------------------------------------------------
# region comes from the mapped microgrids where available and from the nearest mapped cell
# otherwise; depth bins are quantiles of depth over water. Both are always defined, so every
# cell has somewhere to borrow from.
fn.cell_strata <- function(template, grid.reef, n.dbin=5){
  v <- raster::getValues(template)
  water <- which(!is.na(v))
  xy <- raster::xyFromCell(template, water)
  reg <- rep(NA_character_, length(water))
  m <- match(water, grid.reef$cell)
  reg[!is.na(m)] <- grid.reef$region[m[!is.na(m)]]
  known <- which(!is.na(reg))
  if(length(known) && any(is.na(reg))){
    nn <- FNN::get.knnx(xy[known, , drop=FALSE], xy[is.na(reg), , drop=FALSE], k=1)$nn.index[,1]
    reg[is.na(reg)] <- reg[known][nn]
  }
  if(all(is.na(reg))) reg[] <- 'ALL'
  br <- unique(stats::quantile(abs(v[water]), seq(0, 1, length.out=n.dbin+1), na.rm=TRUE))
  dbin <- cut(abs(v[water]), br, include.lowest=TRUE, labels=FALSE)
  data.frame(cell = water, lon = xy[,1], lat = xy[,2], depth = v[water], region = reg,
             stratum = paste0(reg, '_d', dbin), stringsAsFactors = FALSE)
}

#--- estimate the shrinkage constant from the variance structure ------------------------------
# The empirical-Bayes k is the ratio of sampling variance to between-cell variance. Writing
# d_i for the density observed in cell i from scanned area M_i,
#     Var(d_i - t_stratum) = tau2 + sigma2 / M_i
# so regressing the squared deviations on 1/M_i gives sigma2 (slope) and tau2 (intercept),
# and k = sigma2 / tau2. Large k means cells differ little relative to how noisily each is
# measured, so borrow heavily; small k means real between-cell structure worth preserving.
# Returns NA when the fit is degenerate (too few mapped cells, or a non-positive component),
# and the caller falls back to the median scanned area.
fn.estimate_k <- function(M, A, stratum){
  i <- which(M > 0)
  if(length(i) < 30) return(NA_real_)
  d <- A[i]/M[i]
  num <- tapply(A[i], stratum[i], sum); den <- tapply(M[i], stratum[i], sum)
  t.s <- (num/den)[stratum[i]]; t.s[is.na(t.s)] <- mean(d)
  y <- (d - t.s)^2; x <- 1/M[i]
  fit <- try(stats::lm(y ~ x), silent=TRUE)
  if(inherits(fit, 'try-error')) return(NA_real_)
  tau2 <- unname(coef(fit)[1]); sig2 <- unname(coef(fit)[2])
  if(!is.finite(tau2) || !is.finite(sig2) || tau2 <= 0 || sig2 <= 0) return(NA_real_)
  sig2/tau2
}

#--- smooth borrowing target ------------------------------------------------------------------
# The quantity a cell shrinks TOWARD. Binning region x depth into strata makes that target
# piecewise constant: on this grid 61% of water cells ended up holding one of just 10 distinct
# values, which renders as flat plateaus with hard edges on the stratum boundaries -- a mapping
# artifact, not a biological pattern.
#
# Instead, fit the observed cell densities as a smooth surface,
#     d ~ s(depth) + s(lon, lat),   quasibinomial, weights = scanned area,
# and predict everywhere. Depth enters continuously rather than binned, space varies smoothly
# rather than by region block, and the weighting by scanned area means well-mapped cells drive
# the fit. Predictions are clamped to the observed density range so the spatial smooth cannot
# extrapolate wildly into ground far from any mapping.
#
# Returns NULL when the class is too sparse to fit, and the caller falls back to stratum means.
fn.smooth_target <- function(strata, M, A, k.sp=30, k.dep=5){
  if(!requireNamespace('mgcv', quietly=TRUE)) return(NULL)
  i <- which(M > 0)
  if(length(i) < 50) return(NULL)
  d <- A[i]/M[i]
  if(sum(d > 0) < 20 || max(d) <= 0) return(NULL)
  df <- data.frame(y = pmin(pmax(d, 0), 1), w = M[i],
                   depth = abs(strata$depth[i]), lon = strata$lon[i], lat = strata$lat[i])
  k.sp <- min(k.sp, max(10, floor(nrow(df)/10)))
  fit <- try(suppressWarnings(mgcv::gam(
           y ~ s(depth, k=k.dep) + s(lon, lat, k=k.sp),
           data=df, weights=df$w, family=stats::quasibinomial())), silent=TRUE)
  if(inherits(fit, 'try-error')) return(NULL)
  nd <- data.frame(depth = abs(strata$depth), lon = strata$lon, lat = strata$lat)
  p <- try(as.numeric(mgcv::predict.gam(fit, newdata=nd, type='response')), silent=TRUE)
  if(inherits(p, 'try-error') || any(!is.finite(p))) return(NULL)
  pmin(pmax(p, 0), max(d))          # never extrapolate above the densest cell observed
}

#--- local IDW borrowing target ----------------------------------------------------------------
# Inverse-distance weighted mean of the OBSERVED densities in the nearest mapped cells, with a
# hard search radius. Two properties the GAM surface does not have:
#   - it never predicts habitat further than `maxdist` from a cell that was actually scanned;
#     beyond that radius the target is 0, because there is no evidence either way and inventing
#     reef there is worse than declining to.
#   - it cannot manufacture a nearshore gradient out of offshore observations: the GAM's
#     s(lon,lat) term is unanchored in unscanned water and happily extended NH right to the
#     shoreline, which is not where high-relief natural reef lives.
# Neighbours are weighted by 1/dist^idp AND by scanned area, so a well-mapped neighbour counts
# for more than a barely-scanned one.
# maxdist is in decimal degrees on the model grid (1 deg is roughly 110 km). It is by far the
# strongest control on how far predictions extend: 0.5/1/2 deg fill 27%/49%/80% of unmapped
# cells. idp is a weak lever by comparison -- moving it from 2 to 0.5 shifts the mean ~5%.
# Defaults are tuned so the domain mean stays close to the density actually observed in mapped
# cells; widening maxdist much beyond 1 deg pushes it above that, since IDW then fills distant
# cells from the reef-targeted mapped set.
fn.idw_target <- function(strata, M, A, idp=1.5, nmax=16, maxdist=1){
  if(!requireNamespace('FNN', quietly=TRUE)) return(NULL)
  i <- which(M > 0)
  if(length(i) < 5) return(NULL)
  d <- A[i]/M[i]
  src <- cbind(strata$lon[i], strata$lat[i])
  tgt <- cbind(strata$lon,    strata$lat)
  kk <- min(nmax, length(i))
  nn <- FNN::get.knnx(src, tgt, k=kk)
  w  <- 1/pmax(nn$nn.dist, 1e-6)^idp
  w[nn$nn.dist > maxdist] <- 0                 # hard radius: no evidence, no habitat
  w <- w * matrix(M[i][nn$nn.index], nrow(w), ncol(w))   # precision-weight by scanned area
  sw <- rowSums(w)
  out <- ifelse(sw > 0, rowSums(w * matrix(d[nn$nn.index], nrow(w), ncol(w))) / sw, 0)
  # a mapped cell's own observation is at distance 0 and dominates its own target, which is
  # fine -- the caller blends it with that same observation anyway.
  pmin(pmax(out, 0), max(d))
}

#--- empirical-Bayes reef density per class --------------------------------------------------
# reef_h(i) = (M_i * d_h(i) + k * damp * t_h(s_i)) / (M_i + k)
#   M_i = scanned area in cell i, d_h(i) = reef_h area / M_i, t_h(s) = pooled stratum density,
#   k   = the shrinkage constant, in m2 of scanned area. A cell weights its own measurement
#         M/(M+k), so k is literally "how much scanning it takes to half-trust a cell".
# Unscanned cells (M_i = 0) collapse to damp * t_h(s_i).
#
# k accepts:
#   'auto'    per-class estimate from fn.estimate_k (falls back to 'median' where degenerate)
#   'median'  the median scanned area over mapped cells -- median cell gets a 50/50 split
#   numeric   a fixed value in m2, applied to every class
#
# target accepts:
#   'idw'     local inverse-distance mean of observed densities, zero beyond maxdist (default)
#   'smooth'  GAM surface s(depth)+s(lon,lat) -- smooth, but extrapolates into unscanned water
#   'stratum' shrink toward the (region x depth-bin) pooled mean  -- piecewise constant
fn.shrink_reef <- function(strata, grid.reef, damp.borrowed=1, k='median', target='stratum',
                           idp=1.5, nmax=16, maxdist=1){
  n <- nrow(strata)
  M <- rep(0, n); Aa <- matrix(0, n, length(REEF.CLASSES), dimnames=list(NULL, REEF.CLASSES))
  m <- match(grid.reef$cell, strata$cell); ok <- !is.na(m)
  M[m[ok]] <- grid.reef$scanned[ok]
  Aa[m[ok], ] <- as.matrix(grid.reef[ok, REEF.CLASSES])
  k.med <- stats::median(M[M > 0])
  message('  ', sum(M > 0), ' of ', n, ' water cells have mapping; median scanned area = ',
          format(round(k.med), big.mark=','), ' m2')
  out <- matrix(NA_real_, n, length(REEF.CLASSES), dimnames=list(NULL, REEF.CLASSES))
  kk  <- setNames(rep(NA_real_, length(REEF.CLASSES)), REEF.CLASSES)
  for(h in REEF.CLASSES){
    k.h <- if(is.numeric(k)) k[[1]] else if(identical(k, 'auto'))
             fn.estimate_k(M, Aa[, h], strata$stratum) else k.med
    if(!is.finite(k.h) || k.h < 0){
      if(identical(k, 'auto')) message('    ', h, ': variance fit degenerate, using median')
      k.h <- k.med
    }
    kk[h] <- k.h
    t.s <- NULL
    if(identical(target, 'smooth')){
      t.s <- fn.smooth_target(strata, M, Aa[, h])
      if(is.null(t.s)) message('    ', h, ': smooth target would not fit, using stratum means')
    } else if(identical(target, 'idw')){
      t.s <- fn.idw_target(strata, M, Aa[, h], idp=idp, nmax=nmax, maxdist=maxdist)
      if(is.null(t.s)) message('    ', h, ': idw target would not fit, using stratum means')
    }
    if(is.null(t.s)){
      num <- tapply(Aa[, h], strata$stratum, sum)
      den <- tapply(M,       strata$stratum, sum)
      t.s <- ifelse(den > 0, num/den, 0)[strata$stratum]
      t.s[is.na(t.s)] <- 0
    }
    d.i <- ifelse(M > 0, Aa[, h]/M, 0)
    out[, h] <- (M * d.i + k.h * damp.borrowed * t.s) / (M + k.h)
  }
  message('    k by class (m2, and as % of median scanned area):')
  message('      ', paste(sprintf('%s=%s (%.0f%%)', REEF.CLASSES,
          format(round(kk[REEF.CLASSES]), big.mark=','), 100*kk[REEF.CLASSES]/k.med), collapse='  '))
  attr(out, 'k') <- kk
  pmin(pmax(out, 0), 1)
}

#--- area-average a native-resolution raster onto the template -------------------------------
# Averages the native cells whose centres fall in each template cell, ignoring NA. Done by
# hand rather than with resample() so NODATA stays NODATA instead of being interpolated over.
fn.area_average <- function(r, template){
  xy <- raster::xyFromCell(r, seq_len(raster::ncell(r)))
  v  <- raster::getValues(r)
  cell <- raster::cellFromXY(template, xy)
  ok <- !is.na(cell) & !is.na(v)
  tapply(v[ok], cell[ok], mean)
}

#--- rock + grain-size triangle from the RAW dbSeabed grids ----------------------------------
# Values < 0 are the -99 NODATA flag and are kept as NA all the way through, then filled from
# the stratum mean so that "unmeasured" is never silently asserted as "no rock".
fn.dbseabed_layers <- function(dir.raw, template, strata){
  fs <- file.path(dir.raw, paste0('Gmf_', c('RCK', TRIANGLE)),
                  paste0('gmf_', c('RCK', TRIANGLE), '_val.asc'))
  if(any(!file.exists(fs))) stop('raw dbSeabed layer(s) missing under ', dir.raw)
  n <- nrow(strata)
  get1 <- function(f){
    r <- raster::raster(f)
    r <- raster::crop(r, raster::extent(template))
    r[r < 0] <- NA
    agg <- fn.area_average(r, template)
    out <- rep(NA_real_, n)
    j <- match(as.integer(names(agg)), strata$cell)   # land cells return NA and are dropped
    out[j[!is.na(j)]] <- as.numeric(agg)[!is.na(j)]
    out/100
  }
  M <- sapply(fs, get1)
  colnames(M) <- c('RCK', TRIANGLE)
  message('  raw dbSeabed: rock measured in ', sum(!is.na(M[,'RCK'])), ' of ', n, ' water cells (',
          round(100*mean(!is.na(M[,'RCK'])), 1), '%)')
  # Fill NODATA from the stratum mean, but SHRINK that mean toward the domain mean in
  # proportion to how many cells the stratum actually measured. Without the shrinkage a
  # stratum holding three rocky measurements dictates the fill for two hundred unmeasured
  # cells -- on this grid that alone inflated mean rock from 0.102 to 0.145.
  fill <- function(x){
    n.s <- tapply(!is.na(x), strata$stratum, sum)
    m.s <- tapply(x,         strata$stratum, mean, na.rm=TRUE)
    g   <- mean(x, na.rm=TRUE)
    kk  <- stats::median(n.s[n.s > 0])
    t.s <- (n.s * ifelse(is.na(m.s), g, m.s) + kk * g) / (n.s + kk)
    ifelse(is.na(x), t.s[strata$stratum], x)
  }
  for(k in colnames(M)) M[, k] <- fill(M[, k])
  tri <- M[, TRIANGLE, drop=FALSE]
  rs  <- rowSums(tri)
  tri[rs > 0, ] <- tri[rs > 0, , drop=FALSE] / rs[rs > 0]
  tri[rs <= 0, ] <- 1/length(TRIANGLE)
  list(rock = pmin(pmax(M[, 'RCK'], 0), 1), triangle = tri)
}

#--- orchestrator -----------------------------------------------------------------------------
# file.gdb  : GFISHER_EAST_Universe_20XX.gdb
# dir.raw   : .../EnvironmentalDrivers2EwE/data/dbSEABED
# file.sgr  : optional seagrass raster on any grid; area-averaged onto the template.
# split.unconsolidated : FALSE writes one UNC layer; TRUE writes GVL/SND/MUD instead.
# depth.max.reef : reef is forced to 0 below this depth (m). 300 m is the deepest cell in which
#             GFISHER observed any reef class; note only 0.1% of scanning went deeper than 200 m,
#             so the 200-300 m band is filled by IDW from very little evidence. The old pipeline
#             used 200 m.
# k         : shrinkage constant -- 'auto' (per-class, from the variance structure), 'median'
#             (median scanned area), or a fixed value in m2. See fn.shrink_reef. k controls
#             almost entirely the UPPER TAIL of the reef layers, not their means.
fn.make_habitat_basemaps <- function(depth, file.gdb, dir.raw, dir.out, file.sgr=NULL,
                                     k='auto', target='stratum', idp=1.5, nmax=16, maxdist=1,
                                     damp.borrowed=1, n.dbin=5,
                                     split.unconsolidated=FALSE,
                                     depth.max.reef=300, write=TRUE){
  if(!requireNamespace('FNN', quietly=TRUE)) stop("package 'FNN' is required")
  if(inherits(depth, 'SpatRaster')) depth <- raster::raster(depth)
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
  res.min <- round(raster::res(depth)[1]*60, 0)
  dims <- paste0(nrow(depth), 'x', ncol(depth))
  cat('Building habitat basemaps at ', res.min, ' min (', dims, ')\n', sep='')

  cat('Reading GFISHER microgrids...\n');  mgr <- fn.microgrid_reef(file.gdb)
  cat('Aggregating to the model grid...\n'); gr <- fn.reef_to_grid(mgr, depth)
  st <- fn.cell_strata(depth, gr, n.dbin=n.dbin)
  cat('Shrinking reef density toward the ', target, ' target (k = ',
      if(is.numeric(k)) format(round(k), big.mark=',') else k, ')...\n', sep='')
  reef <- fn.shrink_reef(st, gr, damp.borrowed=damp.borrowed, k=k, target=target,
                         idp=idp, nmax=nmax, maxdist=maxdist)
  k.used <- attr(reef, 'k')
  if(is.finite(depth.max.reef)) reef[abs(st$depth) > depth.max.reef, ] <- 0

  cat('Reading raw dbSeabed...\n'); db <- fn.dbseabed_layers(dir.raw, depth, st)

  sgr <- rep(0, nrow(st))
  if(!is.null(file.sgr) && file.exists(file.sgr)){
    a <- fn.area_average(raster::raster(file.sgr), depth)
    j <- match(as.integer(names(a)), st$cell)         # land cells return NA and are dropped
    sgr[j[!is.na(j)]] <- as.numeric(a)[!is.na(j)]
    sgr[is.na(sgr)] <- 0; sgr <- pmin(pmax(sgr, 0), 1)
    cat('  seagrass: mean', round(mean(sgr), 5), '\n')
  }

  # raw fractions, then normalized independently (see header note 4). dbSeabed already
  # partitions the seabed into rock + unconsolidated, so those two carry the bulk of the
  # weight and reef/seagrass dilute them proportionally.
  unc <- 1 - db$rock
  raw <- cbind(reef, RCK = db$rock, UNC = unc, SGR = sgr)
  tot <- rowSums(raw)
  tot[tot <= 0] <- 1
  P <- raw / tot

  if(split.unconsolidated){
    for(k in TRIANGLE) P <- cbind(P, db$triangle[, k] * P[, 'UNC'])
    colnames(P)[(ncol(P)-length(TRIANGLE)+1):ncol(P)] <- TRIANGLE
    P <- P[, setdiff(colnames(P), 'UNC')]
  }

  lyr <- colnames(P)
  out <- list()
  for(k in lyr){
    r <- raster::raster(depth); r[] <- NA_real_
    r[st$cell] <- P[, k]
    names(r) <- k
    out[[k]] <- r
    if(write){
      f <- file.path(dir.out, sprintf('habitat_%s_%dmin_%s.asc', k, res.min, dims))
      raster::writeRaster(r, f, format='ascii', overwrite=TRUE)
    }
  }

  qc <- data.frame(layer = lyr,
                   k_used = c(round(k.used[REEF.CLASSES]), rep(NA, length(lyr)-length(REEF.CLASSES))),
                   mean_prop = round(colMeans(P), 5),
                   max_prop  = round(apply(P, 2, max), 4),
                   cells_gt_10pct = as.integer(colSums(P > 0.10)),
                   stringsAsFactors = FALSE, row.names = NULL)
  s <- rowSums(P)
  cat('\n--- QC ---\n'); print(qc, row.names=FALSE)
  cat('\nsum over layers: min', round(min(s), 6), ' max', round(max(s), 6),
      ' (should be 1)\n')
  cat('water cells:', nrow(st), '  with GFISHER mapping:', sum(st$cell %in% gr$cell), '\n')
  if(write){
    write.csv(qc, file.path(dir.out, sprintf('habitat_basemap_QC_%dmin.csv', res.min)),
              row.names=FALSE)
    cat('\nWrote', length(lyr), 'layers to', dir.out, '\n')
  }
  invisible(list(layers = raster::stack(out), P = P, strata = st, grid.reef = gr, qc = qc))
}

#--- plot a finished basemap set ---------------------------------------------------------------
# One panel per layer, each stretched to its own maximum -- without that the reef classes are
# invisible next to UNC, which occupies ~86% of the average cell. Reads the .asc files written
# by fn.make_habitat_basemaps rather than taking the stack, so it can be re-run on saved output.
# (The old fn.plot_GFISHER_habitats in R/legacy/ expects the superseded 'GFISHER_<CLS>_prop_*'
# filenames and a microgrid footprint raster, and does not work here.)
fn.plot_habitat_basemaps <- function(dir.maps, file.png=NULL, order=NULL){
  fs <- list.files(dir.maps, pattern='^habitat_.*\\.asc$', full.names=TRUE)
  fs <- fs[!grepl('aux\\.xml$', fs)]
  if(!length(fs)) stop('no habitat_*.asc layers found in ', dir.maps)
  code <- sub('^habitat_([A-Z]+)_.*$', '\\1', basename(fs))
  if(is.null(order)) order <- c(REEF.CLASSES, 'RCK', 'UNC', TRIANGLE, 'SGR')
  o <- order(match(code, order), na.last=TRUE)
  fs <- fs[o]; code <- code[o]
  res.min <- round(raster::res(raster::raster(fs[1]))[1]*60, 0)
  if(is.null(file.png)) file.png <- file.path(dir.maps, paste0('habitat_basemaps_', res.min, 'min.png'))
  n <- length(fs); nc <- ceiling(sqrt(n)); nr <- ceiling(n/nc)
  png(file.png, width=500*nc, height=500*nr, res=110)
  op <- par(mfrow=c(nr, nc), mar=c(2,2,3,5))
  for(i in seq_along(fs)){
    r <- raster::raster(fs[i]); v <- raster::getValues(r)
    zmax <- max(v, na.rm=TRUE); if(!is.finite(zmax) || zmax <= 0) zmax <- 1
    raster::plot(r, col=colorRamps::matlab.like(50), colNA='grey85', zlim=c(0, zmax),
                 main=sprintf('%s   mean=%.5f  max=%.3f', code[i], mean(v, na.rm=TRUE), zmax))
  }
  par(op); dev.off()
  message('  wrote ', file.png)
  invisible(file.png)
}
