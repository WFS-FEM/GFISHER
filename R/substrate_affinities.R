# substrate_affinities.R ----------------------------------------------------------
# Substrate affinities for the four dbSeabed classes (RCK, GVL, SND, MUD), estimated
# ENTIRELY INDEPENDENTLY of the reef-type analysis in R/site_level_affinities.R.
# Nothing here conditions on HAB_STRAT, and affinities are normalized only within the
# substrate families below, so the two tables can be produced and defended separately.
#
# ---------------------------------------------------------------------------------
# WHAT dbSEABED ACTUALLY IS -- and why this module reads the RAW grids
# ---------------------------------------------------------------------------------
# The four dbSeabed layers are NOT a single 4-part composition. Measured on the raw
# grids over the WFS domain:
#
#   GVL + SND + MUD  sums to 99-100 on its own (median 99)  <- a closed grain-size triangle
#   RCK              is a separate hardness axis, 0-100, essentially uncorrelated with the
#                    triangle (r = -0.11, -0.05, +0.14) and still leaving the triangle at 99
#                    even where RCK > 50.
#
# The pre-made 5-min layers in maps/dbSEABED/ renormalize all FOUR to sum to 1, which
# compresses rock and deflates the triangle wherever rock occurs. That manufactured closure
# forces the availability-weighted mean of every selection ratio to exactly 1 and leaves no
# room for a ratio to move. Worse, scripts/04_make_ascii_dbSEABED.R maps the -99 NODATA flag
# to 0 before resampling, so RCK's 27,797 genuinely-observed cells become ~90,000 cells of
# mostly fabricated zeros.
#
# So this module reads the RAW dbSEABED grids (EnvironmentalDrivers2EwE/data/dbSEABED/
# Gmf_<CLS>/gmf_<CLS>_val.asc), treats values < 0 as NA rather than 0, and keeps the two
# families apart:
#
#   'grainsize' GVL / SND / MUD, renormalized to sum to 1. A closed composition, so the
#               availability-weighted mean of w is 1 by construction and A ranks the three
#               against one another.
#   'rock'      RCK / 100 as its own axis, scored against its complement (1 - RCK). A_RCK is
#               therefore "how much more rock-covered are the sites this group occupies than
#               the sites the survey visited", NOT a rank against gravel/sand/mud.
#
# Native resolution is 0.02 deg (1.2 arc-min) versus the 5-min grid, which matters far more
# than it sounds: at 5 min the 14,565 stations collapse onto 1,131 cells and only 13 stations
# (in ONE cell) sit on >=50% rock. At native resolution they occupy 3,047 cells and 914
# stations (in 162 cells) sit on >=50% rock. The rock gradient only exists at native scale.
# The cost is coverage: ~75% of stations fall on a cell with valid data on all four layers.
#
# ---------------------------------------------------------------------------------
# ESTIMATOR: use-vs-availability SELECTION RATIO, evaluated at the GFISHER sites.
#   H_s      = substrate at station s (triangle renormalized; rock as its own pair)
#   avail_h  = mean H_h over the stations in the frame     (what the survey encountered)
#   used_h   = sum(MaxN_s * H_sh) / sum(MaxN_s)            (MaxN-weighted use)
#   w_h      = used_h / avail_h
#   A_h      = w_h / max(w_h within family)
#
# WHY the station frame: the survey is reef-targeted and that shows in substrate too, so
# scoring against the whole grid would recover the survey design rather than fish preference.
# Scoring against the stations actually visited asks the answerable question: among the sites
# the cameras went to, did this group concentrate on the rockier or the muddier ones?
#
# Two depth/region adjustments are reported side by side, as in the reef module:
#   'strat' : w computed within each SPACE_STRAT, pooled with W_k = the stratum's share of
#             the mapped microgrid universe.
#   'depth' : w computed within depth-quantile bins, pooled equally, so every class is scored
#             over a common depth distribution.

suppressPackageStartupMessages({library('raster'); library('sp')})

GRAINSIZE.CODES <- c('GVL','SND','MUD')   # closed triangle
ROCK.CODE       <- 'RCK'                  # independent hardness axis
SUBSTRATE.CODES <- c(ROCK.CODE, GRAINSIZE.CODES)

#--- load the RAW dbSeabed grids ------------------------------------------------------------
# dir.raw : .../EnvironmentalDrivers2EwE/data/dbSEABED  (subfolders Gmf_RCK, Gmf_GVL, ...)
# Values < 0 are the -99 NODATA flag and become NA. Optionally cropped to `ext`.
fn.load_dbseabed_raw <- function(dir.raw, codes=SUBSTRATE.CODES, ext=NULL){
  fs <- file.path(dir.raw, paste0('Gmf_', codes), paste0('gmf_', codes, '_val.asc'))
  miss <- fs[!file.exists(fs)]
  if(length(miss)) stop('raw dbSeabed layer(s) not found:\n  ', paste(miss, collapse='\n  '))
  s <- stack(lapply(fs, raster)); names(s) <- codes
  if(!is.null(ext)) s <- crop(s, ext)
  s <- reclassify(s, cbind(-Inf, 0 - .Machine$double.eps, NA))   # -99 NODATA -> NA, keep 0
  message('  raw dbSeabed: ', ncol(s), 'x', nrow(s), ' cells at ',
          round(res(s)[1]*60, 2), ' arc-min')
  s
}

#--- pre-made 5-min layers, for comparison only ---------------------------------------------
fn.load_dbseabed_5min <- function(dir.db, res.min=5, kind='prop', codes=SUBSTRATE.CODES){
  fs <- file.path(dir.db, sprintf('gmf_%s_%s_%dmin.asc', codes, kind, res.min))
  if(any(!file.exists(fs))) stop('layer(s) not found in ', dir.db)
  s <- stack(fs); names(s) <- codes; s
}

#--- query at station coordinates, then split into the two families -------------------------
# Returns list(tri = [station x 3] renormalized grain-size, rock = numeric proportion 0-1,
#              ok = logical, which stations carry valid data on everything).
fn.dbseabed_at_stations <- function(db.stack, lon, lat, scale=100){
  pts <- sp::SpatialPoints(cbind(lon, lat),
                           proj4string=sp::CRS('+proj=longlat +datum=WGS84 +no_defs'))
  if(!is.na(raster::crs(db.stack))) pts <- sp::spTransform(pts, raster::crs(db.stack))
  X <- raster::extract(db.stack, pts)
  colnames(X) <- names(db.stack)
  tri <- X[, GRAINSIZE.CODES, drop=FALSE] / scale
  rs  <- rowSums(tri)
  good.tri <- !is.na(rs) & rs > 0
  tri[good.tri, ] <- tri[good.tri, , drop=FALSE] / rs[good.tri]   # close the triangle to 1
  tri[!good.tri, ] <- NA_real_
  rock <- X[, ROCK.CODE] / scale
  rock[!is.na(rock) & rock > 1] <- 1
  list(tri=tri, rock=rock, ok = good.tri & !is.na(rock))
}

#--- coverage diagnostic --------------------------------------------------------------------
fn.substrate_coverage <- function(tri, rock, db.stack, q=0.5){
  G <- getValues(db.stack)
  gt <- G[, GRAINSIZE.CODES, drop=FALSE] / 100
  rs <- rowSums(gt); keep <- !is.na(rs) & rs > 0
  gt <- gt[keep, , drop=FALSE] / rs[keep]
  gr <- G[, ROCK.CODE] / 100; gr <- gr[!is.na(gr)]
  a.st <- c(colMeans(tri, na.rm=TRUE), RCK=mean(rock, na.rm=TRUE))
  a.gr <- c(colMeans(gt),              RCK=mean(gr))
  cn <- c(GRAINSIZE.CODES, ROCK.CODE)
  tg <- a.st[cn] / a.gr[cn]
  flag <- function(r) if(is.na(r)) 'none' else if(r >= 0.7) 'good' else
                      if(r >= 0.3) 'partial' else 'poor'
  H <- cbind(tri, RCK=rock)
  data.frame(code = cn,
             family = ifelse(cn == ROCK.CODE, 'rock', 'grainsize'),
             mean_at_stations = round(a.st[cn], 4),
             mean_over_grid   = round(a.gr[cn], 4),
             targeting        = round(tg, 3),
             n_stations_gt    = as.integer(colSums(H[, cn, drop=FALSE] >= q, na.rm=TRUE)),
             coverage         = sapply(tg, flag),
             row.names = NULL, stringsAsFactors = FALSE)
}

#--- selection ratio within groups of stations, pooled with weights -------------------------
# Groups where the model group has zero total MaxN carry no use information, so they drop out
# of BOTH the used and the availability pool and the weights renormalize -- otherwise used and
# avail would be built over different station sets and w would be biased.
fn.substrate_w <- function(H, E, grp, W, min.n=10){
  lev <- intersect(names(W), levels(droplevels(as.factor(grp))))
  num <- rep(0, ncol(H)); den <- rep(0, ncol(H)); tot <- 0
  for(k in lev){
    i <- which(grp == k)
    if(length(i) < min.n || sum(E[i]) <= 0) next
    wk  <- W[[k]]
    num <- num + wk * (colSums(E[i] * H[i, , drop=FALSE]) / sum(E[i]))
    den <- den + wk * colMeans(H[i, , drop=FALSE])
    tot <- tot + wk
  }
  if(tot <= 0) return(setNames(rep(NA_real_, ncol(H)), colnames(H)))
  w <- (num / tot) / (den / tot)
  w[den == 0] <- NA_real_
  setNames(w, colnames(H))
}

#--- ratios for both families, returned as one named vector over c(GVL,SND,MUD,RCK) ---------
# grain size: the closed triangle, A normalized across the three.
# rock:       scored against its complement, so A_RCK = w_RCK / max(w_RCK, w_notRCK).
fn.substrate_ratios <- function(tri, rock, E, grp, W, min.n=10){
  w.t <- fn.substrate_w(tri, E, grp, W, min.n)
  w.r <- fn.substrate_w(cbind(RCK=rock, notRCK=1-rock), E, grp, W, min.n)
  A.t <- if(all(is.na(w.t))) w.t else w.t / max(w.t, na.rm=TRUE)
  A.r <- if(all(is.na(w.r))) NA_real_ else unname(w.r['RCK'] / max(w.r, na.rm=TRUE))
  list(w = c(w.t, RCK=unname(w.r['RCK'])),
       A = c(A.t, RCK=A.r))
}

#--- bootstrap over stations, within stratum ------------------------------------------------
fn.boot_substrate <- function(tri, rock, E, strat, dbin, Ws, Wd, n.boot=1000, seed=1, min.n=10){
  set.seed(seed)
  cn <- c(GRAINSIZE.CODES, ROCK.CODE)
  idx <- split(seq_along(E), strat)
  Bs <- matrix(NA_real_, n.boot, length(cn), dimnames=list(NULL, cn))
  Bd <- Bs
  for(b in seq_len(n.boot)){
    i <- unlist(lapply(idx, function(ii) if(length(ii)) sample(ii, length(ii), TRUE) else ii),
                use.names=FALSE)
    Bs[b, ] <- fn.substrate_ratios(tri[i, , drop=FALSE], rock[i], E[i], strat[i], Ws, min.n)$A[cn]
    Bd[b, ] <- fn.substrate_ratios(tri[i, , drop=FALSE], rock[i], E[i], dbin[i],  Wd, min.n)$A[cn]
  }
  qs <- function(B) t(apply(B, 2, stats::quantile, c(0.025, 0.975), na.rm=TRUE))
  list(strat=qs(Bs), depth=qs(Bd))
}

#--- batch every model group ----------------------------------------------------------------
# Requires R/site_level_affinities.R to have been sourced: the station frame, group matrix and
# stanza pooling are shared with the reef analysis. The DATA is shared; the ESTIMATE is not.
fn.batch_substrate_affinities <- function(maxn.long, file.env, bbox, dir.raw, dir.out,
                                          W=NULL, file.gdb=NULL, n.boot=1000, seed=1,
                                          n.dbin=5, min.n=10, min.pos=30, pool.stanzas=TRUE,
                                          tag=''){
  for(f in c('fn.site_stations','fn.site_group_matrix','fn.pool_site_stanzas','fn.stratum_weights'))
    if(!exists(f)) stop('source R/site_level_affinities.R first (missing ', f, ')')
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)

  cat('Building station frame...\n')
  st <- fn.site_stations(file.env, bbox)
  cat('Loading RAW dbSeabed...\n')
  db <- fn.load_dbseabed_raw(dir.raw,
          ext=raster::extent(bbox[['lonW']], bbox[['lonE']], bbox[['latS']], bbox[['latN']]))
  q  <- fn.dbseabed_at_stations(db, st$lon, st$lat)
  cat('  ', sum(q$ok), ' of ', nrow(st), ' stations carry valid data on all four layers',
      sprintf(' (%.1f%%)\n', 100*mean(q$ok)), sep='')
  st <- st[q$ok, ]; tri <- q$tri[q$ok, , drop=FALSE]; rock <- q$rock[q$ok]

  cov <- fn.substrate_coverage(tri, rock, db)
  cat('\nSubstrate coverage (targeting = station mean / grid mean):\n')
  print(cov, row.names=FALSE)

  if(is.null(W)){
    if(!is.null(file.gdb) && file.exists(file.gdb)) W <- fn.stratum_weights(file.gdb)
    else { message('  no gdb supplied -- station-count stratum weights')
           W <- prop.table(table(st$strat)) }
  }
  Ws <- W[intersect(names(W), unique(st$strat))]
  br <- unique(stats::quantile(st$depth, seq(0, 1, length.out=n.dbin+1), na.rm=TRUE))
  st$dbin <- factor(paste0('d', cut(st$depth, br, include.lowest=TRUE, labels=FALSE)))
  Wd <- setNames(rep(1/nlevels(st$dbin), nlevels(st$dbin)), levels(st$dbin))
  cat('\ndepth bins (m):', paste(round(br,1), collapse=' | '), '\n')

  gm <- fn.site_group_matrix(maxn.long, st)
  if(isTRUE(pool.stanzas)){
    cat('Pooling sparse stanzas...\n')
    gm <- fn.pool_site_stanzas(gm$y, gm$groups, min.pos=min.pos, use=rep(TRUE, nrow(st)))
  }
  y.all <- gm$y; groups <- gm$groups
  cn <- c(GRAINSIZE.CODES, ROCK.CODE)

  long <- list(); panels <- list()
  for(g in seq_len(nrow(groups))){
    E <- y.all[, g]
    cat(sprintf('  mod%-3d %-22s n_pos=%d\n', groups$modnumber[g], groups$modname[g], sum(E > 0)))
    rs <- fn.substrate_ratios(tri, rock, E, st$strat, Ws, min.n)
    rd <- fn.substrate_ratios(tri, rock, E, st$dbin,  Wd, min.n)
    ci <- fn.boot_substrate(tri, rock, E, st$strat, st$dbin, Ws, Wd, n.boot, seed, min.n)
    row <- data.frame(
      modnumber = groups$modnumber[g], modname = groups$modname[g], code = cn,
      family    = ifelse(cn == ROCK.CODE, 'rock', 'grainsize'),
      n_pos     = sum(E > 0),
      avail     = round(c(colMeans(tri), RCK=mean(rock))[cn], 5),
      w_strat   = round(as.numeric(rs$w[cn]), 3),
      A_strat   = round(as.numeric(rs$A[cn]), 3),
      A_strat_lo= round(ci$strat[cn, 1], 3), A_strat_hi = round(ci$strat[cn, 2], 3),
      w_depth   = round(as.numeric(rd$w[cn]), 3),
      A_depth   = round(as.numeric(rd$A[cn]), 3),
      A_depth_lo= round(ci$depth[cn, 1], 3), A_depth_hi = round(ci$depth[cn, 2], 3),
      coverage  = cov$coverage[match(cn, cov$code)],
      pooled_with = groups$pooled_with[g], stringsAsFactors = FALSE)
    long[[g]] <- row; panels[[g]] <- row
  }
  long <- do.call(rbind, long)
  ord <- c(ROCK.CODE, GRAINSIZE.CODES)   # report as RCK, GVL, SND, MUD

  mk.wide <- function(col){
    w <- reshape(long[, c('modnumber','modname','code',col)],
                 idvar=c('modnumber','modname'), timevar='code', direction='wide')
    names(w) <- sub(paste0('^', col, '\\.'), '', names(w))
    w <- w[, c('modnumber','modname', ord[ord %in% names(w)])]
    w[order(w$modnumber), ]
  }
  wide.s <- mk.wide('A_strat'); wide.d <- mk.wide('A_depth')

  sfx <- if(nzchar(tag)) paste0('_', tag) else ''
  f.long <- file.path(dir.out, paste0('substrate_affinity_long', sfx, '.csv'))
  f.ws   <- file.path(dir.out, paste0('substrate_affinity_A_wide_strat', sfx, '.csv'))
  f.wd   <- file.path(dir.out, paste0('substrate_affinity_A_wide_depth', sfx, '.csv'))
  f.cov  <- file.path(dir.out, paste0('substrate_coverage', sfx, '.csv'))
  f.pdf  <- file.path(dir.out, paste0('substrate_affinity_fits', sfx, '.pdf'))
  write.csv(long,   f.long, row.names=FALSE)
  write.csv(wide.s, f.ws,   row.names=FALSE, quote=FALSE)
  write.csv(wide.d, f.wd,   row.names=FALSE, quote=FALSE)
  write.csv(cov,    f.cov,  row.names=FALSE)

  pdf(f.pdf, width=10, height=7.5, onefile=TRUE)
  op <- par(mfrow=c(2,2), mar=c(4,4,3,1))
  for(p in panels){
    p <- p[match(ord, p$code), ]
    if(all(is.na(p$A_strat)) && all(is.na(p$A_depth))){
      plot.new(); title(paste0(p$modname[1], '\n(no fit)')); next
    }
    M <- rbind(strat=p$A_strat, depth=p$A_depth); M[is.na(M)] <- 0
    yhi <- max(1.05, p$A_strat_hi, p$A_depth_hi, na.rm=TRUE)
    bp <- barplot(M, beside=TRUE, names.arg=p$code, las=2, ylim=c(0, yhi),
                  col=c('steelblue4','darkorange3'), border=NA,
                  ylab='affinity  (selection ratio, rescaled)',
                  main=sprintf('mod%d  %s\n(%d stations with MaxN>0)',
                               p$modnumber[1], p$modname[1], p$n_pos[1]))
    suppressWarnings({
      arrows(bp[1, ], p$A_strat_lo, bp[1, ], p$A_strat_hi, angle=90, code=3, length=0.02)
      arrows(bp[2, ], p$A_depth_lo, bp[2, ], p$A_depth_hi, angle=90, code=3, length=0.02)
    })
    abline(v=(bp[2,1]+bp[1,2])/2, lty=3, col='grey40')   # rock | grain-size divider
    box(bty='l')
  }
  par(op)
  plot.new()
  legend('center', bty='n', fill=c('steelblue4','darkorange3'),
         legend=c('stratified by SPACE_STRAT', 'standardized over depth bins'),
         title='adjustment')
  legend('bottom', bty='n', lty=c(3,NA), col=c('grey40',NA),
         legend=c('divider: RCK is an independent hardness axis, scored against its complement',
                  'GVL/SND/MUD are a closed triangle - A ranks them against one another only'))
  dev.off()

  cat('\nWrote:\n  ', f.long, '\n  ', f.ws, '\n  ', f.wd, '\n  ', f.cov, '\n  ', f.pdf, '\n', sep='')
  invisible(list(long=long, wide.strat=wide.s, wide.depth=wide.d, coverage=cov,
                 stations=st, tri=tri, rock=rock))
}
