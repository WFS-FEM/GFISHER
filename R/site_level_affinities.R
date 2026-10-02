# site_level_affinities.R ---------------------------------------------------------
# Habitat affinities estimated at the SITE (video station) level, from the paired GFISHER
# camera + habitat survey, instead of from 5-min raster cells.
#
# WHY a separate module: the raster-based selection ratios in R/selection_ratio_affinities.R
# compare MaxN-weighted habitat use against habitat availability averaged over ~9 km Ecospace
# cells. Three problems come with that: (1) the camera saw one ~200 m spot, not a 9 km average;
# (2) outside the mapped microgrid footprint the reef layers are IDW guesses (natural) or
# fabricated zeros (artificial), so part of the "availability" was never measured; and (3) the
# survey deliberately targeted geoforms, so use/availability partly recovers the survey design
# rather than fish preference -- high-relief natural reef is 0.85% of the mapped universe by
# area but 6.1% of an effort-weighted surveyed cell, a 7x enrichment.
#
# The station data already carry the pairing, in two independent families:
#
#   'reef'      HAB_STRAT -- the habitat class read AT the site, in the SAME taxonomy as the
#               side-scan mapping ({N,A} x {L,M,H} x {S,M,L}). Collapsed to the first two
#               characters it is exactly the six Ecospace classes (AL,AM,AH,NL,NM,NH).
#   'substrate' ARTI_PER / ROCK_PER / SED_PER -- percent cover scored from the video itself.
#               These three close to 100 (median 100, 9657 of 9665 WFS stations within 95-105)
#               and are available on 9665 stations, 9376 of which also carry HAB_STRAT, so the
#               two families describe essentially the same set of sites. Each station is
#               assigned its dominant substrate (ART/RCK/SED).
#
# NOTE on substrate resolution: the finer in-situ reads (SHELLGRAVEL_PER, SILTSANDCLAY_PER)
# exist on only 3215 WFS stations, of which just 154 also carry HAB_STRAT, so gravel/sand/mud
# cannot be separated here. The in-situ family is hard (RCK) vs soft (SED) vs artificial (ART).
# If Ecospace needs the GVL/SND/MUD split, that ratio still has to come from dbSeabed.
#
# ESTIMATOR: because availability is controlled by the sampler in a targeted stratified design,
# this module does NOT form a use/availability ratio. It estimates the RESPONSE per station in
# each habitat class -- mean MaxN with zeros retained, which is relative habitat capacity, the
# quantity Ecospace actually wants -- and normalizes it to [0,1] WITHIN EACH FAMILY, so A=1
# marks the best reef class and, separately, the best substrate class. Two adjustments are
# produced side by side because class is badly confounded with depth (mean station depth is
# 75.8 m in NM vs 26.8 m in AM):
#
#   'strat' : design-based. Mean MaxN within each (class x SPACE_STRAT) cell, pooled across
#             strata with weights W_k = the stratum's share of the mapped microgrid universe.
#             Non-parametric; costs precision in thin cells.
#   'depth' : model-based. quasi-Poisson GLM  maxn ~ class + poly(depth, 2), then G-computation
#             (predict every station as if it were each class in turn and average), so each
#             class is evaluated over the SAME depth distribution. Uses all stations; assumes
#             the depth response is common across classes.
#
# Agreement between the two is the signal that a class ranking is real rather than a depth
# artifact. Both are bootstrapped over stations (the actual independent unit, ~10k of them,
# versus 1,131 raster cells).

suppressPackageStartupMessages(library('stats'))

HAB.CLASSES <- c('AL','AM','AH','NL','NM','NH')   # family 'reef'      (from HAB_STRAT)
SED.CLASSES <- c('ART','RCK','SED')               # family 'substrate' (from *_PER cover)

#--- station frame: one row per video station, with both class families, depth and stratum --
# file.env : env3LABS_*.csv (station table). Column names are lowercased on read.
# bbox     : c(latN, latS, lonW, lonE), as in the main driver.
# min.dom  : a station is assigned a dominant substrate only if that substrate covers at least
#            this percent; below it the substrate read is treated as ambiguous (NA).
# Returns a data.frame: reference, lon, lat, depth, strat, cls (reef class or NA),
# sed (substrate class or NA). Rows lacking a stratum or a depth are dropped, since every
# estimator here needs both.
fn.site_stations <- function(file.env, bbox, min.dom=50,
                             lon.col='lon_dd', lat.col='lat_dd', id.col='reference',
                             class.col='hab_strat', strat.col='space_strat', depth.col='depth',
                             sed.cols=c(ART='arti_per', RCK='rock_per', SED='sed_per')){
  e <- read.csv(file.env, header=TRUE, stringsAsFactors=FALSE)
  names(e) <- tolower(names(e))
  need <- c(id.col, lon.col, lat.col, class.col, strat.col, depth.col)
  if(any(!need %in% names(e)))
    stop('env file missing column(s): ', paste(need[!need %in% names(e)], collapse=', '))
  miss.sed <- sed.cols[!sed.cols %in% names(e)]
  if(length(miss.sed)) message('  substrate column(s) not found, family will be empty: ',
                               paste(miss.sed, collapse=', '))

  num <- function(x) suppressWarnings(as.numeric(x))
  d <- data.frame(reference = e[[id.col]],
                  lon = num(e[[lon.col]]), lat = num(e[[lat.col]]),
                  depth = num(e[[depth.col]]),
                  strat = trimws(e[[strat.col]]),
                  cls = substr(trimws(e[[class.col]]), 1, 2),   # NLL/NHM/AMS -> 6 classes
                  stringsAsFactors = FALSE)

  # dominant in-situ substrate, from the three percent-cover fields that close to 100
  P <- sapply(sed.cols, function(k) if(k %in% names(e)) num(e[[k]]) else rep(NA_real_, nrow(e)))
  P <- matrix(P, nrow=nrow(e), dimnames=list(NULL, names(sed.cols)))
  full <- stats::complete.cases(P)
  d$sed <- NA_character_
  if(any(full)){
    win <- max.col(P[full, , drop=FALSE], ties.method='first')
    top <- P[full, , drop=FALSE][cbind(seq_len(sum(full)), win)]
    lab <- colnames(P)[win]
    lab[!is.na(top) & top < min.dom] <- NA_character_   # no clear dominant cover
    d$sed[full] <- lab
  }

  d <- d[!is.na(d$lon) & !is.na(d$lat), ]
  d <- d[d$lat >= bbox[2] & d$lat <= bbox[1] & d$lon >= bbox[3] & d$lon <= bbox[4], ]
  d <- d[!duplicated(d$reference), ]
  d <- d[nzchar(d$strat) & !is.na(d$depth), ]

  d$cls[!d$cls %in% HAB.CLASSES] <- NA_character_
  d$sed[!d$sed %in% SED.CLASSES] <- NA_character_
  d
}

#--- station x model-group MaxN matrix, WITH the zeros put back ----------------------------
# maxn.long : the data.frame from fn.make_gfisher_videodataset() -- presence records only
#             (its melt step drops maxn==0), carrying reference / modnumber / modname / maxn.
# A station absent from maxn.long for a group genuinely observed zero of it, so the matrix is
# filled with 0 rather than NA; that zero is what makes a mean MaxN a density instead of a
# conditional-on-presence mean.
fn.site_group_matrix <- function(maxn.long, stations){
  need <- c('reference','modnumber','maxn')
  if(any(!need %in% names(maxn.long)))
    stop('maxn.long missing column(s): ', paste(need[!need %in% names(maxn.long)], collapse=', '))
  d <- maxn.long[!is.na(maxn.long$modnumber) & !is.na(maxn.long$maxn), ]
  d$modnumber <- as.integer(as.character(d$modnumber))
  d <- d[!is.na(d$modnumber), ]

  groups <- unique(d[, c('modnumber','modname')])
  groups <- groups[order(groups$modnumber), ]
  rownames(groups) <- NULL
  groups$pooled_with <- NA_character_

  d <- d[d$reference %in% stations$reference, ]
  agg <- aggregate(maxn ~ reference + modnumber, data=d, FUN=sum)

  y <- matrix(0, nrow(stations), nrow(groups),
              dimnames=list(stations$reference, paste0('mod', groups$modnumber)))
  i <- match(agg$reference, stations$reference)
  j <- match(agg$modnumber, groups$modnumber)
  ok <- !is.na(i) & !is.na(j)
  y[cbind(i[ok], j[ok])] <- agg$maxn[ok]
  list(y=y, groups=groups)
}

#--- split a modname into species base + stanza index --------------------------------------
# Multi-stanza groups arrive from the species list as '<species> <k>', with a trailing '+' on
# the plus group ('gag 0' ... 'gag 5+'). Single-stanza groups ('sharks') get stanza NA and are
# never pooled. (The raster module sees the same groups after filename sanitizing, as
# 'gag-5-', and parses them with its own fn.parse_stanza.)
fn.parse_stanza_site <- function(modname){
  has <- grepl(' [0-9]+[+]?$', modname)
  data.frame(base   = ifelse(has, sub(' [0-9]+[+]?$', '', modname), modname),
             stanza = suppressWarnings(as.integer(
                        ifelse(has, sub('.* ([0-9]+)[+]?$', '\\1', modname), NA))),
             stringsAsFactors = FALSE)
}

#--- pool sparse stanzas with their nearest same-species neighbours -------------------------
# A stanza seen at fewer than min.pos stations cannot support a density comparison across six
# classes -- red grouper 1 (4 positive stations) returned a degenerate 0/0/0/0.255/0/1 fit.
# Rather than leave that to be patched downstream, merge its per-station MaxN with the nearest
# stanza(s) of the SAME species (adding one neighbour at a time, nearest stanza first) until
# the pooled column clears min.pos, and estimate from that. Neighbours keep their own
# independent fits; only the sparse stanza's column changes, and `pooled_with` records what
# went into it.
#
# `use` restricts the positive-station count to the stations the estimator will actually run
# on. It must be the REEF subset: that family has to separate six classes rather than three,
# so it is the binding constraint, and counting over the whole frame lets a group slip through
# (gag 0 has 34 positive stations frame-wide but only 19 carrying a HAB_STRAT). Pooling is
# decided once and both families inherit it, so they describe the same pooled group.
fn.pool_site_stanzas <- function(y, groups, min.pos=30, use=NULL){
  nm  <- groups$modname
  stz <- fn.parse_stanza_site(nm)
  if(is.null(use)) use <- rep(TRUE, nrow(y))
  npos <- colSums(y[use, , drop=FALSE] > 0)
  for(i in which(npos < min.pos & !is.na(stz$stanza))){
    cand <- setdiff(which(stz$base == stz$base[i] & !is.na(stz$stanza)), i)
    if(length(cand) == 0) next
    cand <- cand[order(abs(stz$stanza[cand] - stz$stanza[i]), stz$stanza[cand])]
    v <- y[, i]; part <- character(0); np <- npos[i]
    for(j in cand){
      v    <- v + y[, j]
      part <- c(part, nm[j])
      np   <- sum(v[use] > 0)
      if(np >= min.pos) break
    }
    message(sprintf("    mod%d %s: only %d station(s) with MaxN>0 -- pooled with %s (now %d)",
                    groups$modnumber[i], nm[i], npos[i], paste(part, collapse=' + '), np))
    y[, i] <- v
    groups$pooled_with[i] <- paste(part, collapse=' + ')
  }
  list(y=y, groups=groups)
}

#--- stratum weights from the mapped microgrid universe ------------------------------------
# Counts mapped microgrids per SpaceStrat, so pooling across strata reflects how much of the
# mapped shelf each stratum represents rather than how hard it was surveyed. Reads only the
# attribute table of the microgrid layer.
fn.stratum_weights <- function(file.gdb, layer.pattern='Microgrid', strat.field='SpaceStrat'){
  if(!requireNamespace('sf', quietly=TRUE)) stop("package 'sf' is required")
  lyrs <- sf::st_layers(file.gdb)
  lyr  <- lyrs$name[grep(layer.pattern, lyrs$name, ignore.case=TRUE)]
  if(length(lyr) != 1) stop('could not identify a unique microgrid layer in ', file.gdb)
  m <- sf::st_drop_geometry(sf::st_read(file.gdb, layer=lyr, quiet=TRUE))
  if(!strat.field %in% names(m)) stop(strat.field, ' not found in layer ', lyr)
  w <- table(trimws(m[[strat.field]]))
  w <- w[nzchar(names(w))]
  message('  stratum weights from ', sum(w), ' mapped microgrids across ', length(w), ' strata')
  w / sum(w)
}

#--- estimator 1: design-based, stratified by SPACE_STRAT ----------------------------------
# Mean MaxN per station within each (class x stratum) cell, pooled across strata with weights
# W. Cells holding fewer than min.cell stations are dropped and the weights renormalized over
# the strata that survive FOR THAT CLASS, so a class present in only part of the shelf is
# scored on the part where it occurs rather than being penalized for absent strata.
fn.density_stratified <- function(y, cls, strat, W, min.cell=5){
  lev <- levels(cls)
  n   <- table(cls, strat)
  s   <- tapply(y, list(cls, strat), mean)
  out <- setNames(rep(NA_real_, length(lev)), lev)
  for(h in lev){
    k  <- colnames(s)[!is.na(s[h, ]) & n[h, ] >= min.cell]
    wk <- W[k]; wk <- wk[!is.na(wk)]
    if(!length(wk) || sum(wk) == 0) next
    out[h] <- sum(wk * s[h, names(wk)]) / sum(wk)
  }
  out
}

#--- estimator 2: model-based, depth as a continuous covariate -----------------------------
# quasi-Poisson GLM (log link) of MaxN on class + a quadratic in depth, then G-computation:
# every station is predicted as if it belonged to each class in turn and the predictions are
# averaged, so all classes are compared over the same observed depth distribution.
fn.density_depth <- function(y, cls, depth, df.depth=2){
  na.out <- setNames(rep(NA_real_, nlevels(cls)), levels(cls))
  ok <- !is.na(depth) & is.finite(depth)
  if(sum(ok) < 50 || sum(y[ok]) <= 0) return(na.out)
  d <- data.frame(y=y[ok], cls=droplevels(cls[ok]), depth=depth[ok])
  if(nlevels(d$cls) < 2) return(na.out)
  fit <- try(suppressWarnings(
           glm(y ~ cls + poly(depth, df.depth), data=d, family=quasipoisson())), silent=TRUE)
  if(inherits(fit, 'try-error')) return(na.out)
  out <- na.out
  for(h in levels(d$cls)){
    nd <- d; nd$cls <- factor(h, levels=levels(d$cls))
    p  <- try(predict(fit, newdata=nd, type='response'), silent=TRUE)
    if(!inherits(p, 'try-error')) out[h] <- mean(p, na.rm=TRUE)
  }
  out
}

#--- normalize a density vector to an affinity in [0,1] ------------------------------------
fn.affinity <- function(v){
  mx <- suppressWarnings(max(v, na.rm=TRUE))
  if(!is.finite(mx) || mx <= 0) return(setNames(rep(NA_real_, length(v)), names(v)))
  v / mx
}

#--- stratified bootstrap over stations ----------------------------------------------------
# Resamples stations WITHIN each SPACE_STRAT, preserving the survey's stratum allocation, and
# recomputes both estimators. Returns 2.5%/97.5% percentile bounds on each affinity.
fn.boot_site <- function(y, cls, strat, depth, W, n.boot=500, seed=1, min.cell=5, df.depth=2){
  set.seed(seed)
  lev <- levels(cls); K <- length(lev)
  idx.by.k <- split(seq_along(cls), strat)
  Bs <- matrix(NA_real_, n.boot, K, dimnames=list(NULL, lev))
  Bd <- matrix(NA_real_, n.boot, K, dimnames=list(NULL, lev))
  for(b in seq_len(n.boot)){
    i <- unlist(lapply(idx.by.k, function(ii) if(length(ii)) sample(ii, length(ii), TRUE) else ii),
                use.names=FALSE)
    Bs[b, ] <- fn.affinity(fn.density_stratified(y[i], cls[i], strat[i], W, min.cell))[lev]
    Bd[b, ] <- fn.affinity(fn.density_depth(y[i], cls[i], depth[i], df.depth))[lev]
  }
  qs <- function(B) t(apply(B, 2, stats::quantile, c(0.025, 0.975), na.rm=TRUE))
  list(strat=qs(Bs), depth=qs(Bd))
}

#--- one family (reef or substrate) for one group -------------------------------------------
fn.site_one <- function(yg, st, class.col, W, n.boot, seed, min.cell, df.depth){
  keep <- !is.na(st[[class.col]])
  s <- st[keep, ]; v <- yg[keep]
  s$cf <- factor(s[[class.col]],
                 levels=intersect(if(class.col=='cls') HAB.CLASSES else SED.CLASSES,
                                  unique(s[[class.col]])))
  lev <- levels(s$cf)
  d.s <- fn.density_stratified(v, s$cf, s$strat, W, min.cell)
  d.d <- fn.density_depth(v, s$cf, s$depth, df.depth)
  ci  <- fn.boot_site(v, s$cf, s$strat, s$depth, W, n.boot, seed, min.cell, df.depth)
  data.frame(class      = lev,
             n_stations = as.integer(table(s$cf)[lev]),
             n_pos      = as.integer(tapply(v > 0, s$cf, sum)[lev]),
             mean_maxn  = round(as.numeric(tapply(v, s$cf, mean)[lev]), 4),
             dens_strat = round(as.numeric(d.s[lev]), 4),
             A_strat    = round(as.numeric(fn.affinity(d.s)[lev]), 3),
             A_strat_lo = round(ci$strat[lev, 1], 3),
             A_strat_hi = round(ci$strat[lev, 2], 3),
             dens_depth = round(as.numeric(d.d[lev]), 4),
             A_depth    = round(as.numeric(fn.affinity(d.d)[lev]), 3),
             A_depth_lo = round(ci$depth[lev, 1], 3),
             A_depth_hi = round(ci$depth[lev, 2], 3),
             stringsAsFactors = FALSE)
}

#--- batch every model group, both families -------------------------------------------------
fn.batch_site_affinities <- function(maxn.long, file.env, bbox, dir.out, W=NULL, file.gdb=NULL,
                                     n.boot=500, seed=1, min.cell=5, df.depth=2,
                                     min.class.n=30, min.pos=30, pool.stanzas=TRUE,
                                     min.dom=50, tag=''){
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
  cat('Building station frame...\n')
  st <- fn.site_stations(file.env, bbox, min.dom=min.dom)
  cat('  ', nrow(st), ' stations with a stratum and depth\n', sep='')

  # drop under-represented classes in each family, independently
  for(cc in c('cls','sed')){
    tb <- table(st[[cc]])
    drop <- names(which(tb < min.class.n))
    if(length(drop)){
      message('  dropping ', cc, ' class(es) with <', min.class.n, ' stations: ',
              paste(drop, collapse=', '))
      st[[cc]][st[[cc]] %in% drop] <- NA_character_
    }
  }
  cat('  reef  (HAB_STRAT):     ', sum(!is.na(st$cls)), ' stations\n', sep='')
  print(table(factor(st$cls, levels=HAB.CLASSES)))
  cat('  substrate (in-situ %): ', sum(!is.na(st$sed)), ' stations\n', sep='')
  print(table(factor(st$sed, levels=SED.CLASSES)))
  cat('  carrying both families:', sum(!is.na(st$cls) & !is.na(st$sed)), '\n')

  if(is.null(W)){
    if(!is.null(file.gdb) && file.exists(file.gdb)) W <- fn.stratum_weights(file.gdb)
    else {
      message('  no gdb supplied -- falling back to station-count stratum weights')
      W <- prop.table(table(st$strat))
    }
  }
  miss <- setdiff(unique(st$strat), names(W))
  if(length(miss)) message('  strata with no weight (dropped from pooling): ',
                           paste(miss, collapse=', '))

  gm <- fn.site_group_matrix(maxn.long, st)
  if(isTRUE(pool.stanzas)){
    cat('Pooling sparse stanzas...\n')
    gm <- fn.pool_site_stanzas(gm$y, gm$groups, min.pos=min.pos, use=!is.na(st$cls))
  }
  y.all <- gm$y; groups <- gm$groups

  fams <- c(reef='cls', substrate='sed')
  fams <- fams[sapply(fams, function(cc) any(!is.na(st[[cc]])))]

  long <- list(); panels <- list()
  for(g in seq_len(nrow(groups))){
    cat(sprintf('  mod%-3d %-22s n_pos=%d%s\n', groups$modnumber[g], groups$modname[g],
                sum(y.all[!is.na(st$cls), g] > 0),   # reef-family stations, as used for pooling
                if(is.na(groups$pooled_with[g])) '' else paste0('  [pooled: ', groups$pooled_with[g], ']')))
    rows <- lapply(names(fams), function(fm){
      r <- fn.site_one(y.all[, g], st, fams[[fm]], W, n.boot, seed, min.cell, df.depth)
      cbind(modnumber=groups$modnumber[g], modname=groups$modname[g], family=fm,
            r, pooled_with=groups$pooled_with[g], stringsAsFactors=FALSE)
    })
    r <- do.call(rbind, rows)
    long[[g]] <- r; panels[[g]] <- r
  }
  long <- do.call(rbind, long)
  ord  <- c(HAB.CLASSES, SED.CLASSES)

  mk.wide <- function(col){
    w <- reshape(long[, c('modnumber','modname','class',col)],
                 idvar=c('modnumber','modname'), timevar='class', direction='wide')
    names(w) <- sub(paste0('^', col, '\\.'), '', names(w))
    w <- w[, c('modnumber','modname', ord[ord %in% names(w)])]
    w[order(w$modnumber), ]
  }
  wide.s <- mk.wide('A_strat'); wide.d <- mk.wide('A_depth')

  sfx <- if(nzchar(tag)) paste0('_', tag) else ''
  f.long <- file.path(dir.out, paste0('site_affinity_long', sfx, '.csv'))
  f.ws   <- file.path(dir.out, paste0('site_affinity_A_wide_strat', sfx, '.csv'))
  f.wd   <- file.path(dir.out, paste0('site_affinity_A_wide_depth', sfx, '.csv'))
  f.pdf  <- file.path(dir.out, paste0('site_affinity_fits', sfx, '.pdf'))
  write.csv(long,   f.long, row.names=FALSE)
  write.csv(wide.s, f.ws,   row.names=FALSE, quote=FALSE)
  write.csv(wide.d, f.wd,   row.names=FALSE, quote=FALSE)

  # one panel per group: the two adjustments side by side, with bootstrap CIs and a divider
  # between the families (A is normalized within each). Where the two adjustments disagree,
  # the class ranking is being driven by depth rather than by habitat.
  pdf(f.pdf, width=10, height=7.5, onefile=TRUE)
  op <- par(mfrow=c(2,2), mar=c(4,4,3,1))
  for(p in panels){
    p <- p[order(match(p$class, ord)), ]
    if(all(is.na(p$A_strat)) && all(is.na(p$A_depth))){
      plot.new(); title(paste0(p$modname[1], '\n(no fit)')); next
    }
    M <- rbind(strat=p$A_strat, depth=p$A_depth); M[is.na(M)] <- 0
    yhi <- max(1.05, p$A_strat_hi, p$A_depth_hi, na.rm=TRUE)
    sub <- if(is.na(p$pooled_with[1])) '' else paste0('  [pooled: ', p$pooled_with[1], ']')
    bp  <- barplot(M, beside=TRUE, names.arg=p$class, las=2, ylim=c(0, yhi),
                   col=c('steelblue4','darkorange3'), border=NA,
                   ylab='affinity  (relative density)',
                   main=sprintf('mod%d  %s%s\n(%d stations with MaxN>0)',
                                p$modnumber[1], p$modname[1], sub,
                                sum(p$n_pos[p$family=='reef'], na.rm=TRUE)))
    suppressWarnings({
      arrows(bp[1, ], p$A_strat_lo, bp[1, ], p$A_strat_hi, angle=90, code=3, length=0.02)
      arrows(bp[2, ], p$A_depth_lo, bp[2, ], p$A_depth_hi, angle=90, code=3, length=0.02)
    })
    fb <- which(diff(as.integer(factor(p$family, levels=unique(p$family)))) != 0)
    if(length(fb)) abline(v=(bp[2, fb] + bp[1, fb+1])/2, lty=3, col='grey40')
    box(bty='l')
  }
  par(op)
  plot.new()
  legend('center', bty='n', fill=c('steelblue4','darkorange3'),
         legend=c('stratified by SPACE_STRAT (design-based)',
                  'depth as continuous covariate (quasi-Poisson GLM, G-computation)'),
         title='depth/region adjustment')
  legend('bottom', bty='n', lty=c(NA,3), col=c(NA,'grey40'),
         legend=c('bars agreeing across both adjustments = habitat signal; disagreeing = depth artifact',
                  'family divider - affinity is normalized within reef and within substrate'))
  dev.off()

  cat('\nWrote:\n  ', f.long, '\n  ', f.ws, '\n  ', f.wd, '\n  ', f.pdf, '\n', sep='')
  invisible(list(long=long, wide.strat=wide.s, wide.depth=wide.d, stations=st, W=W, groups=groups))
}
