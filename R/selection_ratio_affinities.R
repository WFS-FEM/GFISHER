# selection_ratio_affinities.R --------------------------------------------------
# Estimate habitat affinities for each model group from empirical MaxN heatmaps using
# a continuous use-vs-availability SELECTION RATIO (Manly-style), over the sum-to-1 habitat
# basemaps written by R/habitat_basemaps.R: 6 GFISHER reef classes (AL,AM,AH,NL,NM,NH) plus
# rock, unconsolidated bottom and seagrass (RCK,UNC,SGR). See BASEMAP.SPEC below.
#
# This is the RASTER route, and it is one of three. It compares MaxN-weighted habitat use
# against availability averaged over model cells, so it inherits their scale: at 5 min a cell
# is ~9 km across while the camera saw one ~200 m spot. R/site_level_affinities.R estimates
# the same reef affinities at the station instead, and R/substrate_affinities.R handles the
# rock/gravel/sand/mud split from raw dbSeabed. Where they disagree, prefer the one whose
# support scale matches the question.
#
# Method (per group, restricted to SURVEYED cells, effort>0):
#   avail_h = mean H_h over surveyed cells                 (availability)
#   used_h  = sum(MaxN * H_h) / sum(MaxN)                  (MaxN-weighted use)
#   w_h     = used_h / avail_h                             (selection ratio)
#             w>1 selected FOR, w<1 selected AGAINST, w==1 used in proportion to availability.
#   A_h     = w_h / max(w_h WITHIN FAMILY)                 (affinity in [0,1], best in family = 1)
#
# WHY within-family: w is bounded above by 1/avail_h, and the two families sit on very different
# availability scales (reef layers average 1e-4..3e-2 of a cell, sediment 3e-2..6e-1). The four
# dbSeabed layers are also a CLOSED COMPOSITION -- they sum to ~1 per cell, which pins their
# availability-weighted mean w at ~1, so no group can select FOR all sediments. Normalizing over
# all 10 layers at once therefore handed A=1 to a rare reef layer in every group and rescaled the
# sediment affinities by an unrelated reef ratio. Normalizing inside each family removes that
# artifact. The cost: reef A and sediment A are each on their own scale, so A=1 reads as "best
# reef layer for this group" / "best sediment layer for this group", NOT "reef beats sediment".
#
# WHY surveyed-only: most water cells were never surveyed but are stored as 0 in the heatmaps;
# including them would compare use against availability the cameras never visited.
#
# IDENTIFIABILITY: a selection ratio is only meaningful where the survey actually sampled the
# habitat's gradient. The GFISHER video survey is reef-targeted, so it spans gravel/sand well,
# rock only at low-moderate cover, and essentially never sampled mud. fn.availability_coverage()
# computes a per-layer coverage flag (good/partial/poor) so under-sampled layers (mud) self-flag
# instead of fabricating "avoidance". Bootstrap CIs on w give the same warning statistically:
# poorly-covered layers come back with wide, 1-spanning intervals.

suppressPackageStartupMessages(library('raster'))

#--- layer spec: code -> filename pattern -> family ---------------------------------------
# Row order sets the column order of the wide affinity table and the bar order in the PDF:
# low -> medium -> high relief within the reef family, then bottom cover.
# `family` also defines the normalization groups for A (see header).
#
# BASEMAP.SPEC (default) reads the sum-to-1 layers written by R/habitat_basemaps.R. Those
# nine layers partition the cell, so all of them sit in ONE closed composition rather than
# the old open reef set plus a separate closed sediment set. Within-family normalization
# still matters, because reef layers average 3e-5..1.3e-2 of a cell against UNC's 0.86 and
# a single global max would hand A=1 to a rare reef layer every time.
#
# SEAGRASS IS DELIBERATELY EXCLUDED. The basemap still carries an SGR layer -- it has to, for
# the layers to sum to 1 -- but it is not identifiable from this survey: only 7.5% of the
# top-decile seagrass cells were ever surveyed and the sampled range reaches 67% of the layer
# maximum, so fn.availability_coverage flags it 'poor'. Included, it produced A=1.000 for six
# groups and A=0.000 for three red grouper stanzas off a handful of cells. Add 'SGR' back only
# if the survey footprint changes.
#
# Detailed substrate (rock / gravel / sand / mud) is NOT obtained here -- the basemap carries
# unconsolidated bottom as a single UNC layer. That split comes from R/substrate_affinities.R,
# which reads raw dbSeabed at its native 1.2 arc-min where the contrast actually exists.
#
# CAVEAT on the bottom family: with only RCK and UNC in it, and the two together covering ~97%
# of a cell, they are close to a closed two-part composition. w_RCK and w_UNC are therefore
# strongly anti-correlated and one of them is A=1 almost by construction. Read the bottom
# family as "rock vs unconsolidated preference", not as two independent affinities.
BASEMAP.SPEC <- data.frame(
  code    = c('AL','AM','AH','NL','NM','NH','RCK','UNC'),
  pattern = paste0('^habitat_', c('AL','AM','AH','NL','NM','NH','RCK','UNC'),
                   '_.*\\.asc$'),
  family  = c(rep('reef',6), rep('bottom',2)),
  stringsAsFactors = FALSE)

# LEGACY.SPEC reads the retired input_ascii_sum1/ layers -- the cell-area reef proportions
# and the processed gmf_*_val sediment layers. Kept so an older run can be reproduced, but
# see R/habitat_basemaps.R and R/substrate_affinities.R for why both are superseded: the reef
# proportions divide by cell area rather than scanned area, and the sediment layers
# renormalize rock against the grain-size triangle and treat NODATA as zero.
LEGACY.SPEC <- data.frame(
  code    = c('AL','AM','AH','NL','NM','NH','RCK','GVL','SND','MUD'),
  pattern = c('AL_prop.*\\.asc$','AM_prop.*\\.asc$','AH_prop.*\\.asc$',
              'NL_prop.*\\.asc$','NM_prop.*\\.asc$','NH_prop.*\\.asc$',
              'gmf_RCK_val.*\\.asc$','gmf_GVL_val.*\\.asc$',
              'gmf_SND_val.*\\.asc$','gmf_MUD_val.*\\.asc$'),
  family  = c(rep('reef',6), rep('sediment',4)),
  stringsAsFactors = FALSE)

LAYER.SPEC <- BASEMAP.SPEC   # default for every fn.* below

#--- load the habitat layers named in `spec`, labelled by code, as a RasterStack ----------
fn.load_layer_stack <- function(dir.hab, spec=LAYER.SPEC){
  out <- stack()
  for(i in seq_len(nrow(spec))){
    f <- list.files(dir.hab, pattern=spec$pattern[i], full.names=TRUE)
    f <- f[!grepl('\\.aux\\.xml$', f)]
    if(length(f)!=1) stop(sprintf("Expected exactly one '%s' in %s (found %d)",
                                  spec$pattern[i], dir.hab, length(f)))
    out <- addLayer(out, raster(f[1]))
  }
  names(out) <- spec$code
  out
}

#--- survey-effort raster: count of unique video stations per grid cell -------------------
# (Adapted from the retired estimate_habitat_affinities.R, removed in Oct 2026; this module runs standalone.)
fn.build_effort_raster <- function(file.env, depth, lon.col='lon_dd', lat.col='lat_dd',
                                   id.col='reference', save.as=NULL){
  e <- read.csv(file.env, header=TRUE, stringsAsFactors=FALSE)
  names(e) <- tolower(names(e))
  need <- c(id.col, lon.col, lat.col)
  if(any(!need %in% names(e))) stop('env file missing column(s): ',
                                    paste(need[!need %in% names(e)], collapse=', '))
  e <- e[!is.na(e[[lon.col]]) & !is.na(e[[lat.col]]), need]
  e <- unique(e)                                  # one point per station event
  sp::coordinates(e) <- as.formula(paste0('~',lon.col,'+',lat.col))
  sp::proj4string(e) <- sp::CRS('+proj=longlat +datum=WGS84 +no_defs')
  if(inherits(depth,'SpatRaster')) depth <- raster::raster(depth)
  e <- sp::spTransform(e, crs(depth))
  eff <- raster::rasterize(e, depth, field=1, fun='count', background=0)
  eff[is.na(depth)] <- NA
  names(eff) <- 'effort'
  if(!is.null(save.as)) writeRaster(eff, save.as, overwrite=TRUE)
  eff
}

#--- per-layer availability coverage (group-independent: depends only on habitat + effort)--
# surveyed_top_frac : fraction of the highest-availability cells (>= q quantile, & >0) that
#                     were surveyed -- the key identifiability metric.
# sampled_range_frac: max availability among surveyed cells / max availability overall.
# coverage flag     : good (top_frac>=.30) / partial (>=.10) / poor (<.10) / none (no habitat).
fn.availability_coverage <- function(hab.stack, effort, q=0.90){
  H <- getValues(hab.stack); ev <- getValues(effort)
  ok <- stats::complete.cases(H) & !is.na(ev)
  H <- H[ok,,drop=FALSE]; surv <- ev[ok] > 0
  do.call(rbind, lapply(colnames(H), function(c){
    x <- H[,c]
    hi <- x >= stats::quantile(x, q) & x > 0
    p.top <- if(any(hi)) mean(surv[hi]) else NA_real_
    rng   <- if(max(x) > 0) max(x[surv]) / max(x) else NA_real_
    flag  <- if(is.na(p.top)) 'none' else if(p.top>=0.30) 'good' else
             if(p.top>=0.10) 'partial' else 'poor'
    data.frame(code=c, surveyed_top_frac=round(p.top,3),
               sampled_range_frac=round(rng,3), coverage=flag, stringsAsFactors=FALSE)
  }))
}

#--- apply prior habitat constraints to a named selection-ratio vector --------------------
# w      : named numeric selection ratios (names = layer codes, e.g. RCK/GVL/SND/MUD/...).
# con    : list controlling the constraints (all optional):
#   $zero.codes  layers forced to 0 (default 'MUD' -- not identifiable from the reef survey).
#   $sand.code   layer that must sit strictly below every $hard.codes layer (default 'SND').
#   $hard.codes  hard-substrate layers sand must not exceed (default c('RCK','GVL')).
#   $sand.margin fractional gap kept below the smaller hard layer (default 0.05).
# Only lowers sand (rock/gravel keep their estimated values); groups already satisfying the
# ordering are untouched. Returns the constrained w; affinity A is then recomputed from it.
fn.constrain_w <- function(w, con=list()){
  zero.codes  <- if(!is.null(con$zero.codes))  con$zero.codes  else 'MUD'
  sand.code   <- if(!is.null(con$sand.code))   con$sand.code   else 'SND'
  hard.codes  <- if(!is.null(con$hard.codes))  con$hard.codes  else c('RCK','GVL')
  sand.margin <- if(!is.null(con$sand.margin)) con$sand.margin else 0.05
  w[intersect(zero.codes, names(w))] <- 0
  if(!is.na(sand.code) && sand.code %in% names(w) && !is.na(w[[sand.code]])){
    hard <- w[intersect(hard.codes, names(w))]; hard <- hard[!is.na(hard)]
    if(length(hard) > 0){
      cap <- min(hard) * (1 - sand.margin)
      if(w[[sand.code]] > cap) w[[sand.code]] <- cap
    }
  }
  w
}

#--- affinity from constrained selection ratios, normalized WITHIN each family ------------
# w   : named numeric selection ratios (already constrained).
# fam : family label per element of w ('reef' / 'sediment' / ...), same length as w.
# Each family is rescaled by its own maximum, so the best reef layer and the best sediment
# layer both come back as 1. See the header for why a single global max is not used.
# A family whose ratios are all NA, all zero, or non-finite comes back as NA.
# normalize = 'global' (default) rescales every layer by ONE maximum. 'family' rescales each
# family by its own.
#
# 'family' was the right call for the retired input_ascii_sum1 layers, where a global max
# handed A=1 to a rare reef layer in all 18 fitted groups and never once to a sediment layer.
# It is the WRONG call for the sum-to-1 basemaps. Those layers are one partition of the cell
# and their realized ratios are all modest -- AL 0.15-2.14, AH 0.82-2.92, NH 0.77-2.78,
# RCK 0.62-1.68, UNC 0.93-1.02 -- so the rarity artifact stays theoretical: under a global max
# the winning layer is spread across AH/AL/AM/NH/NL, never RCK or UNC.
#
# Normalizing per family instead forced the bottom pair up to A=1 for most groups, telling
# Ecospace that reef fish prefer unconsolidated bottom as strongly as they prefer reef. Under
# a global max the reef columns are UNCHANGED (the maximum always falls in the reef family
# anyway) and mean bottom affinity drops from 0.91/0.89 to 0.56/0.53, which is the intended
# behaviour for reef-associated groups.
fn.affinity_from_w <- function(w, fam, normalize='global'){
  A <- rep(NA_real_, length(w))
  grp <- if(identical(normalize, 'family')) fam else rep('all', length(w))
  for(fm in unique(grp)){
    j  <- which(grp == fm)
    mx <- suppressWarnings(max(w[j], na.rm=TRUE))
    if(is.finite(mx) && mx > 0) A[j] <- w[j] / mx
  }
  names(A) <- names(w)
  A
}

#--- selection ratios for ONE empirical MaxN layer ----------------------------------------
# Returns a per-layer data.frame: code, family, avail, used, w, w_lo, w_hi, w_con, A, sig.
# When con$apply is TRUE the constraints above are applied to w -> w_con and A is derived from
# w_con; the raw w / w_lo / w_hi / sig are left untouched as the empirical evidence.
# n.boot bootstrap resamples of surveyed cells give a 95% CI on w; sig = 'for'/'against'/'ns'
# from whether that CI clears 1. min.pos guards groups with too few non-zero MaxN cells.
fn.selection_ratios <- function(hab.stack, emp.ras, effort, n.boot=1000, seed=1,
                                min.pos=5, spec=LAYER.SPEC, con=NULL, normalize='global'){
  if(!compareRaster(hab.stack, emp.ras, extent=TRUE, rowcol=TRUE, crs=FALSE, stopiffalse=FALSE))
    stop("habitat stack and empirical raster do not share the same grid")
  H <- getValues(hab.stack); E <- getValues(emp.ras); ev <- getValues(effort)
  ok <- stats::complete.cases(H) & !is.na(E) & !is.na(ev) & ev > 0
  H <- H[ok,,drop=FALSE]; E <- E[ok]
  codes <- colnames(H); K <- ncol(H); n <- length(E)
  fam <- spec$family[match(codes, spec$code)]

  na.df <- function(){
    data.frame(code=codes, family=fam, avail=NA_real_, used=NA_real_, w=NA_real_,
               w_lo=NA_real_, w_hi=NA_real_, w_con=NA_real_, A=NA_real_, sig='na',
               n=n, n_pos=sum(E>0), stringsAsFactors=FALSE)
  }
  if(n < 10 || sum(E) <= 0 || sum(E>0) < min.pos) return(na.df())

  avail <- colMeans(H)
  used  <- colSums(E * H) / sum(E)
  w <- used / avail; w[avail==0] <- NA_real_

  set.seed(seed)
  B <- matrix(NA_real_, n.boot, K)
  for(b in seq_len(n.boot)){
    idx <- sample.int(n, n, replace=TRUE)
    sb  <- sum(E[idx]); if(sb <= 0) next
    ab  <- colMeans(H[idx,,drop=FALSE])
    ub  <- colSums(E[idx] * H[idx,,drop=FALSE]) / sb
    r   <- ub/ab; r[ab==0] <- NA_real_
    B[b,] <- r
  }
  w.lo <- apply(B, 2, stats::quantile, 0.025, na.rm=TRUE)
  w.hi <- apply(B, 2, stats::quantile, 0.975, na.rm=TRUE)
  # constrain the selection ratios (MUD=0, sand < rock/gravel) before forming affinities
  names(w) <- codes
  w.con <- if(!is.null(con) && isTRUE(con$apply)) fn.constrain_w(w, con) else w
  A <- fn.affinity_from_w(w.con, fam, normalize=normalize)
  sig  <- ifelse(is.na(w), 'na', ifelse(w.lo > 1, 'for', ifelse(w.hi < 1, 'against', 'ns')))

  data.frame(code=codes, family=fam,
             avail=round(avail,5), used=round(used,5), w=round(w,3),
             w_lo=round(w.lo,3), w_hi=round(w.hi,3), w_con=round(w.con,3), A=round(A,3), sig=sig,
             n=n, n_pos=sum(E>0), stringsAsFactors=FALSE)
}

#--- split a modname into species base + stanza index -------------------------------------
# Multi-stanza groups are named '<species>-<k>', with a trailing '-' on the plus group
# ('gag-0' ... 'gag-5-'). Single-stanza groups ('sharks') get stanza NA and are never pooled.
fn.parse_stanza <- function(modname){
  has <- grepl('-\\d+-?$', modname)
  data.frame(base   = ifelse(has, sub('-\\d+-?$', '', modname), modname),
             stanza = suppressWarnings(as.integer(
                        ifelse(has, sub('.*-(\\d+)-?$', '\\1', modname), NA))),
             stringsAsFactors = FALSE)
}

#--- pool sparse stanzas with their nearest same-species neighbours -----------------------
# A stanza observed in fewer than min.pos surveyed cells cannot support a selection ratio --
# fn.selection_ratios returns an all-NA row for it. Rather than leave a hole to be patched by
# hand downstream, merge its MaxN with the nearest stanza(s) of the SAME species (adding one
# neighbour at a time, nearest stanza first) until the pooled layer clears min.pos, and fit
# that instead. Neighbours keep their own independent fits -- only the sparse stanza's row
# changes -- and `pooled_with` records what went into it so pooled rows stay identifiable.
# info    : list of per-group records (see fn.batch_selection_ratios); each has $ras/$npos.
# npos.fn : counts surveyed cells with MaxN>0 in a raster.
fn.pool_sparse_stanzas <- function(info, npos.fn, min.pos=5){
  nm  <- vapply(info, `[[`, character(1), 'modname')
  stz <- fn.parse_stanza(nm)
  sparse <- which(vapply(info, `[[`, integer(1), 'npos') < min.pos & !is.na(stz$stanza))
  for(i in sparse){
    cand <- setdiff(which(stz$base == stz$base[i] & !is.na(stz$stanza)), i)
    if(length(cand) == 0) next
    cand <- cand[order(abs(stz$stanza[cand] - stz$stanza[i]), stz$stanza[cand])]
    ras <- info[[i]]$ras; part <- character(0); np <- info[[i]]$npos
    for(j in cand){
      ras  <- ras + info[[j]]$ras          # same grid & mask, so NA cells stay NA
      part <- c(part, nm[j])
      np   <- npos.fn(ras)
      if(np >= min.pos) break
    }
    message(sprintf("    mod%d %s: only %d surveyed cell(s) with MaxN>0 -- pooled with %s (now %d)",
                    info[[i]]$modnumber, nm[i], info[[i]]$npos, paste(part, collapse='+'), np))
    info[[i]]$ras <- ras
    info[[i]]$pooled_with <- paste(part, collapse='+')
    info[[i]]$npos <- np
  }
  info
}

#--- batch over every empirical MaxN layer in dir.emp -------------------------------------
# min.pos      : minimum surveyed cells with MaxN>0 required to fit a group.
# pool.stanzas : TRUE (default) pools sparse stanzas with same-species neighbours first;
#                FALSE reproduces the old behaviour of returning an all-NA row for them.
fn.batch_selection_ratios <- function(hab.stack, dir.emp, dir.out, effort,
                                      n.boot=1000, seed=1, spec=LAYER.SPEC, con=NULL, normalize='global',
                                      min.pos=5, pool.stanzas=TRUE){
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
  files <- list.files(dir.emp, pattern='\\.asc$', full.names=TRUE)
  files <- files[grepl('GFISHER_maxn_mod', basename(files))]
  if(length(files)==0) stop('no GFISHER_maxn_mod*.asc files found in ', dir.emp)
  files <- files[order(as.integer(sub('.*_mod(\\d+)_.*','\\1', basename(files))))]

  # resolution tag read off the habitat grid rather than hardcoded, so a 15min run labels its
  # outputs (and parses its modnames) correctly.
  tag <- paste0(round(raster::res(hab.stack)[1]*60, 0), 'min')

  cov <- fn.availability_coverage(hab.stack, effort)   # per-layer, group-independent

  # the cells fn.selection_ratios will actually use; MaxN>0 counted over these drives pooling
  ev  <- getValues(effort)
  ok0 <- stats::complete.cases(getValues(hab.stack)) & !is.na(ev) & ev > 0
  npos.fn <- function(r){ v <- getValues(r)[ok0]; as.integer(sum(!is.na(v) & v > 0)) }

  info <- lapply(files, function(f){
    bn <- basename(f); r <- raster(f)
    list(file=f, bn=bn,
         modnumber   = as.integer(sub('.*_mod(\\d+)_.*', '\\1', bn)),
         modname     = sub('.*_mod\\d+_(.+)_\\d+min_.*', '\\1', bn),
         ras         = r,
         npos        = npos.fn(r),
         pooled_with = NA_character_)
  })
  if(isTRUE(pool.stanzas)) info <- fn.pool_sparse_stanzas(info, npos.fn, min.pos=min.pos)

  long <- list(); plots <- list()
  for(p in info){
    cat(sprintf("  mod%-3d %s ...\n", p$modnumber, p$modname))
    sr  <- fn.selection_ratios(hab.stack, p$ras, effort, n.boot=n.boot, seed=seed,
                               min.pos=min.pos, spec=spec, con=con, normalize=normalize)
    sr$coverage  <- cov$coverage[match(sr$code, cov$code)]
    sr$modnumber <- p$modnumber; sr$modname <- p$modname
    sr$pooled_with <- p$pooled_with
    long[[p$bn]] <- sr
    plots[[p$bn]] <- list(modname=p$modname, modnumber=p$modnumber, sr=sr,
                          pooled_with=p$pooled_with)
  }
  long <- do.call(rbind, long)
  long <- long[, c('modnumber','modname','code','family','avail','used',
                   'w','w_lo','w_hi','w_con','A','sig','coverage','n','n_pos','pooled_with')]

  # wide affinity table (rows = groups, cols = layers, value = A) for Ecospace input
  ord  <- spec$code
  wide <- reshape(long[,c('modnumber','modname','code','A')],
                  idvar=c('modnumber','modname'), timevar='code', direction='wide')
  names(wide) <- sub('^A\\.','', names(wide))
  wide <- wide[, c('modnumber','modname', ord[ord %in% names(wide)])]
  wide <- wide[order(wide$modnumber),]

  f.long <- file.path(dir.out, paste0('selection_ratios_long_', tag, '.csv'))
  f.wide <- file.path(dir.out, paste0('affinity_A_wide_', tag, '.csv'))
  f.cov  <- file.path(dir.out, paste0('availability_coverage_', tag, '.csv'))
  f.pdf  <- file.path(dir.out, paste0('selection_ratio_fits_', tag, '.pdf'))
  write.csv(long, f.long, row.names=FALSE)
  write.csv(wide, f.wide, row.names=FALSE, quote=FALSE)  # modnames are already [A-Za-z0-9-] only
  write.csv(cov,  f.cov,  row.names=FALSE)

  # multipage PDF: per group, barplot of selection ratio w with bootstrap CI, line at w=1,
  # bars colored by coverage flag.
  cov.col <- c(good='grey35', partial='darkorange', poor='red3', none='grey80', na='grey80')
  pdf(f.pdf, width=10, height=7.5, onefile=TRUE)
  op <- par(mfrow=c(2,2), mar=c(4,4,3,1))
  for(p in plots){
    s <- p$sr; s <- s[match(ord, s$code),]
    if(all(is.na(s$w))){ plot.new(); title(paste0(p$modname,'\n(no fit)')); next }
    yhi <- max(s$w_hi, s$w, 1.05, na.rm=TRUE)
    sub <- if(!is.na(p$pooled_with)) paste0('  [pooled with ', p$pooled_with, ']') else ''
    bp <- barplot(s$w, names.arg=s$code, las=2, ylim=c(0, yhi),
                  col=cov.col[s$coverage], border=NA,
                  ylab='selection ratio  (used / available)',
                  main=sprintf('mod%d  %s%s\n(n=%d surveyed, %d with MaxN>0)',
                               p$modnumber, p$modname, sub, s$n[1], s$n_pos[1]))
    # family divider, for reading only -- with normalize='global' A is scaled by one maximum
    fb <- which(diff(as.integer(factor(s$family, levels=unique(s$family)))) != 0)
    if(length(fb)) abline(v=(bp[fb]+bp[fb+1])/2, lty=3, col='grey60')
    suppressWarnings(arrows(bp, s$w_lo, bp, s$w_hi, angle=90, code=3, length=0.03, col='black'))
    abline(h=1, lty=2, col='blue')
    # red dash = constrained selection ratio used for affinity (MUD=0, sand < rock/gravel)
    chg <- which(!is.na(s$w_con) & (is.na(s$w) | abs(s$w_con - s$w) > 1e-9))
    if(length(chg)) points(bp[chg], s$w_con[chg], pch=45, col='red', cex=2.4, lwd=2)
    sigp <- s$sig %in% c('for','against')
    if(any(sigp)) text(bp[sigp], s$w_hi[sigp], '*', pos=3, offset=0.2, col='black', cex=1.3)
  }
  par(op)
  plot.new(); legend('center', title='availability coverage', bty='n',
                     fill=cov.col[c('good','partial','poor')],
                     legend=c('good','partial','poor (e.g. mud — not identifiable)'))
  legend('bottom', bty='n', pch=c(45,NA), lty=c(NA,3), col=c('red','grey60'), pt.cex=2.4,
         legend=c('constrained ratio used for affinity (MUD=0, sand < rock & gravel)',
                  'family divider (reef | bottom); A is scaled by one global maximum'))
  dev.off()

  cat('\nWrote:\n  ', f.long, '\n  ', f.wide, '\n  ', f.cov, '\n  ', f.pdf, '\n', sep='')
  invisible(list(long=long, wide=wide, coverage=cov))
}

#================================================================================
# DRIVER -- mode 'batch' (default, all groups) or 'pilot' (one group)
#================================================================================
if(sys.nframe()==0){
  args <- commandArgs(trailingOnly=TRUE)
  mode    <- if(length(args)>=1) args[1] else 'batch'
  dir.hab <- if(length(args)>=2) args[2] else
    "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/input_ascii_sum1/5min"
  dir.emp <- if(length(args)>=3) args[3] else
    "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/GFISHER/5min/maxn"
  dir.out <- if(length(args)>=4) args[4] else file.path(dir.emp,'affinity_selratio')
  file.env <- if(length(args)>=5) args[5] else
    "C:/Users/dchagaris/Github/WFS-FEM/GFISHER/data/April2026/env3LABS_93to24.csv"
  n.boot   <- if(length(args)>=6) as.integer(args[6]) else 1000
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)

  cat("Loading 10-layer habitat stack...\n")
  hab <- fn.load_layer_stack(dir.hab)
  cat("Building survey-effort raster from", basename(file.env), "...\n")
  eff <- fn.build_effort_raster(file.env, hab[[1]],
           save.as=file.path(dir.out, sprintf('GFISHER_survey_effort_%dmin_%dx%d.asc',
                                              round(raster::res(hab)[1]*60,0), nrow(hab), ncol(hab))))
  cat("  surveyed cells:", sum(getValues(eff)>0, na.rm=TRUE),
      " total stations:", sum(getValues(eff), na.rm=TRUE), "\n")

  cat("\nPer-layer availability coverage:\n")
  print(fn.availability_coverage(hab, eff), row.names=FALSE)

  if(mode=='batch'){
    cat("\nEstimating selection ratios for ALL groups (n.boot=", n.boot, ")...\n", sep='')
    res <- fn.batch_selection_ratios(hab, dir.emp, dir.out, effort=eff, n.boot=n.boot)
    cat("\n--- AFFINITY A WIDE (rows=groups, cols=layers) ---\n")
    print(res$wide, row.names=FALSE)

  } else if(mode=='pilot'){
    pilot <- if(length(args)>=7) args[7] else 'mod44_hogfish'
    f.emp <- list.files(dir.emp, pattern=paste0(pilot,'.*\\.asc$'), full.names=TRUE)
    if(length(f.emp)!=1) stop(paste("pilot empirical layer not uniquely found:", pilot))
    emp <- raster(f.emp[1])
    sr  <- fn.selection_ratios(hab, emp, eff, n.boot=n.boot)
    sr$coverage <- fn.availability_coverage(hab, eff)$coverage[match(sr$code, LAYER.SPEC$code)]
    cat("\n--- PILOT:", pilot, "---\n"); print(sr, row.names=FALSE)
  } else stop("unknown mode: ", mode)
}
