# selection_ratio_affinities.R --------------------------------------------------
# Estimate habitat affinities for each model group from empirical MaxN heatmaps using
# a continuous use-vs-availability SELECTION RATIO (Manly-style), over the 10 habitat
# layers: 6 GFISHER reef classes (AH,AM,AL,NH,NM,NL) + 4 dbSeabed sediment classes
# (RCK,GVL,SND,MUD).
#
# Method (per group, restricted to SURVEYED cells, effort>0):
#   avail_h = mean H_h over surveyed cells                 (availability)
#   used_h  = sum(MaxN * H_h) / sum(MaxN)                  (MaxN-weighted use)
#   w_h     = used_h / avail_h                             (selection ratio)
#             w>1 selected FOR, w<1 selected AGAINST, w==1 used in proportion to availability.
#   A_h     = w_h / max(w_h)                               (affinity in [0,1], best habitat = 1)
# Every layer's ratio is on the same scale ("x more than its average availability"), so reef
# and sediment affinities are directly comparable within a group.
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

#--- layer spec: code -> filename pattern -> family. Default = the 10 layers above. -------
# Extend this data.frame (e.g. add seagrass) to bring more layers into the analysis.
LAYER.SPEC <- data.frame(
  code    = c('AH','AM','AL','NH','NM','NL','RCK','GVL','SND','MUD'),
  pattern = c('AH_prop.*\\.asc$','AM_prop.*\\.asc$','AL_prop.*\\.asc$',
              'NH_prop.*\\.asc$','NM_prop.*\\.asc$','NL_prop.*\\.asc$',
              'gmf_RCK_val.*\\.asc$','gmf_GVL_val.*\\.asc$',
              'gmf_SND_val.*\\.asc$','gmf_MUD_val.*\\.asc$'),
  family  = c(rep('reef',6), rep('sediment',4)),
  stringsAsFactors = FALSE)

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
# (Adapted from R/estimate_habitat_affinities.R so this module runs standalone.)
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

#--- selection ratios for ONE empirical MaxN layer ----------------------------------------
# Returns a per-layer data.frame: code, family, avail, used, w, w_lo, w_hi, A, sig.
# n.boot bootstrap resamples of surveyed cells give a 95% CI on w; sig = 'for'/'against'/'ns'
# from whether that CI clears 1. min.pos guards groups with too few non-zero MaxN cells.
fn.selection_ratios <- function(hab.stack, emp.ras, effort, n.boot=1000, seed=1,
                                min.pos=5, spec=LAYER.SPEC){
  if(!compareRaster(hab.stack, emp.ras, extent=TRUE, rowcol=TRUE, crs=FALSE, stopiffalse=FALSE))
    stop("habitat stack and empirical raster do not share the same grid")
  H <- getValues(hab.stack); E <- getValues(emp.ras); ev <- getValues(effort)
  ok <- stats::complete.cases(H) & !is.na(E) & !is.na(ev) & ev > 0
  H <- H[ok,,drop=FALSE]; E <- E[ok]
  codes <- colnames(H); K <- ncol(H); n <- length(E)
  fam <- spec$family[match(codes, spec$code)]

  na.df <- function(){
    data.frame(code=codes, family=fam, avail=NA_real_, used=NA_real_, w=NA_real_,
               w_lo=NA_real_, w_hi=NA_real_, A=NA_real_, sig='na',
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
  A    <- w / max(w, na.rm=TRUE)
  sig  <- ifelse(is.na(w), 'na', ifelse(w.lo > 1, 'for', ifelse(w.hi < 1, 'against', 'ns')))

  data.frame(code=codes, family=fam,
             avail=round(avail,5), used=round(used,5), w=round(w,3),
             w_lo=round(w.lo,3), w_hi=round(w.hi,3), A=round(A,3), sig=sig,
             n=n, n_pos=sum(E>0), stringsAsFactors=FALSE)
}

#--- batch over every empirical MaxN layer in dir.emp -------------------------------------
fn.batch_selection_ratios <- function(hab.stack, dir.emp, dir.out, effort,
                                      n.boot=1000, seed=1, spec=LAYER.SPEC){
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
  files <- list.files(dir.emp, pattern='\\.asc$', full.names=TRUE)
  files <- files[grepl('GFISHER_maxn_mod', basename(files))]
  if(length(files)==0) stop('no GFISHER_maxn_mod*.asc files found in ', dir.emp)
  files <- files[order(as.integer(sub('.*_mod(\\d+)_.*','\\1', basename(files))))]

  cov <- fn.availability_coverage(hab.stack, effort)   # per-layer, group-independent

  long <- list(); plots <- list()
  for(f in files){
    bn <- basename(f)
    modnumber <- as.integer(sub('.*_mod(\\d+)_.*', '\\1', bn))
    modname   <- sub('.*_mod\\d+_(.+)_5min_.*', '\\1', bn)
    cat(sprintf("  mod%-3d %s ...\n", modnumber, modname))
    emp <- raster(f)
    sr  <- fn.selection_ratios(hab.stack, emp, effort, n.boot=n.boot, seed=seed, spec=spec)
    sr$coverage  <- cov$coverage[match(sr$code, cov$code)]
    sr$modnumber <- modnumber; sr$modname <- modname
    long[[bn]] <- sr
    plots[[bn]] <- list(modname=modname, modnumber=modnumber, sr=sr)
  }
  long <- do.call(rbind, long)
  long <- long[, c('modnumber','modname','code','family','avail','used',
                   'w','w_lo','w_hi','A','sig','coverage','n','n_pos')]

  # wide affinity table (rows = groups, cols = layers, value = A) for Ecospace input
  ord  <- spec$code
  wide <- reshape(long[,c('modnumber','modname','code','A')],
                  idvar=c('modnumber','modname'), timevar='code', direction='wide')
  names(wide) <- sub('^A\\.','', names(wide))
  wide <- wide[, c('modnumber','modname', ord[ord %in% names(wide)])]
  wide <- wide[order(wide$modnumber),]

  f.long <- file.path(dir.out, 'selection_ratios_long_5min.csv')
  f.wide <- file.path(dir.out, 'affinity_A_wide_5min.csv')
  f.cov  <- file.path(dir.out, 'availability_coverage_5min.csv')
  write.csv(long, f.long, row.names=FALSE)
  write.csv(wide, f.wide, row.names=FALSE)
  write.csv(cov,  f.cov,  row.names=FALSE)

  # multipage PDF: per group, barplot of selection ratio w with bootstrap CI, line at w=1,
  # bars colored by coverage flag.
  cov.col <- c(good='grey35', partial='darkorange', poor='red3', none='grey80', na='grey80')
  pdf(file.path(dir.out, 'selection_ratio_fits_5min.pdf'), width=10, height=7.5, onefile=TRUE)
  op <- par(mfrow=c(2,2), mar=c(4,4,3,1))
  for(p in plots){
    s <- p$sr; s <- s[match(ord, s$code),]
    if(all(is.na(s$w))){ plot.new(); title(paste0(p$modname,'\n(no fit)')); next }
    yhi <- max(s$w_hi, s$w, 1.05, na.rm=TRUE)
    bp <- barplot(s$w, names.arg=s$code, las=2, ylim=c(0, yhi),
                  col=cov.col[s$coverage], border=NA,
                  ylab='selection ratio  (used / available)',
                  main=sprintf('mod%d  %s\n(n=%d surveyed, %d with MaxN>0)',
                               p$modnumber, p$modname, s$n[1], s$n_pos[1]))
    suppressWarnings(arrows(bp, s$w_lo, bp, s$w_hi, angle=90, code=3, length=0.03, col='black'))
    abline(h=1, lty=2, col='blue')
    sigp <- s$sig %in% c('for','against')
    if(any(sigp)) text(bp[sigp], s$w_hi[sigp], '*', pos=3, offset=0.2, col='black', cex=1.3)
  }
  par(op)
  plot.new(); legend('center', title='availability coverage', bty='n',
                     fill=cov.col[c('good','partial','poor')],
                     legend=c('good','partial','poor (e.g. mud — not identifiable)'))
  dev.off()

  cat('\nWrote:\n  ', f.long, '\n  ', f.wide, '\n  ', f.cov,
      '\n  ', file.path(dir.out,'selection_ratio_fits_5min.pdf'), '\n', sep='')
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
           save.as=file.path(dir.out,'GFISHER_survey_effort_5min_66x78.asc'))
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
