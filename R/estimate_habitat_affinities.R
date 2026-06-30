# estimate_habitat_affinities.R --------------------------------------------------
# Estimate a 6-element habitat-affinity vector A (AH,AM,AL,NH,NM,NL), each in [0,1],
# for each model group with an empirical MaxN heatmap.
#
# Model (no intercept):  P(cell) = sum_h a_h * H_h(cell)
# Objective: choose A to MAXIMIZE the Spearman (rank) correlation between the predicted
#            heatmap P and the empirical MaxN heatmap E over the common valid cells.
#
# Spearman corr is invariant to any positive monotonic rescaling of P, so only the
# DIRECTION of A (relative affinities) is identified, not its magnitude. We optimize over
# the unit simplex (a_h >= 0, sum a_h = 1) via a softmax reparameterization, then rescale
# the reported A so max(a_h) = 1 ("best habitat = 1").

suppressPackageStartupMessages(library('raster'))

HAB.CODES <- c('AH','AM','AL','NH','NM','NL')

#--- load the 6 habitat-proportion layers, labelled by 2-letter code --------------
fn.load_habitat_stack <- function(dir.hab){
  out <- stack()
  for(code in HAB.CODES){
    f <- list.files(dir.hab, pattern=paste0(code,'_prop.*\\.asc$'), full.names=TRUE)
    if(length(f)!=1) stop(paste0("Expected exactly one '",code,"_prop*.asc' in ",dir.hab,
                                 " (found ",length(f),")"))
    out <- addLayer(out, raster(f[1]))
  }
  names(out) <- HAB.CODES
  out
}

#--- survey-effort raster: count of unique video stations per grid cell -----------
# Reads the env table (one row per station event), counts unique stations per cell on the
# depth/habitat grid. Used to control for sampling: most water cells were never surveyed but
# are stored as 0 in the empirical heatmaps, which spuriously rewards near-zero predictors
# (e.g. point-like artificial-high-relief reef) for matching those unsampled zeros.
fn.build_effort_raster <- function(file.env, depth, lon.col='lon_dd', lat.col='lat_dd',
                                   id.col='reference', save.as=NULL){
  e <- read.csv(file.env, header=TRUE, stringsAsFactors=FALSE)
  names(e) <- tolower(names(e))
  need <- c(id.col, lon.col, lat.col)
  if(any(!need %in% names(e))) stop('env file missing column(s): ',
                                    paste(need[!need %in% names(e)], collapse=', '))
  e <- e[!is.na(e[[lon.col]]) & !is.na(e[[lat.col]]), need]
  e <- unique(e)                                  # one point per station event
  coordinates(e) <- as.formula(paste0('~',lon.col,'+',lat.col))
  proj4string(e) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  if(inherits(depth,'SpatRaster')) depth <- raster::raster(depth)
  e <- spTransform(e, crs(depth))
  eff <- raster::rasterize(e, depth, field=1, fun='count', background=0)
  eff[is.na(depth)] <- NA
  names(eff) <- 'effort'
  if(!is.null(save.as)) writeRaster(eff, save.as, overwrite=TRUE)
  eff
}

#--- softmax: unconstrained x in R^6  ->  a on the unit simplex -------------------
.softmax <- function(x){ e <- exp(x - max(x)); e/sum(e) }

#--- estimate A for one empirical layer ------------------------------------------
# effort         : optional survey-effort raster (same grid) from fn.build_effort_raster().
# effort.control : 'none'     -> use all valid cells, plain Spearman (original behaviour);
#                  'restrict' -> keep only surveyed cells (effort>0);
#                  'partial'  -> all valid cells, partial rank corr controlling for effort;
#                  'both'     -> restrict to surveyed cells AND partial out effort (default).
fn.estimate_habitat_affinities <- function(hab.stack, emp.ras, effort=NULL,
                                           effort.control='both', n.random=4000,
                                           n.polish=12, seed=1){
  if(!compareRaster(hab.stack, emp.ras, extent=TRUE, rowcol=TRUE, crs=FALSE, stopiffalse=FALSE))
    stop("habitat stack and empirical raster do not share the same grid")
  H   <- getValues(hab.stack)               # n x 6
  E   <- getValues(emp.ras)                 # n
  eff <- if(!is.null(effort)) getValues(effort) else NULL
  use.eff <- !is.null(eff) && effort.control != 'none'

  # common valid mask (all habitat layers AND empirical non-NA); optionally surveyed-only
  ok <- stats::complete.cases(H) & !is.na(E)
  if(use.eff && effort.control %in% c('restrict','both')) ok <- ok & !is.na(eff) & eff > 0
  H <- H[ok,,drop=FALSE]; E <- E[ok]
  effv <- if(use.eff) eff[ok] else NULL
  n <- length(E)

  # guard: need variation in E and in the predictors
  if(n < 10 || stats::sd(E)==0 || all(colSums(H)==0))
    return(list(A=setNames(rep(NA_real_,6),HAB.CODES), A.simplex=NA, rho=NA_real_,
                spearman=NA_real_, pearson=NA_real_, n=n, control=effort.control,
                pred=NULL, mask=ok))

  # partial control: residualize the *ranks* of both E and P on the rank of effort, then
  # correlate the residuals (a partial Spearman). Pre-compute effort's rank-residual maker.
  partial <- use.eff && effort.control %in% c('partial','both')
  if(partial){
    reff  <- rank(effv); reffc <- reff - mean(reff); denom <- sum(reffc^2)
    if(denom==0) partial <- FALSE
  }
  resid.eff <- function(v){ vc <- v - mean(v); if(partial) vc - (sum(vc*reffc)/denom)*reffc else vc }

  Ey <- resid.eff(rank(E))                  # response (rank-residualized once)
  score <- function(a){                     # partial (or plain) rank correlation of P with E
    p <- as.numeric(H %*% a)
    if(stats::sd(p)==0) return(-Inf)
    py <- resid.eff(rank(p))
    if(stats::sd(py)==0) return(-Inf)
    suppressWarnings(stats::cor(py, Ey))
  }
  obj <- function(x) -score(.softmax(x))    # minimize -score

  # global pre-screen over the simplex (Dirichlet draws + structured corners/uniform),
  # then polish the best few with Nelder-Mead. The objective is piecewise-constant in A, so a
  # broad random screen guards against the plateaus a local optimizer would get stuck on.
  set.seed(seed)
  rand <- matrix(rgamma(n.random*6, shape=0.5), ncol=6); rand <- rand/rowSums(rand)
  corners <- diag(6)*0.95 + 0.05/6          # near each pure habitat
  starts <- rbind(rep(1/6,6), corners, rand)
  sc.s <- apply(starts, 1, score)
  top  <- order(sc.s, decreasing=TRUE)[seq_len(min(n.polish, nrow(starts)))]

  best <- list(val=Inf, par=NULL)
  for(i in top){
    a0 <- pmax(starts[i,], 1e-6); x0 <- log(a0)         # simplex point -> softmax space
    fit <- stats::optim(x0, obj, method='Nelder-Mead',
                        control=list(maxit=2000, reltol=1e-10))
    if(fit$value < best$val) best <- list(val=fit$value, par=fit$par)
  }

  a.simplex <- .softmax(best$par)
  A <- a.simplex / max(a.simplex)            # rescale so best habitat = 1
  names(A) <- names(a.simplex) <- HAB.CODES
  p <- as.numeric(H %*% a.simplex)
  list(A=round(A,4), A.simplex=round(a.simplex,4),
       rho      = round(-best$val, 4),                          # objective (partial if controlled)
       spearman = round(suppressWarnings(stats::cor(rank(p), rank(E))), 4), # plain, reference
       pearson  = round(suppressWarnings(stats::cor(p, E)), 4),            # reference only
       n = n, control = effort.control,
       pred = { r <- emp.ras; r[] <- NA; r[which(ok)] <- p; r },
       mask = ok)
}

#--- batch over every empirical MaxN layer in dir.emp -----------------------------
fn.batch_affinities <- function(hab.stack, dir.emp, dir.out, effort=NULL,
                                effort.control='both', write.asc=TRUE, seed=1){
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)
  files <- list.files(dir.emp, pattern='\\.asc$', full.names=TRUE)
  files <- files[grepl('GFISHER_maxn_mod', basename(files))]
  if(length(files)==0) stop('no GFISHER_maxn_mod*.asc files found in ', dir.emp)
  files <- files[order(as.integer(sub('.*_mod(\\d+)_.*','\\1', basename(files))))]

  rows <- list(); preds <- list()
  for(f in files){
    bn <- basename(f)
    modnumber <- as.integer(sub('.*_mod(\\d+)_.*', '\\1', bn))
    modname   <- sub('.*_mod\\d+_(.+)_5min_.*', '\\1', bn)
    cat(sprintf("  mod%-3d %s ...\n", modnumber, modname))
    emp <- raster(f)
    res <- fn.estimate_habitat_affinities(hab.stack, emp, effort=effort,
                                          effort.control=effort.control, seed=seed)
    rows[[bn]]  <- data.frame(modnumber=modnumber, modname=modname,
                              as.list(res$A), rho=res$rho, spearman=res$spearman,
                              pearson=res$pearson, n=res$n, control=res$control,
                              check.names=FALSE, stringsAsFactors=FALSE)
    preds[[bn]] <- list(modname=modname, rho=res$rho, emp=emp, pred=res$pred)
    if(write.asc && !is.null(res$pred)){
      fout <- file.path(dir.out, paste0('GFISHER_affinity_pred_mod',modnumber,'_',
                        gsub('[^A-Za-z0-9]+','-',modname),'_5min_',nrow(emp),'x',ncol(emp),'.asc'))
      writeRaster(res$pred, fout, overwrite=TRUE)
    }
  }
  out <- do.call(rbind, rows)
  fcsv <- file.path(dir.out,'habitat_affinities_5min.csv')
  write.csv(out, fcsv, row.names=FALSE)

  # multipage PDF: 3 groups/page, each row = empirical | predicted
  pdf(file.path(dir.out,'habitat_affinity_fits_5min.pdf'), width=9, height=11, onefile=TRUE)
  op <- par(mfrow=c(3,2), mar=c(3,3,3,5))
  for(p in preds){
    plot(p$emp, colNA='black', main=paste0(p$modname,'\nempirical MaxN'))
    if(!is.null(p$pred)) plot(p$pred, colNA='black', main=paste0('predicted  rho=',round(p$rho,3)))
    else { plot.new(); title('no fit (degenerate)') }
  }
  par(op); dev.off()
  cat('\nWrote ', fcsv, ' and habitat_affinity_fits_5min.pdf to ', dir.out, '\n', sep='')
  invisible(out)
}

#================================================================================
# DRIVER -- dispatch on mode: 'batch' (default, all groups) or 'pilot' (one group)
#================================================================================
if(sys.nframe()==0){   # only runs when the file is executed, not when sourced
  args   <- commandArgs(trailingOnly=TRUE)
  mode    <- if(length(args)>=1) args[1] else 'batch'
  dir.hab <- if(length(args)>=2) args[2] else
    "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/input_ascii_sum1/5min"
  dir.emp <- if(length(args)>=3) args[3] else
    "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/GFISHER/5min/maxn"
  dir.out <- if(length(args)>=4) args[4] else file.path(dir.emp,'affinity')
  file.env <- if(length(args)>=5) args[5] else
    "C:/Users/dchagaris/Github/WFS-FEM/GFISHER/data/April2026/env3LABS_93to24.csv"
  effort.control <- if(length(args)>=6) args[6] else 'both'
  if(!dir.exists(dir.out)) dir.create(dir.out, recursive=TRUE)

  cat("Loading habitat stack...\n")
  hab <- fn.load_habitat_stack(dir.hab)

  eff <- NULL
  if(effort.control!='none' && file.exists(file.env)){
    cat("Building survey-effort raster from", basename(file.env), "...\n")
    eff <- fn.build_effort_raster(file.env, hab[[1]],
             save.as=file.path(dir.out,'GFISHER_survey_effort_5min_66x78.asc'))
    cat("  surveyed cells:", sum(getValues(eff)>0, na.rm=TRUE),
        " total stations:", sum(getValues(eff), na.rm=TRUE), "\n")
  } else if(effort.control!='none'){
    warning("env file not found; proceeding with effort.control='none'")
    effort.control <- 'none'
  }

  if(mode=='batch'){
    cat("Estimating affinities for ALL groups (effort.control=", effort.control, ")...\n", sep='')
    out <- fn.batch_affinities(hab, dir.emp, dir.out, effort=eff, effort.control=effort.control)
    cat("\n--- AFFINITY TABLE (sorted by modnumber) ---\n")
    print(out, row.names=FALSE)
    cat("\nrho summary:\n"); print(summary(out$rho))

  } else if(mode=='pilot'){
    pilot <- if(length(args)>=7) args[7] else 'mod44_hogfish'
    f.emp <- list.files(dir.emp, pattern=paste0(pilot,'.*\\.asc$'), full.names=TRUE)
    if(length(f.emp)!=1) stop(paste("pilot empirical layer not uniquely found:", pilot))
    emp <- raster(f.emp[1])
    cat("Estimating affinities for pilot group:", pilot, "\n")
    res <- fn.estimate_habitat_affinities(hab, emp, effort=eff, effort.control=effort.control)
    cat("\n--- PILOT RESULT:", pilot, "---\n")
    cat("n valid cells:", res$n, "  rho(", res$control, "):", res$rho,
        "  plain Spearman:", res$spearman, "  Pearson r:", res$pearson, "\n")
    cat("Affinity A (max=1):\n"); print(res$A)
    png(file.path(dir.out, paste0('pilot_',pilot,'.png')), width=10, height=5, units='in', res=150)
    par(mfrow=c(1,2), mar=c(3,3,3,5))
    plot(emp, colNA='black', main=paste0('empirical MaxN\n',pilot))
    plot(res$pred, colNA='black', main=paste0('predicted (rho=',round(res$rho,3),')'))
    dev.off()
  } else stop("unknown mode: ", mode)
}
