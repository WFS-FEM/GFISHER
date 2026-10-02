# Issue #5: how much do the stanza-group affinities move with the stage 2 seed?
# Stages 2, 3 and 4a for seeds 1..10; everything written under the temp folder, nothing in the repo.
suppressPackageStartupMessages(library(raster))
setwd('C:/Repos/WFS-FEM/GFISHER')
source('R/video_dataset.R'); source('R/maxn_maps.R'); source('R/selection_ratio_affinities.R')
E <- 'C:/Users/User/AppData/Local/Temp/gfisher_baseline/exp_issue5'
note <- function(...) { msg <- paste0(format(Sys.time(), '%H:%M:%S'), '  ', paste0(...)); cat(msg, '\n'); cat(msg, '\n', file=file.path(E,'notes.log'), append=TRUE) }
d <- 'data/April2026'; bbox <- c(latN=30.5, latS=25, lonW=-87.5, lonE=-81)
depth <- raster('data/bathymetry/depth 5min 66x78.asc')
seeds <- 1:10; n.boot <- 1000
note('START issue #5 seed experiment: seeds ', paste(seeds, collapse=','), ', n.boot=', n.boot, ', basemaps = committed output/basemaps/5min')
hab <- fn.load_layer_stack('output/basemaps/5min')
eff <- fn.build_effort_raster(file.path(d,'env3LABS_93to24.csv'), hab[[1]])
note('habitat stack: ', nlayers(hab), ' layers; effort raster: ', sum(getValues(eff)>0, na.rm=TRUE), ' surveyed cells, ', sum(getValues(eff), na.rm=TRUE), ' stations')
for (s in seeds) {
  t0 <- Sys.time()
  m <- suppressMessages(fn.make_gfisher_videodataset(file.path(d,'maxn3LABS_93to24.csv'), file.path(d,'env3LABS_93to24.csv'), file.path(d,'lens3LABS_93to24.csv'),
         bbox, 'data/Master Species List.xlsx', col.modnum='modnum_mice', col.modname='modname_mice', col.fg='fg_mice', seed=s))
  dm <- file.path(E, sprintf('seed%02d', s), 'maps'); da <- file.path(E, sprintf('seed%02d', s), 'aff')
  dir.create(dm, recursive=TRUE, showWarnings=FALSE)
  invisible(capture.output(fn.make_GFISHER_maxn_maps(m, depth, plot=FALSE, fun=mean, background=0, dir.out=dm)))
  aff <- NULL; invisible(capture.output(aff <- fn.batch_selection_ratios(hab, dir.emp=dm, dir.out=da, effort=eff, n.boot=n.boot, con=list(apply=FALSE))))
  long <- aff$long; long$seed <- s
  write.csv(long, file.path(E, sprintf('long_seed%02d.csv', s)), row.names=FALSE)
  note(sprintf('seed %2d done in %4.1f min: stage2 rows=%d, total maxn=%d, affinity rows=%d', s, as.numeric(difftime(Sys.time(), t0, units='mins')), nrow(m), sum(m$maxn), nrow(long)))
}
# ---- summarise across seeds ----
L <- do.call(rbind, lapply(list.files(E, 'long_seed.*csv$', full.names=TRUE), read.csv, stringsAsFactors=FALSE))
L$species <- sub('[ -]?[0-9]+[+]?$', '', L$modname)
L$multistanza <- L$species %in% c('gag','red grouper')
agg <- function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x))
S <- aggregate(cbind(w, A) ~ modnumber + modname + code + family + multistanza, L, agg)
S <- do.call(data.frame, S)
# bootstrap half-width from seed 1 (each seed's bootstrap is itself seeded, so this is one realisation)
b1 <- L[L$seed==1, c('modnumber','code','w_lo','w_hi','n_pos','pooled_with')]
S <- merge(S, b1, by=c('modnumber','code'))
S$w_range <- S$w.max - S$w.min
S$boot_width <- S$w_hi - S$w_lo
S$seed_vs_boot <- round(S$w_range / S$boot_width, 2)   # <1: seed spread smaller than the bootstrap interval
S$A_range <- S$A.max - S$A.min
S <- S[order(S$modnumber, match(S$code, c('AL','AM','AH','NL','NM','NH','RCK','UNC','SGR'))), ]
write.csv(S, file.path(E, 'summary_across_seeds.csv'), row.names=FALSE)
note('SUMMARY written: summary_across_seeds.csv (', nrow(S), ' group x layer rows)')
note('control: single-stanza groups max w_range = ', signif(max(S$w_range[!S$multistanza]), 3), ', max A_range = ', signif(max(S$A_range[!S$multistanza]), 3), ' (must be 0)')
ms <- S[S$multistanza, ]
note('multistanza: median seed_vs_boot = ', round(median(ms$seed_vs_boot, na.rm=TRUE), 2), '; share of rows with seed range < bootstrap width = ', round(mean(ms$seed_vs_boot < 1, na.rm=TRUE), 2))
note('multistanza: A_range quantiles (0,50,90,100%) = ', paste(round(quantile(ms$A_range, c(0,.5,.9,1)), 3), collapse=' / '))
note('END')
