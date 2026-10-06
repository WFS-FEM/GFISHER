# _setup.R -------------------------------------------------------------------------
# Session setup for `process GFISHER data.R`: package check, plot device, input manifest.
# Base R only, so it can run before any package is attached. Sourced once at the top of
# the driver, after the repo-root check. Mirrors scripts/_setup.R in RedTideMaps and
# R/data_setup_functions.R in EcospaceBasemap. The inputs EcospaceBasemap produces (raw
# dbSEABED grids, seagrass rasters, depth values) ship here as copies with MD5s in each
# data/*/SOURCE.md; the download and rasterisation code lives only in EcospaceBasemap.

# Packages ---------------------------------------------------------------------------
# Required by the four live stages (R/habitat_basemaps.R, video_dataset.R, maxn_maps.R and
# the three affinity modules). Optional ones are only needed for non-default paths and are
# reported, not required: mgcv for fn.make_habitat_basemaps(target='smooth').
GFISHER.PACKAGES <- c('sf', 'sp', 'raster', 'FNN', 'colorRamps', 'reshape2', 'truncnorm', 'readxl')
GFISHER.PACKAGES.OPTIONAL <- c('mgcv')

#' Stop early, with one install line, if any required package is missing.
#' Runs before any library() call so a missing package fails in the first second of a run
#' rather than partway through stage 2.
fn.check_packages <- function(pkgs = GFISHER.PACKAGES, optional = GFISHER.PACKAGES.OPTIONAL){
  have <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)
  missing <- pkgs[!have(pkgs)]
  if(length(missing) > 0)
    stop('Missing R package(s): ', paste(missing, collapse = ', '),
         '\nInstall with:\n  install.packages(c(', paste0('"', missing, '"', collapse = ', '), '))',
         call. = FALSE)
  opt.missing <- optional[!have(optional)]
  if(length(opt.missing) > 0)
    message('Optional package(s) not installed (only needed for non-default options): ',
            paste(opt.missing, collapse = ', '))
  invisible(TRUE)
}

# Plot device ------------------------------------------------------------------------
#' Open a recording graphics window when that makes sense, and do nothing otherwise.
#' Replaces the old Windows-only recording-device call, which failed under Rscript and on
#' non-Windows systems and warns when .SavedPlots does not exist.
fn.plot_device <- function(){
  if(interactive() && .Platform$OS.type == 'windows'){
    try(rm(.SavedPlots, envir = .GlobalEnv), silent = TRUE)   # RGui plot history
    try(dev.new(record = TRUE), silent = TRUE)
  }
  invisible(NULL)
}

# Input discovery --------------------------------------------------------------------
#' The GFISHER geodatabase inside `dir`, when there is exactly one. Returns NA (not an
#' error) when none is found so fn.check_inputs can report it with the others.
fn.find_gdb <- function(dir, pattern = '^GFISHER_EAST_Universe.*\\.gdb$'){
  if(!dir.exists(dir)) return(NA_character_)
  g <- list.files(dir, pattern = pattern, full.names = TRUE, include.dirs = TRUE)
  g <- g[dir.exists(g)]
  if(length(g) == 1) return(g)
  if(length(g) > 1) stop('More than one geodatabase matches ', pattern, ' in ', dir,
                         '; set file.gdb in config.local.R.', call. = FALSE)
  NA_character_
}

#' The four raw dbSEABED grids the basemap and substrate stages read.
fn.dbseabed_files <- function(dir.dbseabed, codes = c('RCK', 'GVL', 'SND', 'MUD'))
  file.path(dir.dbseabed, paste0('Gmf_', codes), paste0('gmf_', codes, '_val.asc'))

# Input manifest ---------------------------------------------------------------------
#' One row per input the driver reads, with its resolved path, which stages need it,
#' whether it is required, and where it comes from. `how` is one of:
#'   repo    ships with the repository
#'   auto    downloaded by the driver on first run
#'   manual  must be obtained by hand (FWRI data that cannot be redistributed on GitHub)
#'   derived produced by a sibling repo (EcospaceBasemap) or supplied by the author
fn.data_manifest <- function(dir.data, file.gdb, file.spplist, file.depth, dir.dbseabed,
                             file.seagrass, res){
  gdb <- if(is.null(file.gdb) || is.na(file.gdb)) file.path(dir.data, 'GFISHER_EAST_Universe_<year>.gdb') else file.gdb
  rbind(
    data.frame(key = 'bathymetry', stage = '1,3,4', required = TRUE, how = 'repo',
               path = file.depth,
               source = 'data/bathymetry/ (ships with the repo; made by EcospaceBasemap)'),
    data.frame(key = 'species_list', stage = '2', required = TRUE, how = 'repo',
               path = file.spplist,
               source = 'data/Master Species List.xlsx (ships with the repo)'),
    data.frame(key = 'geodatabase', stage = '1,4b,4c', required = TRUE, how = 'manual',
               path = gdb,
               source = 'FWRI GFISHER side-scan habitat mapping (Sean Keenan), ~275 MB. Request from FWRI or the repo author; place in dir.data or set file.gdb.'),
    data.frame(key = 'survey_maxn', stage = '2', required = TRUE, how = 'manual',
               path = file.path(dir.data, 'maxn3LABS_93to24.csv'),
               source = 'FWRI 3LABS video survey 1993-2024 (maxn), ~23 MB. Request with the geodatabase; place in dir.data.'),
    data.frame(key = 'survey_env', stage = '2,4', required = TRUE, how = 'manual',
               path = file.path(dir.data, 'env3LABS_93to24.csv'),
               source = 'FWRI 3LABS station/environment file, ~6 MB. Same source.'),
    data.frame(key = 'survey_lens', stage = '2', required = TRUE, how = 'manual',
               path = file.path(dir.data, 'lens3LABS_93to24.csv'),
               source = 'FWRI 3LABS length file, ~13 MB. Same source.'),
    data.frame(key = 'dbseabed', stage = '1,4c', required = TRUE, how = 'repo',
               path = dir.dbseabed,
               source = 'CSDMS dbSEABED raw grids, tracked copy in data/dbseabed/ (MD5s in its SOURCE.md). If removed: git checkout -- data/dbseabed, or set dir.dbseabed (or dir.ecospace.basemap) to an EcospaceBasemap clone; EcospaceBasemap fn.pull_dbseabed() is the download path.'),
    data.frame(key = 'seagrass', stage = '1', required = FALSE, how = 'repo',
               path = file.seagrass,
               source = paste0('Seagrass raster on the ', res, '-min grid. The 5- and 15-min ones ship in data/seagrass/ (copied from EcospaceBasemap; MD5s in SOURCE.md), or set dir.ecospace.basemap to read a clone. Other resolutions must be made there. Without it SGR = 0.')),
    stringsAsFactors = FALSE)
}

#' Check every input before any long computation and stop with instructions if a
#' required one is missing. Prints a status table so a run log records what was used.
fn.check_inputs <- function(manifest, stop.on.missing = TRUE){
  m <- manifest
  m$present <- vapply(seq_len(nrow(m)), function(i){
    p <- m$path[i]
    if(is.na(p) || !nzchar(p)) return(FALSE)
    if(m$key[i] == 'dbseabed') return(all(file.exists(fn.dbseabed_files(p))))
    file.exists(p)
  }, logical(1))

  cat('\nInput check\n', strrep('-', 100), '\n', sep = '')
  cat(sprintf('  %-13s %s\n', 'basemap src',
              if(exists('dir.ecospace.basemap') && !is.null(dir.ecospace.basemap))
                paste0('EcospaceBasemap clone at ', dir.ecospace.basemap)
              else 'copies shipped in data/ (dir.ecospace.basemap not set)'))
  cat(sprintf('  %-13s %-8s %-8s %-8s %s\n', 'input', 'stage', 'need', 'how', 'status / path'))
  for(i in seq_len(nrow(m)))
    cat(sprintf('  %-13s %-8s %-8s %-8s %s  %s\n', m$key[i], m$stage[i],
                if(m$required[i]) 'required' else 'optional', m$how[i],
                if(m$present[i]) 'OK     ' else 'MISSING', m$path[i]))
  cat(strrep('-', 100), '\n', sep = '')

  miss <- m[!m$present, ]
  if(nrow(miss) > 0){
    cat('\nMissing input(s):\n')
    for(i in seq_len(nrow(miss)))
      cat(sprintf('  %-13s %s\n', miss$key[i], miss$source[i]))
    cat('  See README.md > Getting the data, and config.local.example.R for the path overrides.\n')
  }

  short <- m$key[!m$present & m$required]
  if(length(short) > 0){
    msg <- paste0('Missing required input(s): ', paste(short, collapse = ', '), '.')
    if(isTRUE(stop.on.missing)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(m)
}
