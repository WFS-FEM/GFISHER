# video_dataset.R -----------------------------------------------------------------
# Build the station x model-group MaxN table from the raw 3LABS video survey files.
# This is the input to R/maxn_maps.R and to the site-level affinity modules.
#
# Extracted verbatim from the former 'R/GFISHER functions.R' when the map-building code
# was split by stage; the estimator is unchanged.
#
# NOTE: the melt step drops maxn == 0, so the returned table holds PRESENCE RECORDS ONLY.
# Anything needing true zeros (a density rather than a conditional-on-presence mean) has to
# reinstate them against the full station list -- see fn.site_group_matrix in
# R/site_level_affinities.R.

suppressPackageStartupMessages({
  library('sp'); library('raster'); library('reshape2'); library('truncnorm')   # readxl is called by namespace
})

fn.make_gfisher_videodataset <- function(file.maxn, file.env, file.len, bbox, file.spplist,
                                         col.modnum='modnum', col.modname='modname', col.fg='fg',
                                         seed=1,
                                         file.lookup=file.path(dirname(dirname(file.maxn)),
                                                               paste0('GFISHER_species_fg_',col.modnum,'.csv'))){
  # seed: the multistanza step draws lengths (rtruncnorm) and pairs observed lengths with
  # individuals (sample.int), so stanza assignment is random. A fixed seed makes stages 2-4
  # reproducible run to run; NULL restores the unseeded behaviour. The draw moves individuals
  # between stanzas of a species only; totals per station and species are unaffected.
  # file.lookup: where the survey-taxon -> model-group key is written (NULL = not written). The
  # default keeps the historical location (the parent of the survey-data folder); the driver
  # points it at output/ instead so nothing is written next to the inputs.
  if(!is.null(seed)) set.seed(seed)
  # Species groupings are set by the caller, not hardcoded, via the three column arguments:
  #   col.modnum / col.modname : which spplist-sheet columns give the model-group number and name
  #   col.fg                   : which spp_stanzas_sizes column maps age stanzas to group numbers
  # Defaults reproduce the original Ecospace scheme (modnum/modname/fg). The MICE scheme is selected
  # by passing modnum_mice/modname_mice/fg_mice. See process GFISHER data.R for the driver knob.
  
  #fxn arguments/inputs
  # bbox <- bbox
  # spplist <- spplist
  
  #--------------------------------import and prepare data------------------------------------------
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
  # Read the whole spplist sheet (colIndex wide enough to include the MICE columns at ~44-45), then
  # promote the requested scheme's columns to the generic names 'modnumber'/'modname' used downstream.
  modspp.raw = as.data.frame(readxl::read_excel(file.spplist, sheet="spplist"))   # readxl: no Java needed
  modspp.raw = modspp.raw[!is.na(modspp.raw$modnum),]
  # crosswalk ORIGINAL modnum -> scheme group number, for the hardcoded taxon fallbacks further down
  cw = unique(modspp.raw[,c('modnum', col.modnum)]); names(cw) = c('orig','scheme')
  modspp = modspp.raw
  modspp$modnumber = modspp.raw[[col.modnum]]
  modspp$modname   = modspp.raw[[col.modname]]
  spplist <- melt(modspp, id.vars=c('modnumber','modname'), measure.vars=c('species','query','og_name','class','order','family','genus'), variable.name='var',value.name='taxon')
  spplist <- spplist[,-3]
  spplist <- data.frame(lapply(spplist, tolower), stringsAsFactors = FALSE)
  spplist$taxon <- ifelse(spplist$taxon=="",NA,spplist$taxon)
  spplist <- unique(spplist)
  spplist <- spplist[!spplist$modnumber %in% c('99','0'),]   # drop 'not included' (99) and 'omitted' (0)
  #spplist <- spplist[complete.cases(spplist),]
  # number -> name lookup for the requested scheme (replaces the original-only 'model groups' sheet)
  modgrps <- unique(spplist[!is.na(spplist$modnumber), c('modnumber','modname')])

  ##multistanza size at age-------------------------------------------------------------------------
  sizeatage <- as.data.frame(readxl::read_excel(file.spplist, sheet="spp_stanzas_sizes"))
  
  #keep species included in the model
  keeptaxa = tolower(sort(unique(spplist$taxon)))
  spp.gfsh = data.frame(taxon=sort(unique(names(dat.maxn)[-c(1:2)])),
                        taxon2 = sub("_sp$","",sort(unique(names(dat.maxn)[-c(1:2)]))),
                        taxon3 = sapply(strsplit(sort(unique(names(dat.maxn)[-c(1:2)])), "_"), `[`, 1))
  spp.gfsh$match = ifelse(spp.gfsh$taxon %in% keeptaxa | spp.gfsh$taxon2 %in% keeptaxa | spp.gfsh$taxon3 %in% keeptaxa,TRUE,FALSE)
  spp.gfsh$modnumber = ifelse(spp.gfsh$taxon %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon,spplist$taxon)],
                              ifelse(spp.gfsh$taxon2 %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon2,spplist$taxon)],
                                     ifelse(spp.gfsh$taxon3 %in% keeptaxa, spplist$modnumber[match(spp.gfsh$taxon3,spplist$taxon)],NA)))
  
  # ambiguous family/genus catch-alls, given as ORIGINAL modnumbers (19='other snapper', 36='other SWG')
  # and translated to the requested scheme via the crosswalk so they land in the right group.
  spp.gfsh$modnumber = ifelse(spp.gfsh$taxon=='lutjanidae_sp', as.character(cw$scheme[match(19, cw$orig)]),
                              ifelse(spp.gfsh$taxon %in% c('epinephelus_sp','epinephelus_striatus'), as.character(cw$scheme[match(36, cw$orig)]),spp.gfsh$modnumber))
  spp.gfsh$modname = modgrps$modname[match(spp.gfsh$modnumber, modgrps$modnumber)]

  if(!is.null(file.lookup)){
    if(!dir.exists(dirname(file.lookup))) dir.create(dirname(file.lookup), recursive=TRUE)
    write.csv(spp.gfsh, file.lookup, row.names=FALSE)
  }
  
  #melt video data----------------------------------------------------------------------------------
  dat.maxn.long <- melt(dat.maxn[,-1],id.vars='reference', variable.name='sciname', value.name='maxn')
  dat.maxn.long <- dat.maxn.long[dat.maxn.long$maxn>0,]
  dat.maxn.long$sciname <- gsub("_"," ",as.character(dat.maxn.long$sciname))
  
  #get size data for multistanza species------------------------------------------------------------
  # A species is multistanza UNDER THIS SCHEME only if its scheme fg column lists more than one group
  # (e.g. gag -> '4-5-6-7-8-9' under MICE). Species whose stanzas collapse to a single group under the
  # scheme (e.g. red snapper -> 'snappers') fall through to the spplist assignment below instead. This
  # replaces the old hardcoded [c(1,3,4)] index, and reproduces it exactly for the original scheme.
  fg.scheme <- as.character(sizeatage[[col.fg]])
  ntok <- sapply(strsplit(fg.scheme, "-"), function(z) sum(nzchar(z)))
  is.multi <- !is.na(sizeatage$stanzas) & sizeatage$stanzas!="0" & ntok > 1
  multistanza.fg  <- as.numeric(unlist(strsplit(fg.scheme[is.multi], "-")))
  multistanza.spp <- tolower(sizeatage$sciname[is.multi])
  dat.lf$sciname <- tolower(gsub("_"," ",dat.lf$sciname))
  
  lf2 = dat.lf[dat.lf$sciname %in% multistanza.spp,c('reference','sciname','length_mm')]
  lf.spp.mean <- aggregate(length_mm~sciname, lf2, mean)
  lf.spp.sd <- aggregate(length_mm~sciname, lf2, sd)
  lf.spp.min <- aggregate(length_mm~sciname, lf2, min)
  lf.spp.max <- aggregate(length_mm~sciname, lf2, max)
  lf.spp.cnt <- aggregate(length_mm~sciname, lf2, length)
  lf.spp.sd$length_mm[is.na(lf.spp.sd$length_mm)] <- lf.spp.mean$length_mm[is.na(lf.spp.sd$length_mm)]* mean(lf.spp.sd[,2]/lf.spp.mean[,2],na.rm=T)

  lf3 <- merge(lf2, dat.maxn.long)

  dat3 <- dat.maxn.long[,c('reference','sciname','maxn')]
  dat3$sciname <- gsub("_"," ",dat3$sciname)
  dat3 <- dat3[dat3$sciname %in% multistanza.spp,]
  # every multistanza species with MaxN records needs length records to draw from
  no.len <- setdiff(unique(dat3$sciname), lf.spp.mean$sciname)
  if(length(no.len)) stop('multistanza species with MaxN records but no length records: ',
                          paste(no.len, collapse=', '))

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
    laa.i = as.numeric(unlist(sizeatage[which(tolower(sizeatage$sciname)==spp.i)[1],which(substr(names(sizeatage),1,3)=='age')]))
    stanzas.i = as.numeric(unlist(strsplit(sizeatage$stanzas[tolower(sizeatage$sciname)==spp.i],"-")))
    groups.i = as.numeric(unlist(strsplit(sizeatage[[col.fg]][tolower(sizeatage$sciname)==spp.i],"-")))
    size.i = dat4$len_mm[i]/10
    age.i = which.min(abs(laa.i-size.i))-1
    stz.i = tail(which(age.i-stanzas.i>=0),1)
    # an age below the first stanza boundary would otherwise give a zero-length index and the
    # cryptic 'replacement has length zero'; fail with the facts instead
    if(length(stz.i)==0) stop(sprintf('%s: length %d mm -> age %d is below the first stanza (%s)',
                                      spp.i, dat4$len_mm[i], age.i, sizeatage$stanzas[tolower(sizeatage$sciname)==spp.i][1]))
    grp.i = groups.i[stz.i]
    grpname.i = modgrps$modname[match(grp.i, modgrps$modnumber)]
    
    dat4$modnumber[i] <- grp.i
    dat4$modname[i] <- grpname.i
    
  }
  dat4$n_at_length = 1
  dat4 <- aggregate(n_at_length~reference+sciname+modnumber+modname, data=dat4, sum)
  dat4$maxn <- dat4$n_at_length
  dat4$n_at_length <- NULL
  
  #--------------------------assign NON-multistanza species to fg -------------------------------------
  #resume here...
  dat5 <-  dat.maxn.long[!dat.maxn.long$sciname %in% multistanza.spp,c('reference','sciname','maxn')]
  spp.key <- spp.gfsh[,c('taxon','modnumber','modname')]
  spp.key$taxon <- gsub("_"," ", spp.key$taxon)          # match dat5$sciname format (spaces, not underscores)
  dat5 <- merge(dat5, spp.key, by.x='sciname',by.y='taxon', all.x=T)
  
  
  #put it back together
  dat.full <- rbind(dat4, dat5)

  #check counts: every raw MaxN must be accounted for after the stanza split and group assignment.
  #A mismatch means records were lost, so stop rather than print a message that scrolls past.
  sum1 <- aggregate(maxn~sciname, data=dat.full, sum)
  names(sum1)[2] <- 'final_maxn'
  sum2 <- aggregate(maxn~sciname, data=dat.maxn.long, sum)
  names(sum2)[2] <- 'raw_maxn'
  chk <- merge(sum1,sum2)
  err <- which(chk$final_maxn != chk$raw_maxn)
  if(length(err)>0){
    stop('MaxN counts did not sum back up after processing for: ',
         paste(sprintf('%s (%d raw -> %d final)', chk$sciname[err], chk$raw_maxn[err], chk$final_maxn[err]), collapse='; '))
  }

  #merge back with env data
  dat.full <- merge(dat.full, stations[,-which(names(stations) %in% c('dup','complete'))], all.x=T, by='reference')
  return(dat.full)
} #eof
