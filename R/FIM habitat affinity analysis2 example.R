rm(list=ls());rm(.SavedPlots);gc();graphics.off();windows(record=T)
#.libPaths("C:\\R\\win-library")
library('raster')
install.packages('rgdal')
library('USAboundaries')
library('gridExtra')
library('concaveman')
library('dplyr')
library('sf')

################################################################################
#           setup, read and prepare data
################################################################################
setwd("C:\\Users\\dchagaris\\OneDrive - University of Florida\\Suwannee NAS Project\\Ecosystem Modelling\\SREM\\FIM\\habitat affinity analysis")

#read fim data
dir.fim = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Suwannee NAS Project\\Ecosystem Modelling\\SREM\\FIM\\SREM FIM 2020"
filename.fim = paste0(dir.fim,"/data/FIM 3 biomassses QA/FIM CK biomass density FG 1996-2020 QA 20230808.csv")
dat = read.csv(filename.fim, stringsAsFactors = F)
dat$date = as.Date(dat$date,format="%d%b%Y")
names(dat)
phys = unique(dat[,c(1:24,28:32)])

################################################################################
#           logistic habitat preference models
################################################################################
#clasify each sample as a single habitat type-----------------------------------
table(phys$bveg); table(phys$bottom); table(phys$shore); table(phys$zone)
phys$habitat = NA
phys$habitat = ifelse(phys$bveg%in%c('SAV','SAVAlg') & phys$bottomvegcover>=30,'seagrass',
                ifelse(phys$shore=='oyster'& phys$bottom=='oyster','oyster',
                  ifelse(phys$shore=='marsh_mangrove'& phys$bottom!='oyster','marsh',NA)))

phys$habitat = ifelse(is.na(phys$habitat),
                 ifelse(phys$bveg=='SAV','seagrass',
                 ifelse(phys$bottom=='oyster','oyster',
                 ifelse(phys$shore=='oyster','oyster',
                 ifelse(phys$shore=='marsh_mangrove','marsh',
                 ifelse(phys$zone=='F' | phys$salinity<3,'river','base'))))),
               phys$habitat)

#check habitat assignments------------------------------------------------------
table(phys$habitat,useNA='always')
table(phys$bveg,phys$habitat,useNA='always')
table(phys$bottom,phys$habitat,useNA='always')
table(phys$shore,phys$habitat,useNA='always')
table(phys$bveg2,phys$habitat,useNA='always')
table(phys$bottom2,phys$habitat,useNA='always')
table(phys$shore2,phys$habitat,useNA='always')

#bring habitat variable into full dataset---------------------------------------
dat$habitat = phys$habitat[match(dat$reference,phys$reference)]

#loop through species and fit logistic regression-------------------------------
modgroups = sort(unique(dat$functional.group))
hpf.out = data.frame(matrix(NA,nrow=length(modgroups),ncol=length(unique(dat$habitat))))
rownames(hpf.out) = modgroups; names(hpf.out) = sort(unique(dat$habitat))
for(g in 1:length(modgroups)){
  #g=20
  modname = modgroups[g] 
  dsub = subset(dat,functional.group==modname)
  dsub.pos = subset(dsub,response>0)
  
  #check factor levels
  t.gr = prop.table(table(dsub.pos$gr2))
  t.zn = prop.table(table(dsub.pos$zone))

  #keep factor levels
  keep.gr = names(t.gr)[which(t.gr>0.01)]
  keep.zn = names(t.zn)[which(t.zn>0.01)]

  #subset and make positive only and binomial datasets
  dsub2 = subset(dsub,gr2%in%keep.gr & zone%in%keep.zn) 
  d.pos = subset(dsub2,response>0)
  d.pos$response = log(d.pos$response+1)
  d.bin = dsub2
  d.bin$response = ifelse(d.bin$response>0,1,0)
  
  #fit logistic model and predict probabilities for each habitat
  hpf.mod = glm(response~habitat,data=d.bin,family=binomial(link='logit'))
  hpf.pred = data.frame(habitat=sort(unique(d.bin$habitat)))
  hpf.pred$pred = predict(hpf.mod, newdata = hpf.pred, type = "response")
  
  #store output
  hpf.out[g,] = hpf.pred$pred[match(names(hpf.out),hpf.pred$habitat)]
}
getwd()
#write.csv(hpf.out,'FIM logistic habitat preferences 20230808.csv',row.names=T)

#hpf.out = read.csv('FIM logistic habitat preferences.csv')
#names(hpf.out)[1] = 'FG'
poolcodes = read.csv(paste0(dir.fim,"/data/basic inputs/poolcodes.csv"))
poolcodes$FG2 = tolower(poolcodes$FG2)

hpf.out2 = merge(poolcodes[,c(3,1)],hpf.out,by.x='FG2',by.y=0,all.x=T)
hpf.out2 = hpf.out2[order(hpf.out2$poolcode),c(2,1,5,4,7,6,3)]

write.csv(hpf.out2,paste0('FIM logistic habitat preferences ',gsub(":","",Sys.time()),'.csv'),row.names=F)

################################################################################
# get polygon for fim samples
################################################################################
#create polygon of fim sample area----------------------------------------------
#there was a typo in a longitude entry, remove it
phys2 = phys[-which.min(phys$longitude),]
coordinates(phys2) = ~longitude+latitude
crs(phys2) = "+proj=robin"
proj4string(phys2) = CRS("+proj=longlat +datum=WGS84")
phys2 = spTransform(phys2,CRS="+proj=utm +zone=16 +datum=WGS84 +units=m")
#fwri.grid = projectRaster(depth.4k,crs="+proj=utm +zone=16 +datum=WGS84 +units=m")

#make a polygon
fim.pnts <- phys2 %>% st_as_sf(coords = c("longitude", "latitude"),crs=4326)
fim.poly <- concaveman(fim.pnts,concavity = 4)
plot(fim.poly, reset = T)
plot(fim.pnts, add = TRUE)

#convert to sp object
fim.poly2 = as(fim.poly,'Spatial')
plot(fim.poly2)

#area of fim universe
fim.area.m2 = raster::area(fim.poly2)
getwd()
sf::st_write(obj=fim.poly,dsn=paste0(getwd(),'/fim_poly.kml'),driver='kml')

################################################################################
# area of habitat
################################################################################
#From website: All data are provided in the Florida State Plane Projection, North 
#Zone, Units US Feet, Datum HPGN (83/90). In some cases data may be Florida 
#State Plane Projection, North Zone, Units US Feet, Datum NAD83. Currently, 
#most data are available as district wide layers. No customized clipping 
#or joining will be performed.

#get and set proj4string based on description above; FL state plan projection is 
#in the USAboundaries package data
state_proj = data.frame(state_proj)
state_proj$proj4_string[state_proj$state=='FL' & state_proj$zone=='north']
FLproj = state_proj$proj4_string[state_proj$state=='FL' & state_proj$zone=='north']
#change units to ft
FLproj = "+proj=lcc +lat_1=30.75 +lat_2=29.58333333333333 +lat_0=29 +lon_0=-84.5 +x_0=600000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=ft +no_defs"

#store results
hab.prop.out = data.frame()

#base directory for habitats
dir.hab = 'C:\\Users\\dchagaris\\OneDrive - University of Florida\\Suwannee NAS Project\\Ecosystem Modelling\\SREM\\habitat'

#seagrass-----------------------------------------------------------------------
hab.shp <- shapefile(paste0(dir.hab,'\\seagrass\\seagrass2001_SRWMD\\seagrass01.shp'))
hab.label = 'Seagrass'

#define projection of shapefile
proj4string(hab.shp) = CRS(FLproj)

#get seagrass and fim polys on same coordinate reference system.................
if(hab.label=='Seagrass'){
  hab2 = subset(hab.shp,substr(hab.shp$FLUCSDESC,1,8)=='Seagrass')
} else{
  hab2 = hab.shp
}
hab2$ID = 1:nrow(hab2)# sg2$SRW_SAV01_
hab2 = spTransform(hab2,CRS=projection(fim.poly))

plot(hab2,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),col='green',border='darkgreen',axes=T)
plot(fim.poly,add=T,border='blue',lwd=2)

#determine overlap of habitat polygons with fim sampling universe..............
#test: polygon 573 is fully inside the fim area, polygon 1 is fully outside
hab2 = st_as_sf(hab2)
area.out = data.frame()
for(i in 1:nrow(hab2)){
  #i=227
  hab.sub = hab2[i,]
  #plot(sg.sub,ylim=c(min(extent(sg2)[3],extent(fim.poly)[3]),max(extent(sg2)[4],extent(fim.poly)[4])),col='green',axes=T,main=i)
  #plot(fim.poly,add=T,border='red')
  polyint = st_intersection(hab.sub,fim.poly)
  polyint.area <- polyint %>% mutate(area = st_area(.) %>% as.numeric())
  if(nrow(polyint.area)==0){
    area.tmp = data.frame(id=i,area.tot.ft2=hab.sub$AREA, area.tot.m2=0.092903*hab.sub$AREA, area.olap.m2=0, geometry=hab.sub$geometry)
  } else{
    area.tmp = data.frame(id=i,area.tot.ft2=hab.sub$AREA, area.tot.m2=0.092903*hab.sub$AREA, area.olap.m2=polyint.area$area, geometry=polyint.area$geometry)
  }
  area.out = rbind(area.out,area.tmp)
}

#get total habitat area and percent coverage
hab.area.m2 = sum(area.out$area.olap.m2)
hab.prop = hab.area.m2/fim.area.m2

#summary plot
tiff(paste0(getwd(),"/",hab.label,' coverage map.tiff'),width=5,height=5,units='in',res=300,compression='lzw')
plot(hab2$geometry,col='gray',border='gray',axes=T,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),
  main=hab.label)
  plot(area.out$geometry[area.out$area.olap.m2>0],col='green',border='green',add=T)
  plot(fim.poly,border='blue',lwd=2,add=T)
  text(850000,3220000,paste0(round(hab.prop*100,2),' % coverage'),adj=0,cex=1.5)
dev.off()

#store result
hab.prop.out = rbind(hab.prop.out,data.frame(habitat=hab.label,proportion=hab.prop))

#oyster-------------------------------------------------------------------------
dir.hab = 'C:\\Users\\dchagaris\\OneDrive - University of Florida\\Suwannee NAS Project\\Ecosystem Modelling\\SREM\\habitat'
hab.shp <- shapefile(paste0(dir.hab,'\\Oyster\\oyster2001_SRWMD\\oyster01.shp'))
hab.label = 'Oyster'
#area is in feet
names(hab.shp)
str(hab.shp$AREA)
sum(hab.shp$ACRES)
hab.shp$AREA/hab.shp$ACRES
tapply(hab.shp$ACRES,hab.shp$HABITAT,sum)


#define projection of shapefile
proj4string(hab.shp) = CRS(FLproj)

#get habitat and fim polys on same coordinate reference system.................
if(hab.label=='Seagrass'){
  hab2 = subset(hab.shp,substr(hab.shp$FLUCSDESC,1,8)=='Seagrass')
} else{
  hab2 = hab.shp
}
hab2$ID = 1:nrow(hab2)# sg2$SRW_SAV01_
hab2 = spTransform(hab2,CRS=projection(fim.poly))

plot(hab2,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),col='green',border='darkgreen',axes=T)
plot(fim.poly,add=T,border='blue',lwd=2)

#determine overlap of habitat polygons with fim sampling universe..............
hab2 = st_as_sf(hab2)
area.out = data.frame()
for(i in 1:nrow(hab2)){
  #i=227
  print(paste0(i,'/',nrow(hab2)));flush.console();
  hab.sub = hab2[i,]
  #plot(sg.sub,ylim=c(min(extent(sg2)[3],extent(fim.poly)[3]),max(extent(sg2)[4],extent(fim.poly)[4])),col='green',axes=T,main=i)
  #plot(fim.poly,add=T,border='red')
  polyint = st_intersection(hab.sub,fim.poly)
  polyint.area <- polyint %>% mutate(area = st_area(.) %>% as.numeric())
  if(nrow(polyint.area)==0){
    area.tmp = data.frame(id=i,area.tot.ft2=hab.sub$AREA, area.tot.m2=0.092903*hab.sub$AREA, area.olap.m2=0, geometry=hab.sub$geometry)
  } else{
    area.tmp = data.frame(id=i,area.tot.ft2=hab.sub$AREA, area.tot.m2=0.092903*hab.sub$AREA, area.olap.m2=polyint.area$area, geometry=polyint.area$geometry)
  }
  area.out = rbind(area.out,area.tmp)
}

#get total habitat area and percent coverage
sum(area.out$area.tot.ft2)
sum(area.out$area.tot.m2)
sum(area.out$area.olap.m2)
hab.area.m2 = sum(area.out$area.olap.m2)
hab.prop = hab.area.m2/fim.area.m2

#summary plot
tiff(paste0(getwd(),"/",hab.label,' coverage map.tiff'),width=5,height=5,units='in',res=300,compression='lzw')
plot(hab2$geometry,col='gray',border='gray',axes=T,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),
     main=hab.label)
plot(area.out$geometry[area.out$area.olap.m2>0],col='green',border='green',add=T)
plot(fim.poly,border='blue',lwd=2,add=T)
text(850000,3220000,paste0(round(hab.prop*100,2),' % coverage'),adj=0,cex=1.5)
dev.off()

#store result
hab.prop.out = rbind(hab.prop.out,data.frame(habitat=hab.label,proportion=hab.prop))

#saltmarsh----------------------------------------------------------------------
#hab.shp <- shapefile(paste0(dir.hab,'\\saltmarsh\\Suwanee_River_Water_Management_District_(SRWMD)_2016-2017_Land_Use.shp'))
hab.shp <- readOGR(paste0(dir.hab,'\\saltmarsh\\Suwanee_River_Water_Management_District_(SRWMD)_2016-2017_Land_Use.shp'))
hab.label = 'Saltmarsh'
projection(hab.shp)

#define projection of shapefile
#proj4string(hab.shp) = CRS(FLproj)

if(hab.label=='Seagrass'){
  hab2 = subset(hab.shp,substr(hab.shp$FLUCSDESC,1,8)=='Seagrass')
} else if(hab.label=='Saltmarsh'){
  hab2 = hab.shp[hab.shp$LEVEL_3_DE=='Saltwater Marshes',]
}else{
  hab2 = hab.shp
}


hab2$ID = 1:nrow(hab2)# sg2$SRW_SAV01_
hab2 = spTransform(hab2,CRS=projection(fim.poly))
projection(hab2)

plot(hab2,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),col='green',border='darkgreen',axes=T)
plot(fim.poly,add=T,border='blue',lwd=2)

#determine overlap of habitat polygons with fim sampling universe..............
hab2 = st_as_sf(hab2)
projection(hab2)
area.out = data.frame()
for(i in 1:nrow(hab2)){
  #i=1
  hab.sub = hab2[i,]
  #plot(sg.sub,ylim=c(min(extent(sg2)[3],extent(fim.poly)[3]),max(extent(sg2)[4],extent(fim.poly)[4])),col='green',axes=T,main=i)
  #plot(fim.poly,add=T,border='red')
  polyint = st_intersection(hab.sub,fim.poly)
  polyint.area <- polyint %>% mutate(area = st_area(.) %>% as.numeric())
  if(nrow(polyint.area)==0){
    area.tmp = data.frame(id=i,area.tot.ft2=NA, area.tot.m2=hab.sub$SHAPEAREA, area.olap.m2=0, geometry=hab.sub$geometry)
  } else{
    area.tmp = data.frame(id=i,area.tot.ft2=NA, area.tot.m2=hab.sub$SHAPEAREA, area.olap.m2=polyint.area$area, geometry=polyint.area$geometry)
  }
  area.out = rbind(area.out,area.tmp)
}

#get total habitat area and percent coverage
hab.area.m2 = sum(area.out$area.olap.m2)
hab.prop = hab.area.m2/fim.area.m2

#summary plot
tiff(paste0(getwd(),"/",hab.label,' coverage map.tiff'),width=5,height=5,units='in',res=300,compression='lzw')
plot(hab2$geometry,col='gray',border='gray',axes=T,ylim=c(min(extent(hab2)[3],extent(fim.poly)[3]),max(extent(hab2)[4],extent(fim.poly)[4])),
     main=hab.label)
plot(area.out$geometry[area.out$area.olap.m2>0],col='green',border='green',add=T)
plot(fim.poly,border='blue',lwd=2,add=T)
text(800000,3220000,paste0(round(hab.prop*100,2),' % coverage'),adj=0,cex=1)
dev.off()

#store result
hab.prop.out = rbind(hab.prop.out,data.frame(habitat=hab.label,proportion=hab.prop))

write.csv(hab.prop.out,'habitat proportions.csv',row.names=F)

#simple habitat affinity calculation (Monaco et al 2002)------------------------
fimpref = read.csv('FIM logistic habitat preferences.csv')
names(fimpref)[3] = 'saltmarsh'

habprop = read.csv('habitat proportions.csv')
habprop = rbind(habprop,data.frame(habitat='Base',proportion=1-sum(habprop$proportion)))
habprop2 = habprop$proportion
names(habprop2) = tolower(habprop$habitat)
habprop2 = habprop2[match(names(habprop2),names(fimpref)[-1])]

hai = data.frame(matrix(NA,nrow=nrow(fimpref),ncol=length(habprop2)+1))
names(hai) = c('group',names(habprop2))
odds = chess = hai
for(i in 1:nrow(fimpref)){
  #i=1
  tmp.u = as.numeric(fimpref[i,-1]);names(tmp.u) = names(fimpref)[-1]
  tmp.u2 = tmp.u/sum(tmp.u)
  odds[i,1] = hai[i,1] = chess[i,1] = fimpref$group[i]
  for(j in 1:length(tmp.u)){
    #j=2
    hai[i,j+1] = ifelse(tmp.u[j]<=habprop2[j],(tmp.u[j]-habprop2[j])/habprop2[j],
                      (tmp.u[j]-habprop2[j])/(1-habprop2[j]))
    chess[i,j+1] = tmp.u[j]/habprop2[j]
    Aij = habprop2[j]/(1-habprop2[j])
    Gij = tmp.u2[j]/(1-tmp.u2[j])
    Oij = Gij/Aij
    Xij = Oij/(1+Oij)
    odds[i,j+1] = Xij
  }
}
chess[,-1] = prop.table(as.matrix(chess[,2:ncol(chess)]),1)

hai2 = hai
hai2[hai2<0] = 0
hai2[,2:ncol(hai2)] = prop.table(as.matrix(hai2[,2:ncol(hai2)]),1)

hai3 = hai
hai3[,-1] = hai3[,-1]+1
hai3[,2:ncol(hai3)] = prop.table(as.matrix(hai3[,2:ncol(hai3)]),1)
  


#Manly resource selection ratio-------------------------------------------------
library('adehabitatHS')
fimpref = read.csv('FIM logistic habitat preferences.csv')
names(fimpref)[3] = 'saltmarsh'

habprop = read.csv('habitat proportions.csv')
habprop = rbind(habprop,data.frame(habitat='Base',proportion=1-sum(habprop$proportion)))
habprop2 = habprop$proportion
names(habprop2) = tolower(habprop$habitat)
habprop2 = habprop2[match(names(habprop2),names(fimpref)[-1])]

manly.rsr = data.frame(matrix(NA,nrow=nrow(fimpref),ncol=length(habprop2)+1))
names(manly.rsr) = c('group',names(habprop2))
for(i in 1:nrow(fimpref)){
  #i=1
  tmp.u = as.numeric(fimpref[i,-1]);names(tmp.u) = names(fimpref)[-1]
  tmp.rsr = widesI(u=tmp.u,a=habprop2,alpha=0.05)
  manly.rsr[i,1] = fimpref$group[i]
  manly.rsr[i,2:ncol(manly.rsr)] = tmp.rsr$Bi
}

