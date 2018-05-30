library(data.table)
library(raster)
library(rgeos)
library(rgdal)


### OK, to bring all this together at last:
#1) street tree simulation -- "map" is output
#2) Overlay urban forest estimation (edge+interior fractions)
#3) overlay fractional street tree estimation --> Sum(forest+nonforest)
#4) correct final high-biomass "street" pixel fractions to forest.

## street tree results and reconstruction
#### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
## does not collate the genera records, but could (May 28, 2018)
### small biomass pixels first -- have consistent pix id scheme
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.small*")]
npp.dump.chunks <- sub('.*small.', '', npp.dump)

## to get basal area
# ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)
ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)

### for each pixel, get median estimated npp, median # trees
container <- data.frame()
for(c in 1:length(npp.dump)){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "num.trees" object
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  
  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "dbh.stret.small" object
  ba.track <- rep(9999, length(cage.dbh))
  dbh.track <- rep(9999, length(cage.dbh))
  print(paste("patience! Doing BA calculations on chunk", npp.dump.chunks[c]))
  for(b in 1:length(cage.dbh)){
    if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
      ba.track[b] <- NA
      dbh.track[b] <- NA
    } else{
      ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba (=m2) values for each sample, need to convert based on (canopy) area of pixel
      dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
    }
    if(b%%1000==0){print(b)}
  }
  
  ## bind to container
  g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
  container <- rbind(container, g)
  hist(tmp.npp, main=paste(npp.dump.chunks[c]))
  hist(tmp.num, main=paste(npp.dump.chunks[c]))
  hist(dbh.track, main=paste(npp.dump.chunks[c]))
  hist(ba.track, main=paste(npp.dump.chunks[c]))
}

## collect, export, make maps
names(container) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
                      "ba.street.sim", "med.dbh.street.sim")

# write.csv(container, "processed/boston/bos.street.trees.small.npp.simulatorv1.results.csv")
# dim(container) # 98974 -- small only -- compare to 108974 when big trees included
# sum(container$tree.num.street.sim, na.rm=T) ##938k

### tedious reconstruction of the raster comenses here
### this is how the biomass raster was processed prior to the street tree simulator start
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.raw <- biom ### lucky fucker, the 
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,id.proc1:=1:dim(biom.dat)[1]] # 1:354068 --> 
#### !!! you jammy bastard, the groomed data contains the same number of rows as there are cells in the cropped biom raster so the ID's track fine
biom.dat.proc <- biom.dat ## save a copy here

map <- merge(x=biom.dat.proc, y=container, by="id.proc1", all.x=T, all.y=T)
# map[!is.na(bos.biom30m),] ## seems alright
# plot(map[,bos.biom30m], map[,biom.kg]) ## lines up, but doesn't capture anything larger than 20k kg
# range(map[,bos.biom30m], na.rm=T) ## original data still goes up to 50k kg


### big tree processing
### repeat processing for the ~6k "big biomass" pixels
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.big*")]

container.big <- data.frame()
### for each pixel, get median estimated npp, median # trees
for(c in 1){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.big"))
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.big")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.big")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  
  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.big")) ## comes in as "index.track" object
  ba.track <- rep(9999, length(cage.dbh))
  dbh.track <- rep(9999, length(cage.dbh))
  print(paste("patience! Doing BA calculations on chunk big"))
  for(b in 1:length(cage.dbh)){
    if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
      ba.track[b] <- NA
      dbh.track[b] <- NA
    } else{
      ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba values for each sample (in m2, need to compare to canopied area of pixel for proper ba estimation)
      dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
    }
    if(b%%1000==0){print(b)}
  }
  
  ## bind to container
  g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
  container.big <- rbind(container.big, g)
  hist(tmp.npp, main=paste("big"))
  hist(tmp.num, main=paste("big"))
  hist(dbh.track, main=paste("big"))
  hist(ba.track, main=paste("big"))
}

## collect, export, make maps
names(container.big) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
                          "ba.street.sim", "med.dbh.street.sim")

## rebuild container file before final merge into map vectors
container.f <- rbind(container, container.big)
dim(container.f) # 105k cells
map <- merge(x=biom.dat, y=container.f, by="id.proc1", all.x=T, all.y=T)
# map[!is.na(bos.biom30m),] ## seems alright
map[,ba.street.sim:=ba.street.sim/((aoi*bos.can30m)/1E4)] ## convert summed stump BAs to BA m2/ha per canopied area of pixel
# plot(map[,bos.biom30m], map[,biom.kg]) ## lines up, ignores the range above 50k kg
# range(map[,bos.biom30m], na.rm=T) ## original data still goes up to 50k kg
# range(map[, biom.kg], na.rm=T) ## processed only goes up to 30k kg
# map[bos.biom30m>30000, length(bos.biom30m)]/map[,length(bos.biom30m)] ## only 0.5% of pixels greater than 30k kg
write.csv(map, "processed/boston/biom_street/streetsim.v1.results.csv")

# hist(map$ann.npp.street.sim)
# hist(map$tree.num.street.sim)
# hist(map$ba.street.sim) ## peak around 15-20 m2/ha
# hist(map$med.dbh.street.sim) # sharp peak mid 20's, then tiny numbers up to 50cm
# map[,sum(ann.npp.street.sim, na.rm=T)]/(2*(10^3)) ## 12k MgC/yr
# map[,sum(ann.npp.street.sim, na.rm=T)]/(2*(10^3))/(map[,sum(aoi, na.rm=T)]/(10^4)) ## 0.97 MgC/ha/yr, compare Hardiman 1.5 MgC/ha/yr for Boston region
# hist(map[bos.biom30m>0 & is.na(ann.npp.street.sim), bos.biom30m]) ## the ones that failed (3200 total, ~3%) are either very small or above 30k 
# hist(map[bos.biom30m>20000, ann.npp.street.sim]) #450 kgbiomass/yr

# ### export the tifs
# for(j in 1:4){
#   r <- biom
#   r <- setValues(r, map[[(j+5)]])
#   writeRaster(r, paste("processed/boston/", names(map)[5+j], ".v1.tif", sep=""),
#               format="GTiff", overwrite=T)
# }
# 

#### now in-line calculate forest npp edge/interior
library(raster)
library(data.table)
biom <- raster("processed/boston/bos.biom30m.tif")
biom <- crop(biom, aoi)
forest <- raster("processed/boston/bos.forest30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
ed.can <- raster("processed/boston/bos.ed1030m.tif")
can <- raster("processed/boston/bos.can30m.tif")
ed.biom <- raster("processed/boston/bos.biomass.ed10only30m.tif")
forest.biom <- raster("processed/boston/bos.biomass.forestonly30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")

biom.dat <- as.data.table(as.data.frame(biom))
aoi.dat <- as.data.table(as.data.frame(aoi))
forest.dat <- as.data.table(as.data.frame(forest))
ed.can.dat <- as.data.table(as.data.frame(ed.can))
can.dat <- as.data.table(as.data.frame(can))
ed.biom.dat <- as.data.table(as.data.frame(ed.biom))
forest.biom.dat <- as.data.table(as.data.frame(forest.biom))
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat <- cbind(biom.dat, aoi.dat, forest.dat, ed.can.dat, can.dat, ed.biom.dat, forest.biom.dat, isa.dat)
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
names(biom.dat) <- c("biom", "aoi", "forest.area", "ed.can", "can", "ed.biom", "forest.biom", "isa", "pix.ID")

## figure out the relative area of interior forestcanopy in each cell
# biom.dat[, .(range(can,  na.rm=T), range(ed.can, na.rm=T))]
biom.dat[,int.can:=can-ed.can]
# biom.dat[,range(int.can, na.rm=T)] ## includes negative numbers
# biom.dat[int.can<0, length(int.can)] ## 112, some are rounding errors
## kill the obvious bullshit
biom.dat[biom==0, int.can:=0]
biom.dat[biom==0, ed.can:=0]
biom.dat[is.na(biom), ed.can:=NA]
biom.dat[is.na(biom), int.can:=NA]
biom.dat[int.can<0 & int.can>(-0.01), int.can:=0] ## kill rounding errors
# biom.dat[int.can<0, length(int.can)] ## 13
# View(biom.dat[int.can<0,]) ## all areas with 0 or partial forest, places with forest biom>0 are 100% edge biomass
# biom.dat[ed.biom==forest.biom & int.can>0,]
# biom.dat[ed.can>can, ] ## almost all partial pixels, usually forest edge is 100% of forest biomass, seem like minor disagreements in canopy area figuring
biom.dat[ed.can>can, int.can:=0]
# biom.dat[ed.can>can,] ## still 31 records where ed.can doesn't make sense, but fuck it

## now groom biomass
# biom.dat[forest.area==0 & forest.biom>0,] ## 467 where forest biomass is positive with no area
# range(biom.dat[forest.area==0 & forest.biom>0, forest.biom], na.rm=T) ## ~0 to over nine THOUSAND!!!
biom.dat[forest.area==0, forest.biom:=0] ## this eliminates any of the bullshit
biom.dat[forest.area==0, ed.biom:=0] ## also kill bullshit for edge biomass
## also need to eliminate positive forest bioms for invalid out-of-AOI and biom==0
biom.dat[biom==0, forest.biom:=0]
biom.dat[biom==0, ed.biom:=0]
biom.dat[is.na(aoi) | is.na(biom), forest.biom:=0]
biom.dat[is.na(aoi) | is.na(biom), ed.biom:=0]

biom.dat[,int.biom:=forest.biom-ed.biom] # internal forest biomass
# biom.dat[,range(int.biom, na.rm=T)] ## all interior biomass is clear now
biom.dat[,street.biom:=biom-forest.biom] ## street tree biomass
# biom.dat[,range(street.biom, na.rm=T)] ## still have some negatives for street
# View(biom.dat[street.biom<0]) ## all seem to have forest.biom slightly too high
biom.dat[forest.biom>biom, forest.biom:=biom] ## where forest biomass is too high, set to whole-cell biomass
biom.dat[,street.biom:=biom-forest.biom]
# biom.dat[,range(street.biom, na.rm=T)] ## all ok now
# test <- biom.dat[,street.biom+ed.biom+int.biom]
# test2 <- biom.dat[,street.biom+forest.biom]
# plot(biom.dat$biom, test, main="residuals") ## close enough
# plot(biom.dat$biom, test2, main="forest") ## perfect

### How we determined rate of growth based on edge position
andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv") 
andy.bai <- as.data.table(andy.bai)
### set up some basic allometrics to get biomass
## Jenkins, C.J., D.C. Chojnacky, L.S. Heath and R.A. Birdsey. 2003. Forest Sci 49(1):12-35.
## Chojnacky, D.C., L.S. Heath and J.C. Jenkins. 2014. Forestry 87: 129-151.
# Pinus rigida, Pinus strobus, both spg>0.45
# Acer rubrum, Aceraceae <0.50 spg
# Quercus alba, Quercus coccinea, Quercus rubra, Quercus velutina --> deciduous Fagaceae
b0.l <- c(-3.0506, -2.0470, -2.0705)
b1.l <- c(2.6465, 2.3852, 2.4410)
biom.pred <- function(x, b0, b1){exp(b0+(b1*log(x)))} ## dbh in cm, biom. in kg
## artifact: dbh's have hidden NA's that are marked as repeating numbers, basically when the dbh gets too small -- they are higher than the min dbh though
cleanup <- as.matrix(andy.bai[,7:33])
for(r in 1:nrow(cleanup)){
  bust <- which(diff(cleanup[r,], lag = 1)>0) ## tells you where the NA's are located
  if(length(bust)>0){
    cleanup[r, bust:ncol(cleanup)] <- NA
  }
}
cleanup <- as.data.table(cleanup)
andy.bai <- cbind(andy.bai[,1:6], cleanup)
names(andy.bai)[7:33] <- paste0("dbh", 2016:1990)
andy.bai[,incr.ID:=seq(1:dim(andy.bai)[1])]
names(andy.bai)[1:6] <- c("Plot.ID", "Tree.ID", "Spp", "X", "Y", "Can.class")

### biomass change as a % of previous biomass per year from the increment data
taxa <- list(c("PIRI", "PIST"),
             c("ACRU"),
             c("QUAL", "QUCO", "QURU", "QUVE"))
biom.rel <- list()
bop <- data.table()
for(sp in 1:3){
  tmp <- andy.bai[Spp %in% taxa[[sp]],]
  a <- tmp[, lapply(tmp[,7:33], function(x){biom.pred(x, b0.l[sp], b1.l[sp])})]
  a <- cbind(tmp[,.(Tree.ID, incr.ID, Spp, Y, dbh2016)], a)
  bop <- rbind(bop, a)
}
bop$Spp <- as.character(bop$Spp) ### bop is the biomass time series for every core
names(bop)[6:32] <- paste0("biom", 2016:1990)

biom.rel <- data.frame(c(2015:1990))
## get lagged 1yr biomass differences
for(i in 1:dim(bop)[1]){
  te <- unlist(bop[i,6:32])
  te.di <- (-1*diff(te, lag=1))
  biom.rel[,i+1] <- te.di/te[-1]
}
names(biom.rel)[-1] <- paste0("Tree", bop[,incr.ID])
names(biom.rel)[1] <- "Year" ## this has all the relative biomass gains by year (row) for each tree with increment (column)
# plot(biom.rel$Year, biom.rel$Tree2)
# plot(biom.rel$Year, biom.rel$Tree4)
# names(biom.rel)
# for(u in 41:60){  
#   plot(biom.rel$Year, biom.rel[,u], main=paste(colnames(biom.rel)[u]))
# }
## everything looks like a slow decline in relative growth (not clear if that is just a funciton of trees getting bigger or canopy closing)
## most are in the 2-8% range, but some are v. high (20-40%)
mean.na <- function(x){mean(x, na.rm=T)}
m <- apply(biom.rel[,-1], FUN=mean.na, 2)
growth.mean <- data.frame(cbind(c(sub("Tree", '', colnames(biom.rel[,-1]))), m))
names(growth.mean) <- c("incr.ID", "growth.mean")
growth.mean$incr.ID <- as.integer(as.character(growth.mean$incr.ID))
andy.bai <- merge(andy.bai, growth.mean, by="incr.ID")
andy.bai[Y<10, seg:=10]
andy.bai[Y>=10 & Y<20, seg:=20]
andy.bai[Y>=20, seg:=30]
andy.bai[seg==10, seg.F:="A"]
andy.bai[seg==20, seg.F:="B"]
andy.bai[seg==30, seg.F:="C"]
andy.bai$seg.F <- as.factor(andy.bai$seg.F)
andy.bai$growth.mean <- as.numeric(as.character(andy.bai$growth.mean))
andy.bai[,avg.dbh:=apply(as.matrix(andy.bai[,8:34]), FUN=mean.na, 1)]
## note: final interaction model only applies correction coefficients for the interior segments
andy.bai[seg.F=="A", seg.Edge:="E"]
andy.bai[seg.F %in% c("B", "C"), seg.Edge:="I"]
mod.int.edge <- summary(lm(log(growth.mean)~log(avg.dbh)*seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior, no sig. dbh*edge interaction

####
### bring in full andy.dbh data set, groom
ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv")
andy.dbh <- as.data.table(andy.dbh)
names(andy.dbh) <- c("Tree.ID", "Plot.ID", "Spp", "X", "Y", "dbh", "Can.class")
andy.dbh[Y<10, seg:=10]
andy.dbh[Y>=10 & Y<20, seg:=20]
andy.dbh[Y>=20, seg:=30]
andy.dbh[,ba:=ba.pred(dbh)]
setkey(andy.dbh, Tree.ID)
andy.dbh[Spp=="FRAL", Spp:="FRAM"] # Fraxinus americana

### get biomass for all the trees in andy.dbh
biom.pred.key <- data.frame(unique(andy.dbh$Spp))
# For HAVI (Witch hazel), CRSp ?? and TASp ?? we will use Jenkins generalized hardwood coefficients
b0 <- c(-2.0705, -3.0506, -2.6177, -2.0705, -2.0705, -2.2118, -2.0470, -2.0705, -2.6327, -2.3480, -2.48, -2.2271, -1.8384, -2.48, -2.48, -1.8011, -2.2271)
b1 <- c(2.4410, 2.6465, 2.4638, 2.4410, 2.4410, 2.4133, 2.3852, 2.4410, 2.4757, 2.3876, 2.4835, 2.4513, 2.3524, 2.4835, 2.4835, 2.3852, 2.4513)
biom.pred.key <- cbind(biom.pred.key, b0, b1)
biom.pred.key <- as.data.table(biom.pred.key)
names(biom.pred.key) <- c("Spp", "b0", "b1")
for(b in 1:dim(biom.pred.key)[1]){
  andy.dbh[Spp==biom.pred.key[b, Spp], biom:=exp(biom.pred.key[b, b0]+(biom.pred.key[b,b1]*log(dbh)))]
}

### generalized growth coefficients for edge/interior, based on andy.bai
b0 <- mod.int.edge$coefficients[1]
b1 <- mod.int.edge$coefficients[2]
b2 <- mod.int.edge$coefficients[3]

andy.dbh[seg==10, growth:=biom*exp(b0+(b1*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth:=biom*exp(b0+(b1*log(dbh))+b2)]
g <- andy.dbh[, .(sum(growth), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
g[,rel.gain:=V1/V2] ## this deals in forest growth as a function of biomass, not of forest area
plot(g$V2, g$rel.gain, col=as.numeric(g$seg))
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
g[,mean(rel.gain), by=seg.F]





npp.forest.edge <- biom.dat[,ed.biom]*0.0418 ## apply growth of edge trees to edge biomass
npp.forest.int <- biom.dat[,int.biom]*0.0154 ## apply growth of interior trees to interior biomass
range(npp.forest.tot, na.rm=T) ### 0-1093 kg biomass/cell/yr estimate for forest trees based on andy data
(1093/(1000*2)/900)*1E4 ## this is up to 6 MgC/ha/yr for forest fragments in Boston
