library(data.table)
library(raster)
library(rgeos)
library(rgdal)

# setwd("/projectnb/buultra/atrlica/FragEVI/")


##########
### street tree reconstruction of results

vers <- 4 ## which model run are we picking at
#### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = paste("ann.npp.street.v", vers, ".weighted*", sep=""))]
npp.dump.chunks <- sub('.*weighted.', '', npp.dump)
npp.dump.chunks <- sub("\\.sav.*", "", npp.dump.chunks) ## later version get named .sav

## to get basal area
# ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)
ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)

### for each pixel, get median estimated npp, median # trees
container <- data.frame()
med.dbh.rec <- numeric() ## keep a running tally of the dbh of every tree in the nearest-to-median-npp sample in each pixel
dbh.dump <- list() ## a place to put actual dbh samples from targeted retrievals for deeper analysis (e.g. what do the median retrievals look like?)
i=1
for(c in 1:length(npp.dump)){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.num.sims <- sapply(cage.ann.npp, FUN=length) ## number of successful simulations recorded
  tmp.sim.incomp <- rep(0, length(tmp.num.sims))
  tmp.sim.incomp[tmp.num.sims<100] <- 1 ## vector tracking the incomplete simulations
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "num.trees" object
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  load(paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.proc <- proc.track
  load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.biom.sim <- sapply(cage.biom.sim, FUN=median)
  load(paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.attempts <- attempts.track
  load(paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object
  tmp.wts <- sapply(cage.wts, FUN=max)

  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "dbh.stret.small" object
  ba.grand <- rep(9999, length(cage.dbh)) ## BA values for every retreival per cell
  dbh.grand <- rep(9999, length(cage.dbh)) ## dbh values for every retreival per cell
  print(paste("patience! Doing BA calculations on chunk", npp.dump.chunks[c]))
  for(b in 1:length(cage.dbh)){
    if(length(cage.ann.npp[[b]])<5){  ## if too few NPP retreivals were made
      ba.grand[b] <- NA
      dbh.grand[b] <- NA
      dbh.dump[[i]] <- NA
    } else{
      ba.grand[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the total ba (=m2) for each sample, need to convert based on (canopy) area of pixel
      dbh.grand[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
      ### these look at the summary stats for a particular retreival in each cell (ex. the retrieval nearest the median npp)
      dev <- abs(cage.ann.npp[[b]]-tmp.npp[b]) ## deviance of individual NPP retreivals from the median NPP for this cell
      rrr <- which(dev==min(dev)) ## which retreival is the closest to the median
      dbh.dump[[i]] <- cage.dbh[[b]][[rrr[1]]] ## track the tree sample nearest to every npp median in every pixel (just take the first instance, fuck it)
    }
    if(b%%1000==0){print(b)}
    i=i+1 ## keep track of number of pixels processed
  }

  ## bind results in container
  container <- rbind(container, 
                     cbind(tmp.index, tmp.biom, ## basic cell tracking here. each row is 1 pixel
                           tmp.npp, tmp.num, dbh.grand, ba.grand, tmp.biom.sim, ## median npp and tree number for all retrievals
                           tmp.wts, tmp.num.sims, tmp.sim.incomp, tmp.attempts, tmp.proc) ### metrics for how well the simulator performed
  )
}

## collect, export, make maps
names(container) <- c("pix.ID", "biom.kg", 
                      "med.ann.npp.all", "med.tree.num.all", "med.dbh.all", "med.ba.all", "med.biom.all",
                      "max.wts", "num.sims", "sim.incomp", "attempts", "proc.status")

## figure out median dbh, ba, and count for the tree sample in each pixel closest to median npp
pixM.tree.num <- sapply(dbh.dump, FUN=length) ## tree number in median npp retreival per pixel
pixM.dbh <- sapply(dbh.dump, FUN=median) ## the median dbh for the selection of trees nearest the median npp in each cell
pixM.ba <- sapply(sapply(dbh.dump, FUN=ba), FUN=sum) ## summed BA for the selection of trees in the sample nearest median npp 
hist(pixM.tree.num) ## so this is the distribution in tree density per pixel in the most common retreival in each cell
hist(pixM.dbh) ## this is the distribution of MEDIAN dbh for the most common retreival in each cell
hist(pixM.ba) ## this is the distribution of total BA (m2) in each cell (not corrected for canopy) according to the tree sample nearest the median npp
med.dbh.rec <- c(med.dbh.rec, unlist(dbh.dump)) ## append the median dbh record
hist(med.dbh.rec) ## this is the distribution of dbh if you actually went out and counted every simulated tree in the pixels
median(med.dbh.rec, na.rm=T) ### YESSS BITCHESSSSS median is same as street tree records
## could also look at these distributions in e.g. the 25h and 75th percentile retreivals of NPP for each pixel
container <- cbind(container, pixM.ba, pixM.dbh, pixM.tree.num)
sum(duplicated(container$pix.ID)) ## some duplicated pixels, remove
container <- container[!(duplicated(container$pix.ID)),]

## raster reconstruction
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,pix.ID:=1:dim(biom.dat)[1]]
map <- merge(x=biom.dat, y=container, by="pix.ID", all.x=T, all.y=T)

write.csv(map, paste("processed/streettrees.npp.simulator.v", vers, ".results.csv", sep=""))




# ### export the tifs
# for(j in 1:4){
#   r <- biom
#   r <- setValues(r, map[[(j+5)]])
#   writeRaster(r, paste("processed/boston/", names(map)[5+j], ".v1.tif", sep=""),
#               format="GTiff", overwrite=T)
# }
# 


## brief exploratory
### do simulations track cell biomass well enough?
container <- as.data.table(container)
plot(container$biom.kg, container$med.biom.all) ## a few wonkeys but mostly in a tight range -- no successful sims over ~35k
abline(a=0, b=1)
container[,biom.dev:=med.biom.all-biom.kg]
plot(container$biom.kg, container$biom.dev) ## most within +/- 50kg till right at 30k, expands to 100kg

### how prevalent are failures, where are they
sum(container$proc.status==0)/dim(container)[1] # 3.46% failure to simulate
sum(container$biom.kg>30000)/dim(container)[1] # 1.77% are >30k -- so we are getting some but not all of the largest 1%
container[biom.kg>25000 & proc.status==0, length(biom.kg)]/container[biom.kg>25000, length(biom.kg)] ## 66% of pixels >25k failed
container[biom.kg>25000 & proc.status==0, length(biom.kg)]/container[proc.status==0, length(biom.kg)] ## 73% of failures are over 25k
map[proc.status==0, median(bos.can30m, na.rm=T)] ## median 99.8% canopy
hist(map[proc.status==0, (bos.can30m)]) ## most heavily covered except for a small number that must be low biomass
hist(map[proc.status==0, bos.biom30m]) ### bimodal, very low biomass or very high biomass
### might need to just toss anything above ~20k kg -- call it "much more like a forest than a piece of city"

### how is productivity organized in space
plot(container$biom.kg, container$med.ann.npp.all) # tent shape, nearly linear increase up to ~25k then steep decline (sampling large old trees?)
abline(v=22000) ## peaks right at 22k then steep decline
plot(map$bos.can30m, map$med.ann.npp.all) ## everything increases but with inreasing variance up to full canopy where it can be all over
hist(map[bos.biom30m>10, med.ann.npp.all])
hist(map[bos.biom30m>10 & bos.biom30m<22000, med.ann.npp.all]) ## how different are our results if we chuck the weird ass high biomass pixels?

## NPP rate, city wide
map[bos.biom30m>10, npp.MgC.ha:=((med.ann.npp.all/1000*(1/2))/aoi)*1E4]
map[aoi>800, range(npp.MgC.ha, na.rm=T)] ## 0.04-5.2 MgC/ha/yr
hist(map[aoi>800, npp.MgC.ha])
map[bos.biom30m>22000, length(bos.biom30m)]/map[bos.biom30m>10, length(bos.biom30m)] ##5.8% of pix are above 22k
map[is.finite(med.ann.npp.all) & aoi>800 & bos.biom30m>22000, length(med.ann.npp.all)]/map[aoi>800 & bos.biom30m>22000, length(med.ann.npp.all)] ## we got a value for 66% of the large pixels
map[aoi>800, sum(med.ann.npp.all/(2*1000), na.rm=T)] #13161 MgC/yr
map[aoi>800, sum(med.ann.npp.all/(2*1000), na.rm=T)]/(map[, sum(aoi, na.rm=T)]/1E4) ### 1.06 MgC/ha/yr

### excluding large pixels
map[aoi>800 & bos.biom30m<22000, sum(med.ann.npp.all/(2*1000), na.rm=T)] #11789 MgC/yr
map[aoi>800 & bos.biom30m<22000, sum(med.ann.npp.all/(2*1000), na.rm=T)]/(map[, sum(aoi, na.rm=T)]/1E4) ### 0.95 MgC/ha/yr

## average small and large pixel productivity
map[aoi>800 & bos.biom30m<22000, median(npp.MgC.ha, na.rm=T)] ## 0.94 MgC/ha/yr for small pixels
map[aoi>800 & bos.biom30m>22000, median(npp.MgC.ha, na.rm=T)] ## ### 3.65 MgC/ha/yr
 ## contrast to Andy-forest-based estimate, #13.8k tC, 1.1 tC/ha/yr mean estimate (range was ~7k-22k)







##### OLD method for reconstituting street tree simulator results from stored object files
# ### OK, to bring all this together at last:
# #1) street tree simulation -- "map" is output
# #2) Overlay urban forest estimation (edge+interior fractions)
# #3) overlay fractional street tree estimation --> Sum(forest+nonforest)
# #4) correct final high-biomass "street" pixel fractions to forest.
# 
# ## street tree results and reconstruction
# #### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
# ## does not collate the genera records, but could (May 28, 2018)
# ### small biomass pixels first -- have consistent pix id scheme
# obj.dump <- list.files("processed/boston/biom_street/")
# npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.small*")]
# npp.dump.chunks <- sub('.*small.', '', npp.dump)
# 
# ## to get basal area
# # ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)
# ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)
# 
# ### for each pixel, get median estimated npp, median # trees
# container <- data.frame()
# for(c in 1:length(npp.dump)){
#   load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
#   tmp.npp <- sapply(cage.ann.npp, FUN=median)
#   tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
#   
#   ## grab the corresponding other dump files
#   load(paste("processed/boston/biom_street/num.trees.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "num.trees" object
#   tmp.num <- sapply(cage.num.trees, FUN=median)
#   tmp.num[tmp.num==9999] <- NA ## filter the NA flags
#   load(paste("processed/boston/biom_street/index.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "index.track" object
#   tmp.index <- index.track
#   load(paste("processed/boston/biom_street/biom.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "biom.track" object
#   tmp.biom <- biom.track
#   
#   ### figure median dbh and median basal area for each pixel
#   load(paste("processed/boston/biom_street/dbh.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "dbh.stret.small" object
#   ba.track <- rep(9999, length(cage.dbh))
#   dbh.track <- rep(9999, length(cage.dbh))
#   print(paste("patience! Doing BA calculations on chunk", npp.dump.chunks[c]))
#   for(b in 1:length(cage.dbh)){
#     if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
#       ba.track[b] <- NA
#       dbh.track[b] <- NA
#     } else{
#       ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba (=m2) values for each sample, need to convert based on (canopy) area of pixel
#       dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
#     }
#     if(b%%1000==0){print(b)}
#   }
#   
#   ## bind to container
#   g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
#   container <- rbind(container, g)
#   hist(tmp.npp, main=paste(npp.dump.chunks[c]))
#   hist(tmp.num, main=paste(npp.dump.chunks[c]))
#   hist(dbh.track, main=paste(npp.dump.chunks[c]))
#   hist(ba.track, main=paste(npp.dump.chunks[c]))
# }
# 
# ## collect, export, make maps
# names(container) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
#                       "ba.street.sim", "med.dbh.street.sim")
# 
# # write.csv(container, "processed/boston/bos.street.trees.small.npp.simulatorv1.results.csv")
# # dim(container) # 98974 -- small only -- compare to 108974 when big trees included
# # sum(container$tree.num.street.sim, na.rm=T) ##938k
# 
# ### tedious reconstruction of the raster comenses here
# ### this is how the biomass raster was processed prior to the street tree simulator start
# biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
# aoi <- raster("processed/boston/bos.aoi30m.tif")
# biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
# biom.raw <- biom ### lucky fucker, the 
# biom.dat <- as.data.table(as.data.frame(biom))
# biom.dat[,aoi:=as.vector(getValues(aoi))]
# biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
# can <- raster("processed/boston/bos.can30m.tif")
# can.dat <- as.data.table(as.data.frame(can))
# biom.dat <- cbind(biom.dat, can.dat)
# biom.dat[,id.proc1:=1:dim(biom.dat)[1]] # 1:354068 --> 
# #### !!! you jammy bastard, the groomed data contains the same number of rows as there are cells in the cropped biom raster so the ID's track fine
# biom.dat.proc <- biom.dat ## save a copy here
# 
# map <- merge(x=biom.dat.proc, y=container, by="id.proc1", all.x=T, all.y=T)
# 
# ### big tree processing
# ### repeat processing for the ~6k "big biomass" pixels
# obj.dump <- list.files("processed/boston/biom_street/")
# npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.big*")]
# 
# container.big <- data.frame()
# ### for each pixel, get median estimated npp, median # trees
# for(c in 1){
#   load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
#   tmp.npp <- sapply(cage.ann.npp, FUN=median)
#   tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
#   
#   ## grab the corresponding other dump files
#   load(paste("processed/boston/biom_street/num.trees.street.big"))
#   tmp.num <- sapply(cage.num.trees, FUN=median)
#   tmp.num[tmp.num==9999] <- NA ## filter the NA flags
#   load(paste("processed/boston/biom_street/index.track.street.big")) ## comes in as "index.track" object
#   tmp.index <- index.track
#   load(paste("processed/boston/biom_street/biom.track.street.big")) ## comes in as "biom.track" object
#   tmp.biom <- biom.track
#   
#   ### figure median dbh and median basal area for each pixel
#   load(paste("processed/boston/biom_street/dbh.street.big")) ## comes in as "index.track" object
#   ba.track <- rep(9999, length(cage.dbh))
#   dbh.track <- rep(9999, length(cage.dbh))
#   print(paste("patience! Doing BA calculations on chunk big"))
#   for(b in 1:length(cage.dbh)){
#     if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
#       ba.track[b] <- NA
#       dbh.track[b] <- NA
#     } else{
#       ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba values for each sample (in m2, need to compare to canopied area of pixel for proper ba estimation)
#       dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
#     }
#     if(b%%1000==0){print(b)}
#   }
#   
#   ## bind to container
#   g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
#   container.big <- rbind(container.big, g)
#   hist(tmp.npp, main=paste("big"))
#   hist(tmp.num, main=paste("big"))
#   hist(dbh.track, main=paste("big"))
#   hist(ba.track, main=paste("big"))
# }
# 
# ## collect, export, make maps
# names(container.big) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
#                           "ba.street.sim", "med.dbh.street.sim")
# 
# ## rebuild container file before final merge into map vectors
# container.f <- rbind(container, container.big)
# dim(container.f) # 105k cells
# map <- merge(x=biom.dat, y=container.f, by="id.proc1", all.x=T, all.y=T)
# # map[,ba.street.sim:=ba.street.sim/((aoi*bos.can30m)/1E4)] ## convert summed stump BAs to BA m2/ha per canopied area of pixel
# ## don't convert BA by area -- still deciding what area to use
# 
# names(map) <- c("id.proc1", "biom.kg", "aoi", "can", "biom.track", "ann.npp.street.sim", 
#                 "tree.num.street.sim", "ba.street.sim", "med.dbh.street.sim")
# write.csv(map, "processed/boston/biom_street/streetsim.v1.results.csv")
# 
# hist(map$ann.npp.street.sim)
# hist(map$tree.num.street.sim)
# hist(map$ba.street.sim) 
# hist(map$med.dbh.street.sim) # sharp peak mid 20's, then tiny numbers up to 50cm
# map[,sum(ann.npp.street.sim, na.rm=T)]/(2*(10^3)) ## 12k MgC/yr
# map[,sum(ann.npp.street.sim, na.rm=T)]/(2*(10^3))/(map[,sum(aoi, na.rm=T)]/(10^4)) ## 0.97 MgC/ha/yr, compare Hardiman 1.5 MgC/ha/yr for Boston region
# hist(map[bos.biom30m>0 & is.na(ann.npp.street.sim), bos.biom30m]) ## the ones that failed (3200 total, ~3%) are either very small or above 30k
# hist(map[bos.biom30m>20000, ann.npp.street.sim]) #450 kgbiomass/yr in dense cells
# 
# # ### export the tifs
# # for(j in 1:4){
# #   r <- biom
# #   r <- setValues(r, map[[(j+5)]])
# #   writeRaster(r, paste("processed/boston/", names(map)[5+j], ".v1.tif", sep=""),
# #               format="GTiff", overwrite=T)
# # }
# # 

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

### where we get the generalized biomass gain per kg in edge vs. interior forest
andy.dbh[seg==10, growth:=biom*exp(b0+(b1*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth:=biom*exp(b0+(b1*log(dbh))+b2)]
g <- andy.dbh[, .(sum(growth), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
g[,rel.gain:=V1/V2] ## this deals in forest growth as a function of biomass, not of forest area
plot(g$V2, g$rel.gain, col=as.numeric(g$seg))
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
g[,mean(rel.gain), by=seg.F] ## 4.18% edge vs. 1.54% interior, annual gain on biomass per edge category

npp.forest.edge <- biom.dat[,ed.biom]*0.0418 ## apply growth of edge trees to edge biomass
npp.forest.int <- biom.dat[,int.biom]*0.0154 ## apply growth of interior trees to interior biomass
# range(npp.forest.tot, na.rm=T) ### 0-1093 kg biomass/cell/yr estimate for forest trees based on andy data
# (1093/(1000*2)/900)*1E4 ## this is up to 6 MgC/ha/yr for forest fragments in Boston
npp.forest <- cbind(biom.dat$pix.ID, npp.forest.edge, npp.forest.int)
colnames(npp.forest) <- c("pix.ID", "npp.forest.edge", "npp.forest.int")


### now apply npp in the street biomass in the mixed pixels
mix <- biom.dat[forest.area>0.05 & forest.area<0.95,] ## 8783 pix are mixed, should reestimate street tree contribution here
# range(mix[, street.biom], na.rm=T) #- to 33k
# biom.dat[forest.area>0, length(forest.biom)] ## contrast 18k pixels have some forest area
plot(mix$forest.biom, mix$street.biom)
## ok ok, so this looks like the street tree biomass portion of each pixel, re=simulate

# mix[street.biom<20000,] ## <100 truly massive street.biom pixels
# mix[street.biom<5,] ## 21 have itty bitty biomass that can't be fitted with the smallest street trees
mix <- mix[street.biom<20000,] ###  restrict to the fittable street fractions (we will treat v. large street.biom as "forest")
mix <- mix[street.biom>5,] ## eliminate street biomasses too small to fit with street trees
## 8094 total mixed pixels now


### loop the mixed pixels back through to get street tree npp (similar process to street tree script, takes ~8hrs)
library(raster)
library(rgdal)
library(rgeos)

### biomass equation biom~dbh
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## get the street tree record prepped
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
# a <- street[,length(ann.npp)/dim(street)[1], by=genus]
# a <- a[order(a$V1, decreasing = T),]
# # sum(a$V1[1:12]) ### top 12 is 94% of trees present
# focus <- a$genus[1:12]
street[, delta.dbh:=dbh.2014-dbh.2006]
street[, record.good:=0]
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] # exclude fishy looking records
street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
street[record.good==1, at.npp.ann:=at.delta.biom/8]
street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
# street[record.good==1,range(dbh.2006, na.rm=T)]

### biomass increment model based on dbh (exponential function)
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far
# plot(street[record.good==1, log(at.biom.2006)], street[record.good==1, log(at.npp.ann.rel)])

## prep street tree data and biomass data for processing (small biomass first)
clean <- street[record.good==1 & dbh.2006>=5,] # get a good street tree set ready
setkey(clean, at.biom.2006)

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
cage.genus <- list()
index.track <- rep(9999, dim(mix)[1])
biom.track <- rep(9999, dim(mix)[1])

## loop each row (pixel) of the chunk
for(t in 1:dim(mix)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  cage.genus[[t]] <- list()
  cage.dbh[[t]] <- list()
  x <- 0
  q <- 0
  while(x<100 & q<2000){ ## select x workable samples, or quit after q attempts
    ## grab a random number of randomly selected trees out of the street tree database
    grasp <- clean[at.biom.2006<(1.10*mix[t, street.biom]), .(dbh.2006, at.biom.2006, genus)] # don't select any trees that are bigger biomass than the pixel total
    n <- min(c(80, dim(grasp)[1])) ## if pixel biomass is tiny, don't try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## cummulative tally of biomass
    d=1 # track of the number of trees
    while(w<(0.9*mix[t, street.biom])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### check if you've got too much biomass or if you've packed them in too tight
    ## this checks if the total BA (m2/ha) of the street trees is lower than 40 inside the non-forest
    ### Question: what is the right area to evaluate BA for? Is it m2/ha pervious? per ha canopy? per ha ground?
    if((grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]/(mix[t, 1-forest.area]*900)*1E4)<40 & w<(1.10*mix[t, street.biom])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
      cage.dbh[[t]][[x]] <- grasp[1:d, dbh.2006] ## which dbhs did you select
      cage.genus[[t]][[x]] <- grasp[1:d, genus] # running list of which genera you select
      #       print(paste("recorded", x))
    }
    q=q+1 ## record this as an attempt
    #     print(paste("attempted", q))
  }
  if(q<2000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    biom.track[t] <- mix[t, street.biom]
    index.track[t] <- mix[t, pix.ID]
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    biom.track[t] <- mix[t, street.biom]
    index.track[t] <- mix[t, pix.ID]
    print(paste("pixel", t, "error"))
  }
}

## when complete dump everything back into the save file
save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.mix"))
save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.mix"))
save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.mix"))
save(cage.genus, file=paste("processed/boston/biom_street/genus.street.mix"))
save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.mix"))
save(index.track, file=paste("processed/boston/biom_street/index.track.street.mix"))


### process the mixed pixels back into the map database
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.mix*")]

ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)

container.mix <- data.frame()
### for each pixel, get median estimated npp, median # trees
for(c in 1){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.mix"))
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.mix")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.mix")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  
  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.mix")) ## comes in as "index.track" object
  ba.track <- rep(9999, length(cage.dbh))
  dbh.track <- rep(9999, length(cage.dbh))
  print(paste("patience! Doing BA calculations on chunk mix"))
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
  container.mix <- rbind(container.mix, g)
  hist(tmp.npp, main=paste("mix"))
  hist(tmp.num, main=paste("mix"))
  hist(dbh.track, main=paste("mix"))
  hist(ba.track, main=paste("mix"))
}

## collect, export, make maps
names(container.mix) <- c("id.proc1", "biom.kg", "ann.npp.street.mix", "tree.num.street.mix", 
                          "ba.street.mix", "med.dbh.street.mix")

## merge into map vectors
dim(container.mix) # 8094 cells
dim(container.mix[is.na(container.mix$ann.npp.street.mix),]) ## 281 recalcitrants
map <- merge(x=map, y=container.mix, by="id.proc1", all.x=T, all.y=T)
write.csv(map, "processed/boston/biom_street/streetsim.v2MIX.results.csv")

### final collation, grooming, export
gorgor <- merge(x=biom.dat, y=map, by.x="pix.ID", by.y="id.proc1")
gorgor <- gorgor[order(pix.ID),] ## get ready to dump the npp from forest
## first thing: set any 0 biomass pixels to 0, set all things to street sim v1
gorgor[bos.biom30m==0, total.npp:=0]
gorgor[bos.biom30m>0, total.npp:=ann.npp.street.sim]

# r <- raster(biom)
# r <- setValues(r, gorgor$total.npp)
# writeRaster(r, "processed/boston/bos.npp.street.sim.v1.tif", format="GTiff", overwrite=T)
## OK: the stuff that is still missing a value is either deep forest, high biomass, or random tiny biom non-forest

gorgor <- merge(x=gorgor, y=npp.forest, by="pix.ID")
gorgor[!is.na(bos.biom30m), growth.type:=1] ## for street trees
gorgor[forest.area>0, total.npp:=npp.forest.edge+npp.forest.int]
gorgor[forest.area>0, growth.type:=2] ## for forested pixels

## don't NA out anything until you've gotten a forest-only npp
gorgor[forest.area<=0, npp.forest.edge:=NA]
gorgor[forest.area<=0, npp.forest.int:=NA]
gorgor[is.na(forest.area), npp.forest.edge:=NA]
gorgor[is.na(forest.area), npp.forest.int:=NA]
gorgor[,range(npp.forest.edge, na.rm=T)]
gorgor[,range(npp.forest.int, na.rm=T)]
gorgor[forest.area>=0,range(total.npp, na.rm=T)]

r <- raster(biom)
r <- setValues(r, gorgor$npp.forest.edge+npp.forest.int)
writeRaster(r, "processed/boston/bos.npp.forest.edge+int.tif", format="GTiff", overwrite=T)
### this is indeed just the npp for the forested parts, but predicts low npp along edges because it is missing some biomass from the partial street trees 
### the street tree sim parts that fail are usually over the forest parts that show higher npp

# dim(gorgor) # 354k for total raster
# gorgor[!is.na(tree.num.street.mix),] #7813 valid mixed pixels
gorgor[!is.na(total.npp), length(total.npp)] #nearly 137k npp pixels estimated
gorgor[!is.na(bos.biom30m), length(bos.biom30m)] ## almost 138k -- so about 1k pix have biomass but neither forest nor street estimate
hist(gorgor[is.na(total.npp), bos.biom30m]) ### these remaining 1000 pix are either low biomass that didn't sim properly or v. high biomass and (often?) forest adjacent


#### data hygiene
# plot(gorgor$pix.ID)
# gorgor[,.(range(bos.biom30m, na.rm=T), range(biom, na.rm=T), range(biom.kg.x, na.rm=T), range(biom.kg.y, na.rm=T))] ## makes sense: total range x2, range of street sim, range of mixed
# plot(gorgor$aoi.x, gorgor$aoi.y) # same
# plot(gorgor$biom, gorgor$forest.biom+gorgor$street.biom) ## same
# gorgor[!is.na(biom), length(biom)] # 141288 have valid biomass
# gorgor[!is.na(biom.kg.x), length(biom.kg.x)] # 105105 from small street sim
# gorgor[!is.na(biom.kg.y), length(biom.kg.y)] # 8094 mix street trees
# gorgor[!is.na(bos.biom30m), length(bos.biom30m)] # 137963 original biomass records
# (gorgor[!is.na(biom) & is.na(bos.biom30m), range(aoi.x, na.rm=T)]) ## these ~3k are partials (aoi<800) where the biomass was cancelled for the street tree sim

## count high biomass street as edge forest
## count low biomass street as fractional standard median street trees

## we will first count any pixel between 2-98% forest as mixed
gorgor[forest.area<=0.98 & forest.area>=0.02 & street.biom>0, growth.type:=3] #9149 mixed pixels
gorgor[growth.type==2, length(growth.type)] ## contrast 8646 remaining pure forest >98%
gorgor[forest.area<=0.98 & forest.area>=0.02, .(range(street.biom, na.rm=T), range(ann.npp.street.mix, na.rm=T))] ## 0 to 32k

## mixed pixels with no appreciable street biomass anyway (total.npp is captured by forest npp)
View(gorgor[forest.area>=0.02 & is.na(street.biom),]) ##34 of these, no valid biomass
View(gorgor[forest.area>=0.02 & street.biom==0,]) ## 7273 these are entirely forest biomass

## mixed pixels that are >98% forest --> put street biomass over to edge forest npp
View(gorgor[forest.area>0.98 & street.biom>0, ]) ## 348 pixels
gorgor[forest.area>0.98, range(street.biom, na.rm=T)] # 0 867 kg
hist(gorgor[forest.area>0.98 & street.biom>0, street.biom]) ## concentrated below 200kg
gorgor[forest.area>0.98 & street.biom>0, total.npp:=total.npp+(street.biom*0.0418)] ## treat this like edge canopy

## mixed pixels below 2% forest --> npp is street.sim + edge.forest
gorgor[forest.area<0.02 & forest.area>0 & forest.biom>0,] # 740 pixels
gorgor[forest.area<0.02 & forest.area>0 & forest.biom>0, range(forest.biom)] #0-790
gorgor[forest.area<0.02 & forest.area>0 & forest.biom>0, range(bos.biom30m, na.rm=T)]
gorgor[forest.area<0.02 & forest.area>0 & forest.biom>0, total.npp:=ann.npp.street.sim+(forest.biom*0.0418)]

#### mixed pixels between 2 and 98% forest
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.98 & !is.na(ann.npp.street.mix), ] ## we have 7813 with good retreivals on the street trees
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.98 & !is.na(ann.npp.street.mix), total.npp:=total.npp+ann.npp.street.mix]
dude <- gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.98 & is.na(ann.npp.street.mix), ]
### 2009 pixels that are mixed but no street mix retrieval, a bunch have no street biomass
# dude[street.biom>0, unique(growth.type)]
plot(dude$forest.area, dude$street.biom)
#1) below ~5% forest is a large range of up to 30k street biomass 
#2) above ~80 forest is a low range of <1k street biomass and a small handfull >5k
#4) Solution: Treat anything above 5% forest and with positive street tree biomass as edge forest

#####
# for pixels with <5k kg residual street biomass, add npp of standard street trees tot total
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
# a <- street[,length(ann.npp)/dim(street)[1], by=genus]
# a <- a[order(a$V1, decreasing = T),]
# # sum(a$V1[1:12]) ### top 12 is 94% of trees present
# focus <- a$genus[1:12]
street[, delta.dbh:=dbh.2014-dbh.2006]
street[, record.good:=0]
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] # exclude fishy looking records
street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
street[record.good==1, at.npp.ann:=at.delta.biom/8]
street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
street[record.good==1,range(dbh.2006, na.rm=T)]
street[dbh.2006<5, record.good:=0] ## filter for tiny trees, full record is now 2455
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far
#####

std.tree <- biom.pred(street[record.good==1, median(dbh.2006)]) ## about 200kg
std.growth.rate <- exp((mod.biom.rel$coefficients[2]*log(street[record.good==1, median(dbh.2006)]))+mod.biom.rel$coefficients[1]) ## about 8%/yr

## for mixed pixels with <5000 kg residual street biomass & no ann npp retrieval, do standard street trees
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.90 & is.na(ann.npp.street.mix) & street.biom>0 & street.biom<5000,] ## 848 with positive street biomass <5000kg
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.90 & is.na(ann.npp.street.mix) & street.biom>0 & street.biom<5000, total.npp:=total.npp+(street.biom*std.growth.rate)] 

### for mixed pixels >5000kg residual street biomass & no ann.npp retrieval
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.98 & is.na(ann.npp.street.mix) & street.biom>0 & street.biom>=5000,] ## 488 with positive street biomass >=5000kg
gorgor[forest.area>=0.02 & forest.biom>=0 & forest.area<=0.98 & is.na(ann.npp.street.mix) & street.biom>0 & street.biom>=5000, total.npp:=total.npp+(street.biom*0.0418)] ## count the street biomass as forest edge
## for the heavily forested but below 5000kg street biom
gorgor[forest.area>=0.80 & forest.biom>=0 & forest.area<=0.98 & is.na(ann.npp.street.mix) & street.biom>0, total.npp:=total.npp+(street.biom*0.0418)] ## treat these trees like edge trees

### how did we do?
gorgor[bos.biom30m>0 & aoi.x>800, length(bos.biom30m)] #107142 pix with valid biomass
gorgor[bos.biom30m>0 & aoi.x>800 & !is.na(total.npp), length(total.npp)] ## 105791 pix with a valid npp estimate
gorgor[is.finite(bos.biom30m) & aoi.x>800,] ## 135705 complete pix with recorded biomass
map[!is.na(ann.npp.street.sim)] #103998 had a valid street npp sim -- we gained ~1700 unretrieved pixels by simulating forest separeately
# 
# r <- raster(biom)
# r <- setValues(biom, gorgor$total.npp)
# writeRaster(r, filename="processed/boston/bos.npp.forest+sim.tif", format="GTiff", overwrite=T)
# r <- raster(biom)
# r <- setValues(biom, gorgor$growth.type)
# writeRaster(r, filename="processed/boston/bos.growth.type.tif", format="GTiff", overwrite=T)
# 
## looks good, no obvious disjunctions where forest meets street, some linear features of edge enhancement evident
## the remaining gaps in npp seem to be very low or very high/non-forest
hist(gorgor[is.na(total.npp), bos.biom30m])
gorgor[is.na(total.npp), range(forest.area, na.rm=T)] ## all <2% forest
## use a standard street tree for these muffs
gorgor[is.na(total.npp) & is.finite(bos.biom30m) & aoi.x>800, total.npp:=bos.biom30m*std.growth.rate]
## and just to reiterate, if biomass <5, cancel out the npp
gorgor[bos.biom30m<5 & aoi.x>800 & is.na(total.npp), total.npp:=0]
gorgor[bos.biom30m>0 & aoi.x>800 & !is.na(total.npp), length(total.npp)] ## 107142 pix with a valid npp estimate
gorgor[is.finite(bos.biom30m) & aoi.x>800, length(bos.biom30m)] ## 135705 valid biomass pixels still

r <- raster(biom)
r <- setValues(biom, gorgor$total.npp)
writeRaster(r, filename="processed/boston/bos.npp.forest+simV2.tif", format="GTiff", overwrite=T)

(gorgor[,sum(total.npp, na.rm=T)]/(2*1000))/(gorgor[aoi.x>800, sum(aoi.x)]/1E4) ### about 1 Mgc/ha
(gorgor[forest.area>0.8,sum(total.npp, na.rm=T)]/(2*1000))/(gorgor[forest.area>0.8 & aoi.x>800, sum(aoi.x)]/1E4) ### about 2.5 MgC/ha in forest areas
gorgor[forest.area>0.8 & aoi.x>800, mean(can, na.rm=T)] ## 91% canopy avg in forest
gorgor[forest.area<0.2 & aoi.x>800, mean(can, na.rm=T)] ## 26% canopy avg in non-forest
gorgor[forest.area>0.8 & aoi.x>800, mean(bos.biom30m, na.rm=T)] ## 21k kg biomass avg in forest
gorgor[forest.area<0.2 & aoi.x>800, mean(bos.biom30m, na.rm=T)] ## 4k kg biomass avg in non-forest
gorgor[forest.area>0.8 & aoi.x>800, mean(ed.biom, na.rm=T)] ## 21k kg biomass avg in forest

### another big question: What rate do we predict andy plots growing at, what is their biomass density compared to forest and compared to our forest fragments?

plot(gorgor$bos.biom30m, gorgor$total.npp, col=gorgor$growth.type) ## well that's fucking weird and unnatural
### some up to 20000 don't grow at all, wtf. Some grow on a high line, some grow on a low line
plot(gorgor[growth.type==1, bos.biom30m], gorgor[growth.type==1, total.npp], col=gorgor$growth.type) ## street trees, two clouds (<20k, >20k), and a much higher line
View(gorgor[growth.type==1 & total.npp>1000]) ## these all have interior canopy (but NOT biomass), and no forest area
hist(gorgor[growth.type==1 & total.npp>1000, int.can]) ## pretty even, fair amount of low
hist(gorgor[growth.type==1 & total.npp>1000, can]) ## most totally canopied
hist(gorgor[growth.type==1 & total.npp>1000, bos.biom30m]) ## all high, most 30-35k 
## so this population is pixels that are getting classed as street trees but are really more like forests that are getting missed
r <- raster(biom)
gorgor[, bad.street:=0]
gorgor[growth.type==1 & total.npp>1000, bad.street:=1]
r <- setValues(r, gorgor$bad.street)
writeRaster(r, "processed/boston/bos.badstreet.tif", format="GTiff") ## welp no, these are just very dense patches of non-forest canopy
View(gorgor[growth.type==1 & total.npp>1000,])
### on a hunch
gorgor[growth.type==1 & total.npp>1000, total.npp/bos.biom30m] ## YEP these are the high-biomass pixels that I couldn't properly sim so I plugged in the median street tree values for them and bounced

## the forest pixels
plot(gorgor[growth.type==2, bos.biom30m], gorgor[growth.type==2, total.npp]) ## forest, two clouds (<20k, >20k), and a much higher line
View(gorgor[growth.type==2,]) #ok...
## why are some pixels 0 npp?
View(gorgor[growth.type==2 & total.npp<10 & bos.biom30m>10,]) ## these are all either low-biomass forest pixels or have large street biomass that is not represented in the npp calc for unknown reason
## why are some pixels insanely high?
View(gorgor[growth.type==2 & total.npp>1500,]) ## these are mostly non-forest ## something got fucked, a bunch of <2% forest is being called forest
gorgor[growth.type==2 & total.npp>1500, total.npp/street.biom] ## Looks like these are street trees with the "standard tree" plugged in

plot(gorgor[growth.type==3, bos.biom30m], gorgor[growth.type==3, total.npp]) # a spay with a weird population with higher npp per unit biomass


###### Final assessment of street tree+forest simulator v1:
### There are discontinuities between the simulated NPP for "small" and "big" biomass cells -- different distributions selected from
### there are discontinuities between street and forest
### The "mixed" forest/street pixels contain a fair number of artifacts



######
###### Approach 4: Treat all canopy as Andy-like forest
library(raster)
library(data.table)

# ### procecssing 1m edge biomass to 30m cells
# ## identify 1m biomass pixels that are edge canopy (all LULC types)
# ed <- raster("processed/boston/bos.ed10.tif")
# biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# aoi <- raster("processed/boston/bos.aoi.tif")
# # biom <- projectRaster(biom, aoi)
# biom <- crop(biom, aoi)
# 
# ### identify biomass that is edge or not (all LULC included)
# edges.biom <- function(x, y, filename) { # x is edge class, y is biomass, z is forest extent
#   out <- raster(y)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     e <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
#     b <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## biomass
#     b[e==0] <- 0 ## cancel non-edges
#     out <- writeValues(out, b, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# ed.biom <- edges.biom(ed, biom, filename="processed/boston/bos.biomass.ed10only_allLULC.tif")
# r <- raster("processed/boston/bos.biomass.ed10only_allLULC.tif")
# plot(r)

### convert to 1m edge biomass map to 30m aggregate (arcpy)

### now working in the 30m space
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
ed.can <- raster("processed/boston/bos.ed1030m.tif")
can <- raster("processed/boston/bos.can30m.tif")
ed.biom <- raster("processed/boston/bos.biomass.ed10only_allLULC30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")

biom.dat <- as.data.table(as.data.frame(biom))
aoi.dat <- as.data.table(as.data.frame(aoi))
ed.can.dat <- as.data.table(as.data.frame(ed.can))
can.dat <- as.data.table(as.data.frame(can))
ed.biom.dat <- as.data.table(as.data.frame(ed.biom))
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat <- cbind(biom.dat, aoi.dat, ed.can.dat, can.dat, ed.biom.dat, isa.dat)
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
names(biom.dat) <- c("biom", "aoi", "ed.can", "can", "ed.biom", "isa", "pix.ID")

## figure out the relative area of interior forestcanopy in each cell
biom.dat[,int.can:=can-ed.can]
# biom.dat[,range(int.can, na.rm=T)] ## includes negative numbers
# biom.dat[int.can<0, length(int.can)] ## 112, some are rounding errors
## kill the obvious bullshit
biom.dat[biom==0, int.can:=0]
biom.dat[biom==0, ed.can:=0]
biom.dat[is.na(biom), ed.can:=NA]
biom.dat[is.na(biom), int.can:=NA]
biom.dat[int.can<0 & int.can>(-0.01), int.can:=0] ## kill rounding errors
biom.dat[,int.can:=can-ed.can]
# biom.dat[int.can<0, length(int.can)] ## 13
# View(biom.dat[int.can<0,]) ## all areas with 0 or partial forest, places with forest biom>0 are 100% edge biomass
# biom.dat[ed.biom==forest.biom & int.can>0,]
# biom.dat[ed.can>can, ] ## almost all partial pixels, usually forest edge is 100% of forest biomass, seem like minor disagreements in canopy area figuring
biom.dat[ed.can>can, int.can:=0]
# biom.dat[ed.can>can,] ## still 31 records where ed.can doesn't make sense, but fuck it
# biom.dat[ed.can==can,] ## 113k records where canopy is all edge canopy (vast bulk of the pixels)
# biom.dat[aoi>500,] ##138k mostly complete cells

## now ID biomass by edge vs int
biom.dat[,int.biom:=biom-ed.biom] # internal forest biomass
# biom.dat[,range(int.biom, na.rm=T)] ## 0-51k
# biom.dat[,range(ed.biom, na.rm=T)] ## 0-31k (so the peak tips are deep forest somewhere)
# hist(biom.dat[aoi>800, int.biom])
# hist(biom.dat[aoi>800, ed.biom])

### How we determined rate of growth based on edge position
andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv") 
andy.bai <- as.data.table(andy.bai)
### the basic allometrics to get biomass
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
plot(biom.rel$Year, biom.rel$Tree2)
plot(biom.rel$Year, biom.rel$Tree4)
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

### initial edge vs. int growth model was full interactive, but slope modifier term was not significant
# mod.int.edge <- summary(lm(log(growth.mean)~log(avg.dbh)*seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior, no sig. dbh*edge interaction
# b0 <- mod.int.edge$coefficients[1]
# b1 <- mod.int.edge$coefficients[2]
# b2 <- mod.int.edge$coefficients[3]

## removing the slope modifier (i.e. just a different intercept for edge vs. interior)
mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]
### generalized growth coefficients for edge/interior, based on andy.bai
plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean), col=as.numeric(as.factor(andy.bai$seg.Edge)))
abline(b0.bai, b1.bai, col="black")
abline(b0.bai+b2.bai, b1.bai, col="red")

####
### bring in full andy.dbh data set, groom to determine growth rates on area/edge basis
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

### get biomass (kg) as function of dbh (cm) for all the trees in andy.dbh
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
# hist(andy.dbh[biom<2000, biom]) ## a couple of monsters

### where we get the generalized biomass gain per kg in edge vs. interior forest
andy.dbh[seg==10, growth:=biom*exp(b0.bai+(b1.bai*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth:=biom*exp(b0.bai+(b1.bai*log(dbh))+b2.bai)]
g <- andy.dbh[, .(sum(growth), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
g[,rel.gain:=V1/V2] ## this deals in forest growth as a function of biomass, not of forest area
plot(g$V2, g$rel.gain, col=as.numeric(g$seg), xlab="total plot biomass")
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
g[,mean(rel.gain), by=seg.F] ## 4.34% edge vs. 2.85% interior, annual gain on biomass per edge category

biom.dat[,npp.edge:=ed.biom*0.043] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int:=int.biom*0.028] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot:=npp.edge+npp.int]
biom.dat[,range(npp.tot, na.rm=T)] ## 0 1443 kg biomass/cell/yr
biom.dat[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]

biom.dat[,sum(npp.tot, na.rm=T)]*1E-3*(1/2) #13.8k tC ## IDENTICAL TO FIA ESTIMATE?E??E?E?E
(biom.dat[,sum(npp.tot, na.rm=T)]*1E-3*(1/2))/(biom.dat[,sum(aoi, na.rm=T)]*1E-4) ## ie 1.1 tC/ha/yr

hist(biom.dat[aoi>800 & biom>10, npp.tot])
hist(biom.dat[aoi>800 & biom>10, npp.tot.MgC.ha])
## so that's interesting. A lot produce not much, but there's a longer tail of high productivity that starts to look more like real forest -- up to ~8 MgC/ha/yr
## an ancillary question: why does FIA have such a pessimistic idea of productivity, maxes out at 2.5 MgC/ha/yr

## maybe look at how this sorts re. lulc
lulc <- raster("data/LULC/LU_rast_30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
lulc <- crop(lulc, aoi)
lulc <- mask(lulc, aoi)
biom.dat[,lulc:=getValues(lulc)]
hist(biom.dat[,lulc])

## noise testing to look at spread around the model coefficients here
### select randomly the coefficients to within +/- 1 sterr of the model estimates
mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]

b0.bai.range <- seq(mod.int.edge.fin$coefficients[1,1]-mod.int.edge.fin$coefficients[1,2], 
                    mod.int.edge.fin$coefficients[1,1]+mod.int.edge.fin$coefficients[1,2],
                    length.out=50)
b1.bai.range <- seq(mod.int.edge.fin$coefficients[2,1]-mod.int.edge.fin$coefficients[2,2], 
                    mod.int.edge.fin$coefficients[2,1]+mod.int.edge.fin$coefficients[2,2],
                    length.out=50)
b2.bai.range <- seq(mod.int.edge.fin$coefficients[3,1]-mod.int.edge.fin$coefficients[3,2], 
                    mod.int.edge.fin$coefficients[3,1]+mod.int.edge.fin$coefficients[3,2],
                    length.out=50)
dbh.dump <- andy.dbh
biom.dump <- biom.dat
npp.track <- numeric()
beta.track <- list()
edge.track <- numeric()
int.track <- numeric()
for(i in 1:2000){
  b0.rand <- sample(b0.bai.range, 1)
  b1.rand <- sample(b1.bai.range, 1)
  b2.rand <- sample(b2.bai.range, 1)
  dbh.dump[seg==10, growth:=biom*exp(b0.rand+(b1.rand*log(dbh)))]
  dbh.dump[seg%in%c(20,30), growth:=biom*exp(b0.rand+(b1.rand*log(dbh))+b2.rand)]
  g <- dbh.dump[, .(sum(growth), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
  g[,rel.gain:=V1/V2] ## this deals in forest growth as a function of biomass, not of forest area
  edge.factor <- g[seg==10, mean(rel.gain)]
  int.factor <- g[seg!=10, mean(rel.gain)]
  biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
  biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
  biom.dump[,npp.tot:=npp.edge+npp.int]
  biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
  npp.track <- c(npp.track, biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)])
  beta.track[[i]] <- c(b0.rand, b1.rand, b2.rand)
  edge.track <- c(edge.track, edge.factor)
  int.track <- c(int.track, int.factor)
  print(i)
}

npp.track[which(npp.track==min(npp.track))]/(1E3*2) ### minimum of 8.3k tC
npp.track[which(npp.track==max(npp.track))]/(1E3*2) ## max of 22.7k tC
b0.bai #-0.037
b1.bai #-0.916
b2.bai #-0.490
beta.track[[which(npp.track==min(npp.track))]] # -0.28, -0.99 -0.51
beta.track[[which(npp.track==max(npp.track))]] # 0.21, -0.84, -0.45
a <- edge.track-int.track
max(a) ## can be up to 2.84% apart
min(a) # can be as little as 0.7% apart

## at maximum plausible edge vs. interior productivity
edge.factor <- max(edge.track)
int.factor <- min(int.track)
biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
biom.dump[,npp.tot:=npp.edge+npp.int]
biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)]*1E-3*(1/2) ## 19.8k tC at maximum separation

## at minimum plausible edge vs. interior productivity
edge.factor <- min(edge.track)
int.factor <- max(int.track)
biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
biom.dump[,npp.tot:=npp.edge+npp.int]
biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)]*1E-3*(1/2) ## 11.7k tC at minimum separation





######
###### Approach 5: Reconfigure street tree analysis for more effective search
## Look at cell biomass and give a reasonable subset of the dbh data that will not allow excessive density or BA to reach the biomass total
## street tree simulator + results reconstruction
################
#########
### final script for estimating tree distribution and running on the cluster
library(data.table)
library(raster)

### biomass equation biom~dbh
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## get the street tree record prepped
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
street[, record.good:=0] ## ID records with good data quality
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, npp.ann:=(biom.2014-biom.2006)/8]
street[record.good==1, npp.ann.rel:=npp.ann/biom.2006]
street[record.good==1,range(dbh.2006, na.rm=T)] 
street[record.good==1, range(npp.ann, na.rm=T)] ## 202 records show negative growth -- questionable 2006 dbh record, cull
street[record.good==1 & npp.ann<0, record.good:=0] ## 2401
# street[dbh.2006<5, record.good:=0] ## filter for tiny trees

### biomass increment model based on dbh (exponential function)
# mod.biom.rel <- summary(lm(log(npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & npp.ann!=0])) ## exclude 4 no-growth records
# plot(street[record.good==1, log(biom.2006)], street[record.good==1, log(npp.ann.rel)]) # r2 = 0.54, slope significant
## there's <10 records below 5cm with a lot of leverage, cut them
mod.biom.rel <- summary(lm(log(npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & npp.ann!=0 & dbh.2006>=5,])) ## exclude 4 no-growth records
# plot(street[record.good==1 & dbh.2006>=5, log(biom.2006)], street[record.good==1 & dbh.2006>=5, log(npp.ann.rel)]) # r2 = 0.55, looks nicer
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees


## read in biomass and canopy data
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,index:=1:dim(biom.dat)[1]] ## master pixel index for stack of biom/aoi/can, 354068 pix
p <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.biom30m]) ## cdf of cell biomass

## prep street tree data and biomass data for processing (small biomass first)
ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
clean <- street[record.good==1,] # get a good street tree set ready
clean[, ba:=ba.pred(dbh.2006)]
setkey(clean, biom.2006) # 2390 records in final selection, dbh range 5-112 cm
clean <- clean[order(clean$dbh.2006),]
clean[,rank:=seq(from=1, to=dim(clean)[1])] ## all the trees have a fixed size rank now (used in adjusting sampling weights)


## set up the cells to simulate
# hist(biom.dat[bos.biom30m>0, bos.biom30m])
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,] ## ignore extremely low biomass cells, 108k

## some thoughts about how to focus the sampling along the distribution
# ba.pred.cell <- function(x, A, C){(((x/2)^2)*pi)/(A*C)} ## x is dbh (cm), A is AOI (cell size m2), C is canopy fraction (0-1)
# quantile(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.biom30m], probs=c(.10,.50,.90)) ## for the cells we are running, 235, 3740, 17537 kg
# 20000/clean[, median(biom.2006, na.rm=T)] ## this could be 77 median trees in a 20k cell
# ba.pred.cell(clean[,median(dbh.2006, na.rm=T)], 900, 0.6) ## 25cm median dbh, 60% Can, BA is 0.9
# (20000/clean[, median(biom.2006, na.rm=T)])*ba.pred.cell(clean[,median(dbh.2006, na.rm=T)], 900, 0.6) ## 72.7 m2/ha, proposed BA
## so in upper reaches a BA restriction of ~40 m2/ha will help to keep bitty trees down
## maybe start the selection function around what sample will give you a median dbh that will let you get in under the BA cap
# ### what about numerical density (stems/ha)
# ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
# andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv")
# andy.dbh <- as.data.table(andy.dbh)
# names(andy.dbh) <- c("Tree.ID", "Plot.ID", "Spp", "X", "Y", "dbh", "Can.class")
# andy.dbh[Y<10, seg:=10]
# andy.dbh[Y>=10 & Y<20, seg:=20]
# andy.dbh[Y>=20, seg:=30]
# andy.dbh[,ba:=ba.pred(dbh)]
# d <- andy.dbh[,length(dbh), by=.(Plot.ID, seg)]
# d[seg==10, mean(V1)]/(20*10) ## 0.106 trees/m2, about 9.4m2/tree
# d[seg==20, mean(V1)]/(20*10) ## 0.077 trees/m2
# d[seg==30, mean(V1)]/(20*10) ## 0.083 trees/m2
## let's also impose a no-more-than 0.106 trees/m2 rule
# (20000/clean[, median(biom.2006, na.rm=T)])/(900*0.6) ## ok this is 0.143 trees/m2, too high
## 2 prongs: adjust the sampling interval based on BA cap, and the sample size based on density cap

### test samples
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>10 & !is.na(aoi) & !is.na(bos.can30m) & bos.can30m>0,] ## 108k
setkey(runme, index)
# runme <- runme[bos.biom30m>30000,]
runme <- runme[sample(1:dim(runme)[1], size=1000),]

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
cage.genus <- list()
index.track <- integer()
biom.track <- numeric()
cage.wts <- list()
cage.biom.sim <- list()
attempts.track <- integer()
proc.track <- integer()

## for dynamic weighting
incr=300 ## incrememnt to change weighting
k=200 ## constant, arbitrary top of prob weights, higher-->steeper change in sampling weight
D=k/incr ## drop increment, sensitive to number of tries to make while adjusting the sampling weights

## loop each row
for(t in 1:dim(runme)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  biom.sim.track <- numeric()
  wts.track <- integer()
  cage.genus[[t]] <- list()
  cage.dbh[[t]] <- list()
  
  x <- 0 ## count successful simulations
  q <- 0 ## count attempts since last successful sample
  g <- 0 ## weight adjustment start
  z <- 0 ## total number of sample iterations
  jeez <- 0 ## if it's grinding through a slow one
  
  ### this approach uses a moving/shrinking window on the street tree records as biomass starts to get bigger
  ## has the disadvantage of supressing the number of trees at higher biomass (start to oversample very large trees), would need some tuning to get the right sampling window
  # ## 0 figure where and how much to sample
  # cell.q <- p(runme[t, bos.biom30m]) ## this is the percentile of the cell biomass
  # spread <- min((1-cell.q), 0.5) ## interval to sample in dbh records, cap it off if the spread gets bigger than 50%
  # window.l <- max(0, (cell.q-spread))
  # window.h <- cell.q+spread
  # if((spread)==0.5){window.h <- 1} ## basic scheme is that for biomass>median, restrict to (rarer) higher window, but if biom<median, sample whole window
  # dbh.lims <- quantile(clean$dbh.2006, probs=c(window.l, window.h)) ## window of dbh to select based on cell biomass
  # grasp <- clean[biom.2006<runme[t, bos.biom30m] &dbh.2006<=dbh.lims[2] & dbh.2006>=dbh.lims[1], .(dbh.2006, biom.2006, ba)] ## this is a street tree record restricted based on cell biomass
  
  dens.lim <- ceiling(0.106*runme[t,aoi*bos.can30m]) ## maximum number to select, probably only restrictive in low-canopy cells
  
  ### this approach takes the full street trees record and samples it with differing weights if the algorithm fails
  grasp <- clean[biom.2006<(1.1*runme[t, bos.biom30m]), .(dbh.2006, biom.2006, ba, rank, genus)] ## excluding single trees too big to fit cell biomass
  
  ### attempt to sample street trees to fill the cells
  while(x<100 & q<1000){ ## select 100 workable samples, or quit after q attempts with no success
    # ## grab a sample of grasp (no weighting)
    # samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=T),] ## sample grasp with replacement only up to density limit
    
    ## diagnostic plots -- does the weighting scheme really push the distribution north? yes, slowly
    # med <- numeric()
    # for(god in 1:300){
    #   wts=(k-(god*D))+(((god*D)/(dim(grasp)[1]-1))*((grasp$rank)-1))
    #   cock <- numeric()
    #   for(j in 1:500){
    #     samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=T, prob = wts),] 
    #     cock <- c(cock, median(samp$dbh.2006))
    #   }
    #   med <- c(med, median(cock))
    #   # hist(cock, main=paste("weights=", god), xlim=c(20, 45), breaks=13)
    # }
    # plot(1:300, med)
    # plot(grasp$rank, wts)
    # plot(grasp$dbh.2006, wts)

    ## OR: figure out the dynamic weighting to use based on previous failures
    wts=(k-(g*D))+(((g*D)/(dim(grasp)[1]-1))*((grasp$rank)-1))
    # sample grasp with replacement only up to density limit, with specified weights based on iterative reweighting
    samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=F, prob = wts),] 
    w=samp[1, biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme[t, bos.biom30m]) & d<dens.lim){ ## keep adding trees until you just get over the target biomass or run out of records
      w=w+samp[d+1, biom.2006]
      d=d+1
    }
    ### if this is too many trees too tightly packed (over 40m2/ha BA) or not enough biomass in the sample, readjust weights
    if(((samp[1:d, sum(ba)]/runme[t, aoi*bos.can30m])*1E4)>40 | w<(0.90*runme[t, bos.biom30m])){
      if(g<incr){
        g=g+1 ## readjust the sample weights
        if(g==300){print(paste("pixel", runme[t, index], "weights at maximum"))}
        q=0 ## reset attempt timeout clock
      }
    }
    if(((samp[1:d, sum(ba)]/runme[t, aoi*bos.can30m])*1E4)<40 & w<(1.10*runme[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
      x <- x+1 ## record successful sample
      ann.npp <- c(ann.npp, sum(samp[1:d, biom.2006]*exp((mod.biom.rel$coefficients[2]*log(samp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      biom.sim.track <- c(biom.sim.track, w) ## simulated biomass in this sample
      cage.dbh[[t]][[x]] <- samp[1:d, dbh.2006] ## which dbhs did you select
      cage.genus[[t]][[x]] <- samp[1:d, genus] # running list of which genera you select
      wts.track <- c(wts.track, g) ## weights used to get this sample
      if(q>300){print(paste("finally got sample", x, "after", q, "tries")); jeez <- 1} ## to give status reports for difficult/long-process samples
      q=0 ## reset counter if you've found some solution
      if(jeez==1 & x%%10==0){print(paste("patience, just got sample", x))}
      if(x>30 & g>250){q <- (-500)} ## give more room if this is a hard one and you've managed by chance to get a pretty good sample going
    }
    q=q+1 ## record number of attempts since last succesful simulation
    if(q>600 & q%%100==0){print(paste("attempt clock out in", (1000-q)/100))}
    z=z+1 ## record total number of sample iterations
  }
  cage.ann.npp[[t]] <- ann.npp
  cage.num.trees[[t]] <- num.trees
  cage.wts[[t]] <- wts.track
  biom.track <- c(biom.track, runme[t, bos.biom30m])
  index.track <- c(index.track, runme[t, index])
  cage.biom.sim[[t]] <- biom.sim.track
  attempts.track <- c(attempts.track, z)
  if(q<1000){
    proc.track <- c(proc.track, 1)
    print(paste("finished pixel", runme[t, index]))
  } else{
    proc.track <- c(proc.track, 0) ## record as incomplete processing
    print(paste("pixel", runme[t, index], "failed"))
  }
}

# ## basic diagnostics
an <- sapply(cage.ann.npp, FUN=median)
hist(an)
n <- sapply(cage.num.trees, FUN=median)
hist(n)
w <- sapply(cage.wts, FUN=max) ## most seem to punch through without having to reweight the street records
plot(biom.track, sapply(cage.wts, FUN=max)) ## it's real steady from 1-10k, then you start to get a lot more reweights at high and at very low biomass
hist(biom.track)
hist(attempts.track)
plot(attempts.track, sapply(cage.wts, FUN=max))

## for example:
hist(sapply(cage.dbh[[66]], FUN=median))
hist(sapply(cage.ann.npp[[66]], FUN=median))
hist(sapply(cage.num.trees[[66]], FUN=median))
biom.track[66]
cage.genus[[66]]
cage.dbh[[66]]

### which ones fail
biom.track[proc.track==0]
cage.ann.npp[c(which(proc.track==0))] ## are you finding partial retreivals?
cage.dbh[which(proc.track==0)]
cage.ann.npp[which(proc.track==0)]

## save the objects to disk for later reconstituting
save(cage.ann.npp, file=paste("processed/boston/ann.npp.street.v3.weighted"))
save(cage.num.trees, file=paste("processed/boston/num.trees.street.v3.weighted"))
save(cage.dbh, file=paste("processed/boston/dbh.street.v3.weighted"))
save(biom.track, file=paste("processed/boston/biom.track.street.v3.weighted"))
save(cage.ann.npp, file=paste("processed/boston/index.track.street.v3.weighted"))




