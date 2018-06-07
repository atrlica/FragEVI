library(raster)
library(data.table)
library(ggplot2)


# ### Canopy area by cumulative distance from edge
# bos.canF <- raster("processed/boston/bos_can01_filt.tif")
# # bos.can <- raster("data/dataverse_files/bostoncanopy_1m.tif")
# 
# can.sum <- function(x) { # x is canopy 0/1 1m raster object
#   bs <- blockSize(x)
#   y <- integer()
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
#     v[v>1] <- NA  # for some reason some of the NA's are getting labeled as 120
#     y <- c(y, sum(v, na.rm=T))
#     print(paste("finished block", i, "of", bs$n))
#   }
#   return(y)
# }
# 
# ### get the total area of canopy (excluding tiny canopy gaps that are filtered in buffer calculations)
# can.buffs <- list.files("processed/boston/")
# can.buffs <- can.buffs[grep(pattern = "bos.nocan_", x = can.buffs)]
# can.buffs <- can.buffs[grep(pattern=".tif", x=can.buffs)]
# can.buffs <- can.buffs[!grepl(pattern = ".vat", x=can.buffs)]
# can.buffs <- can.buffs[!grepl(pattern = ".aux", x=can.buffs)]
# can.buffs <- can.buffs[!grepl(pattern = ".ovr", x=can.buffs)]
# buff.dist <- as.integer(unlist(rm_between(can.buffs, "nocan_", "mbuff", extract=TRUE)))
# 
# ### prep masking rasters
# ## all of Boston
# # towns <- readOGR(dsn = "F:/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# # bos.AOI <- towns[towns@data$TOWN=="BOSTON",]
# # bos.AOI <- bos.AOI[bos.AOI@data$SHAPE_AREA>1E07,] ## remove Harbor Islands
# # bos.AOI <- spTransform(bos.AOI, crs(bos.canF))
# # bos.AOI@data$include <- 1
# # bos.AOI.r <- rasterize(bos.AOI[bos.AOI@data$include], bos.canF)
# # ba.dat <- as.data.table(as.data.frame(bos.AOI.r))
# # ba.dat[!is.na(OBJECTID) | !is.na(OBJECTID.1), dogshit:=1] ## set all areas to same value
# # bos.AOI.r <- setValues(bos.AOI.r, ba.dat$dogshit)
# # writeRaster(bos.AOI.r, filename="processed/boston/bos.AOI.1m.tif", format="GTiff", overwrite=T)
# 
# ### whole area of canopy first
# print("getting total canopy and area for all Boston")
# ma <- raster("processed/boston/bos.AOI.1m.tif")
# # bos.canF<- mask(bos.canF, ma) ## masking takes awhile and isn't necessary for whole-city -- use mask for full area calc
# tot <- integer()
# dog <- can.sum(bos.canF)
# tot <- c(tot, sum(dog, na.rm=T))
# dist <- 0 ## keep track of the buffer distances as you process the rasters in whatever random order they come in
# area.tot <- can.sum(ma)
# area.tot <- sum(area.tot, na.rm=T)
# 
# ## now load up each 0/1 classed canopy raster (1 represents canopy interior to the indicated buffer distance -- areas diminish with greater buffer)
# for(g in 1:length(can.buffs)){
#   print(paste("working on", can.buffs[g]))
#   r <- raster(paste("processed/boston/", can.buffs[g], sep=""))
#   dog <- can.sum(r)
#   tot <- c(tot, sum(dog, na.rm=T))
#   dist <- c(dist, buff.dist[g])
# }
# results <- cbind(dist, tot)
# results <- results[order(dist),]
# results <- as.data.frame(results)
# colnames(results) <- c("dist", "pix.more.than")
# results$pix.less.than <- results$pix.more.than[1]-results$pix.more.than ## recall that this method completely leaves out any gap areas that are <50m2 -- neither counted as canopy nor as gap area
# results$less.rel <- results$pix.less.than/results$pix.more.than[1]
# results$frac.tot.area <- results$pix.less.than/area.tot
# plot(results$dist, results$less.rel) ## cumulative distribution of canopy edge
# write.csv(results, "processed/bos.can.cummdist.csv")


### pie chart, relative canopy distance fraction
dist <- read.csv("processed/bos.can.cummdist.csv")
library(data.table)
dist <- as.data.table(dist)
ed.all <- dist[dist==0, pix.more.than]
ed.10m <- dist[dist==10, pix.less.than]
ed.20m <- dist[dist==20, pix.less.than]
ed.30m <- dist[dist==30, pix.less.than]
ed.int <- dist[dist==30, pix.more.than]
ed.all
ed.10m/ed.all
ed.20m.only <- ed.20m-ed.10m
ed.20m.only/ed.all
ed.30m.only <- ed.30m-ed.20m
ed.30m.only/ed.all
ed.int/ed.all

slices <- c(ed.10m/ed.all, ed.20m.only/ed.all, ed.30m.only/ed.all, ed.int/ed.all)*100

par(mar=c(1, 1,3, 4))

pie(slices, labels=paste(c("<10m, ", "10-20m, ", "20-30m, ", ">30m, "), round(slices, 0), "%", sep=""),
    main="", font=2, cex=2,
    col = c("salmon", "orange", "lightgoldenrod1", "green3"))
mtext("Fraction of canopy per distance class", side = 3, cex=2.3, font=2)

### do same canopy area calcs for buffers in specific sub-areas
hoods <- c("downtown", "jamaica", "allston", "southie", "dorchester", "hydepark", "common")


### masks for individual test AOIs
for(d in 1:length(hoods)){
  print(paste("rasterizing", hoods[d]))
  clipme <- readOGR(dsn=paste("processed/zones/bos.", hoods[d], ".shp", sep=""), layer=paste("bos.", hoods[d], sep=""))
  # rasterize(clipme, bos.canF, format="GTiff", overwrite=T, filename=paste("processed/zones/bos.", hoods[d], ".tif", sep=""))
  ma <- raster(paste("processed/zones/bos.", hoods[d], ".tif", sep=""))
  ma <- crop(ma, extent(clipme))
  writeRaster(ma, format="GTiff", overwrite=T, filename=paste("processed/zones/bos.", hoods[d], ".tif", sep=""))
}

### get 0 buffer canopy area first, then mask as you go
for(d in 1: length(hoods)){
  print(paste("initializing", hoods[d]))
  rm(results)
  tot <- integer()
  dist <- 0
  ma <- raster(paste("processed/zones/bos.", hoods[d], ".tif", sep=""))
  area.tot <- sum(getValues(ma), na.rm=T) ## total size of the AOI
  r <- raster("processed/boston/bos_can01_filt.tif")
  r <- crop(r, ma)
  r <- mask(r, ma)
  dog <- can.sum(r)
  tot <- c(tot, sum(dog, na.rm=T))
  
  for(g in 1:length(can.buffs)){
    print(paste("working on", can.buffs[g], "in", hoods[d]))
    r <- raster(paste("processed/boston/", can.buffs[g], sep=""))
    r <- crop(r, ma)
    r <- mask(r, ma)
    dog <- can.sum(r)
    tot <- c(tot, sum(dog, na.rm=T))
    dist <- c(dist, buff.dist[g])
  }
  results <- cbind(dist, tot)
  results <- results[order(dist),]
  results <- as.data.frame(results)
  colnames(results) <- c("dist", "pix.more.than")
  results$pix.less.than <- results$pix.more.than[1]-results$pix.more.than ## recall that this method completely leaves out any gap areas that are <50m2 -- neither counted as canopy nor as gap area
  results$less.rel <- results$pix.less.than/results$pix.more.than[1]
  results$frac.tot.area <- results$pix.less.than/area.tot
  write.csv(results, paste("processed/bos.can.cummdist.", hoods[d], ".csv", sep=""))
}

### combined plot, canopy edge area as fraction of total area
results <- read.csv("processed/bos.can.cummdist.csv")
hoods <- c("downtown", "jamaica", "allston", "southie", "dorchester", "hydepark", "common")
cols=rainbow(length(hoods))
plot(results$dist, results$frac.tot.area, pch=1, col="black", type="l", lwd=3, 
     xlab="distance from edge (m)", ylab="area fraction",
     ylim=c(0, 0.65))
for(d in 1:length(hoods)){
  dat <- read.csv(paste("processed/bos.can.cummdist.", hoods[d], ".csv", sep=""))
  lines(dat$dist, dat$frac.tot.area, col=cols[d], type="l", lwd=2)
}
legend(x=60, y=0.4, bty = "n", legend=c("Boston", hoods), fill=c("black", cols))



######
### canopy edge area cumulative, extract by LULC collapsed classes
bos.forest <- raster("processed/boston/bos.forest_only.tif")
bos.dev <- raster("processed/boston/bos.dev_only.tif")
bos.hdres <- raster("processed/boston/bos.hdres_only.tif")
bos.ldres <- raster("processed/boston/bos.ldres_only.tif")
bos.lowveg <- raster("processed/boston/bos.lowveg_only.tif")
bos.water <- raster("processed/boston/bos.water_only.tif")
# 
# for.sum <- sum(getValues(bos.forest), na.rm=T)
# dev.sum <- sum(getValues(bos.dev), na.rm=T)
# hdres.sum <- sum(getValues(bos.hdres), na.rm=T)
# ldres.sum <- sum(getValues(bos.ldres), na.rm=T)
# lowveg.sum <- sum(getValues(bos.lowveg), na.rm=T)
# water.sum <- sum(getValues(bos.water), na.rm=T)
# bos.aoi <- raster("processed/boston/bos.AOI.1m.tif")
# aoi.sum <- sum(getValues(bos.aoi), na.rm=T)
# for.sum/aoi.sum ## 8.2%
# dev.sum/aoi.sum ## 38%
# hdres.sum/aoi.sum ## 38.7%
# ldres.sum/aoi.sum ## 2.0%
# lowveg.sum/aoi.sum # 10.7%
# water.sum/aoi.sum # 2.1%


## sum of canopy area, masked by lulc
can.sum.ma <- function(x,m) { # x is canopy 0/1 1m raster object, m is mask
  bs <- blockSize(x)
  y <- integer()
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
    z <- getValues(m, row=bs$row[i], nrows=bs$nrows[i])
    v[v>1 | v<0] <- NA  # for some reason some of the NA's are getting labeled as 120
    v[z!=1] <- 0 ## cancel values outside mask area
    y <- c(y, sum(v, na.rm=T))
    print(paste("finished block", i, "of", bs$n))
  }
  return(y)
}

### get the total area of canopy (excluding tiny canopy gaps that are filtered in buffer calculations)
can.buffs <- list.files("processed/boston/")
can.buffs <- can.buffs[grep(pattern = "bos.nocan_", x = can.buffs)]
can.buffs <- can.buffs[grep(pattern=".tif", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".vat", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".aux", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".ovr", x=can.buffs)]
buff.dist <- as.integer(unlist(rm_between(can.buffs, "nocan_", "mbuff", extract=TRUE)))


can.master <- raster("processed/boston/bos_can01_filt.tif")
can.tot <- sum(getValues(can.master), na.rm=T) ### 31.7% of boston raster is canopy
aoi.master <- raster("processed/boston/bos.aoi.tif")
aoi.tot <- sum(getValues(aoi.master), na.rm=T)
## cumulative canopy area by edge distance, for each lulc class
lu.classes <- c("dev", "hdres", "ldres", "lowveg", "forest")
for(l in 1:length(lu.classes)){
  print(paste("initializing", lu.classes[l]))
  if(exists("results")){rm(results)} ## clean up previous store file
  tot <- integer()
  dist <- 0
  ma <- raster(paste("processed/boston/bos.", lu.classes[l], "_only.tif", sep=""))
  area.tot <- sum(getValues(ma), na.rm=T) ## total size of the lulc target
  dog <- can.sum.ma(can.master, ma)
  tot <- c(tot, sum(dog, na.rm=T))
  
  for(g in 1:length(can.buffs)){
    print(paste("working on", can.buffs[g], "in", lu.classes[l]))
    r <- raster(paste("processed/boston/", can.buffs[g], sep=""))
    dog <- can.sum.ma(r, ma)
    tot <- c(tot, sum(dog, na.rm=T))
    dist <- c(dist, buff.dist[g])
  }
  results <- cbind(dist, tot)
  results <- results[order(dist),]
  results <- as.data.frame(results)
  colnames(results) <- c("dist", "pix.more.than")
  results$pix.less.than <- results$pix.more.than[1]-results$pix.more.than ## recall that this method completely leaves out any gap areas that are <50m2 -- neither counted as canopy nor as gap area
  results$less.rel <- results$pix.less.than/results$pix.more.than[1]
  results$frac.tot.area <- results$pix.less.than/area.tot
  write.csv(results, paste("processed/bos.can.cummdist.", lu.classes[l], ".csv", sep=""))
}

### combined plot, canopy edge area as fraction of total area (by LULC class)
results <- read.csv("processed/bos.can.cummdist.csv")
lu.classes <- c("forest", "dev", "hdres", "ldres", "lowveg")
cols=rainbow(length(lu.classes))
cols=c("forestgreen", "blue", "red")
par(mar=c(4.5, 5.5, 1, 1), oma=c(0,0,0, 0), xpd=F)

plot(results$dist, results$frac.tot.area, pch=1, col="black", type="l", lwd=7, bty="n", lty=1, 
     xlab="Distance from edge (m)", ylab="Cummulative area fraction",
     ylim=c(0, 0.9), xlim=c(0, 60), yaxt="n", font.lab=2, cex.lab=2, cex.axis=2)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=c("0", "20%", "40%", "60%", "80%"), cex.axis=2)
for(d in 1:length(lu.classes)){
  dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[d], ".csv", sep=""))
  lines(dat$dist, dat$frac.tot.area, col=cols[d], type="l", lwd=3)
}
legend("right", x=40, y=0.80, cex=2, legend=c("Boston", "Forest", "Developed", "Residential"), fill=c("black", cols), bty="n")

### combined plot, canopy edge area as fraction of total area (by LULC class)
results <- read.csv("processed/bos.can.cummdist.csv")
lu.classes <- c("forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other")
cols=rainbow(length(lu.classes))
cols=c("forestgreen", "blue", "red")
par(mar=c(4.5, 5.5, 1, 1), oma=c(0,0,0, 0), xpd=F)

plot(results$dist, results$frac.tot.area, pch=1, col="black", type="l", lwd=7, bty="n", lty=1, 
     xlab="Distance from edge (m)", ylab="Cummulative area fraction",
     ylim=c(0, 0.9), xlim=c(0, 60), yaxt="n", font.lab=2, cex.lab=2, cex.axis=2)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=c("0", "20%", "40%", "60%", "80%"), cex.axis=2)
for(d in 1:length(lu.classes)){
  dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[d], ".csv", sep=""))
  lines(dat$dist, dat$frac.tot.area, col=cols[d], type="l", lwd=3)
}
legend("right", x=40, y=0.80, cex=2, legend=c("Boston", "Forest", "Developed", "Residential"), fill=c("black", cols), bty="n")




lu.classes <- c("forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other")
for.dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[1], ".csv", sep=""))
hd.dat <-  read.csv(paste("processed/bos.can.cummdist.", lu.classes[3], ".csv", sep=""))
md.dat <-  read.csv(paste("processed/bos.can.cummdist.", lu.classes[2], ".csv", sep=""))
dev.dat <-  read.csv(paste("processed/bos.can.cummdist.", lu.classes[4], ".csv", sep=""))
aoi <- raster("processed/boston/bos.AOI.1m.tif")
aoi.dat <- as.data.table(as.data.frame(aoi))


## raw pixel count by LULC at each distance class
bos.forest <- raster("processed/boston/bos.forest.tif")
bos.dev <- raster("processed/boston/bos.dev.tif")
bos.hd.res <- raster("processed/boston/bos.hd.res.tif")
bos.med.res <- raster("processed/boston/bos.med.res.tif")
bos.low.res <- raster("processed/boston/bos.low.res.tif")
bos.lowveg <- raster("processed/boston/bos.lowveg.tif")
bos.other <- raster("processed/boston/bos.other.tif")
bos.water <- raster("processed/boston/bos.water.tif")
bos.canF <- raster("processed/boston/bos_can01_filt.tif")

## sum of canopy area, masked by lulc
can.sum.ma <- function(x,m) { # x is canopy 0/1 1m raster object, m is mask
  bs <- blockSize(x)
  y <- integer()
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
    z <- getValues(m, row=bs$row[i], nrows=bs$nrows[i])
    v[v>1 | v<0] <- NA  # for some reason some of the NA's are getting labeled as 120
    v[z!=1] <- 0 ## cancel values outside mask area
    y <- c(y, sum(v, na.rm=T))
    print(paste("finished block", i, "of", bs$n))
  }
  return(y)
}


contain <- data.frame(distance=integer(), LULC=character(), pix.num=integer())
for(g in 1:length(lu.classes)){
  dist.dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[g], ".csv", sep=""))
  # m <- raster(paste("processed/boston/bos.", lu.classes[g], ".tif", sep=""))
  # r <- raster("processed/boston/bos_can01_filt.tif") ## full canopy layer
  # tot.can <- can.sum.ma(r, m)
  # a <- sum(tot.can, na.rm=T) ## total number of canopy pixels in this LULC class
  pix.d <- 0
  for(e in 2:dim(dist.dat)[1]){
    pix.d <- c(pix.d, dist.dat$pix.less.than[e]-dist.dat$pix.less.than[(e-1)])
  }
  contain <- rbind(contain, cbind(seq(0,100), rep(lu.classes[g], 101), pix.d))
}
colnames(contain) <- c("distance", "LULC", "pix.num")
contain$distance <- as.integer(as.character(contain$distance))
contain$pix.num <- as.integer(as.character(contain$pix.num))

library(ggplot2)
contain <- as.data.table(contain)
ggplot(contain[distance!=0,], aes(x=distance, y=pix.num, fill=LULC)) + 
  geom_area(aes(color=LULC, fill=LULC))+
  xlim(1,50)