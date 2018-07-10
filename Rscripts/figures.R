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

### processing script to get cumulative area by edge distance within LULC
# ## sum of canopy area, masked by lulc
# can.sum.ma <- function(x,m) { # x is canopy 0/1 1m raster object, m is mask
#   bs <- blockSize(x)
#   y <- integer()
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
#     z <- getValues(m, row=bs$row[i], nrows=bs$nrows[i])
#     v[v>1 | v<0] <- NA  # for some reason some of the NA's are getting labeled as 120
#     v[z!=1] <- 0 ## cancel values outside mask area
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
# 
# can.master <- raster("processed/boston/bos_can01_filt.tif")
# can.tot <- sum(getValues(can.master), na.rm=T) ### 31.7% of boston raster is canopy
# aoi.master <- raster("processed/boston/bos.aoi.tif")
# aoi.tot <- sum(getValues(aoi.master), na.rm=T)
# ## cumulative canopy area by edge distance, for each lulc class
# lu.classes <- c("forest", "dev", "hdres", "ldres", "lowveg")
# for(l in 1:length(lu.classes)){
#   print(paste("initializing", lu.classes[l]))
#   if(exists("results")){rm(results)} ## clean up previous store file
#   tot <- integer()
#   dist <- 0
#   ma <- raster(paste("processed/boston/bos.", lu.classes[l], "_only.tif", sep=""))
#   area.tot <- sum(getValues(ma), na.rm=T) ## total size of the lulc target
#   dog <- can.sum.ma(can.master, ma)
#   tot <- c(tot, sum(dog, na.rm=T))
#   
#   for(g in 1:length(can.buffs)){
#     print(paste("working on", can.buffs[g], "in", lu.classes[l]))
#     r <- raster(paste("processed/boston/", can.buffs[g], sep=""))
#     dog <- can.sum.ma(r, ma)
#     tot <- c(tot, sum(dog, na.rm=T))
#     dist <- c(dist, buff.dist[g])
#   }
#   results <- cbind(dist, tot)
#   results <- results[order(dist),]
#   results <- as.data.frame(results)
#   colnames(results) <- c("dist", "pix.more.than")
#   results$pix.less.than <- results$pix.more.than[1]-results$pix.more.than ## recall that this method completely leaves out any gap areas that are <50m2 -- neither counted as canopy nor as gap area
#   results$less.rel <- results$pix.less.than/results$pix.more.than[1]
#   results$frac.tot.area <- results$pix.less.than/area.tot
#   write.csv(results, paste("processed/bos.can.cummdist.", lu.classes[l], ".csv", sep=""))
# }


################
### combined plot, canopy edge area as fraction of total area (by LULC class)
results <- read.csv("processed/bos.can.cummdist.csv")
lu.classes <- c("forest", "dev", "hdres", "ldres", "lowveg")
cols=rainbow(length(lu.classes))
cols=c("forestgreen", "gray55", "red", "yellow", "green")
par(mar=c(4.5, 5.5, 1, 1), oma=c(0,0,0, 0), xpd=F)

plot(results$dist, results$frac.tot.area, pch=1, col="black", type="l", lwd=7, bty="n", lty=1, 
     xlab="Distance from edge (m)", ylab="Cummulative area fraction",
     ylim=c(0, 0.9), xlim=c(0, 60), yaxt="n", font.lab=2, cex.lab=2, cex.axis=2)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=c("0", "20%", "40%", "60%", "80%"), cex.axis=2)
for(d in 1:length(lu.classes)){
  dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[d], ".csv", sep=""))
  lines(dat$dist, dat$frac.tot.area, col=cols[d], type="l", lwd=3)
}
legend("right", x=40, y=0.80, cex=1, legend=c("Boston", "Forest", "Developed", "HDRes", "LDRes", "LowVeg"), fill=c("black", cols), bty="n")






###
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

contain$ha <- contain$pix.num/(1E4)


library(ggplot2)
contain <- as.data.table(contain)
contain <- contain[distance!=0,]
ggplot(contain[distance!=0,], aes(x=distance, y=ha, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  xlim(1,50)+
  scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3"),
                    name="Land Cover",
                    breaks=c("forest", "dev", "hdres", "ldres", "lowveg"),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg."))+
  xlab("Canopy distance from edge (m)")+
  ylab("Total canopy area (ha)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = c(0.8, 0.2))


#### FIA density plots
library(raster)
library(data.table)
a=123.67
b=0.04
AGE=seq(0,150)
y = a*(1-exp(-b*AGE))^3
plot(AGE, y)
z = diff(y)
z.rel <- z/y[2:151]
length(y)
length(z) ## this is growth after 1 year, 2 years, etc
plot(y[2:151], z)
points(y[2:151], z.rel, pch=14, cex=0.5, col="red") #hyperbolic
plot(AGE[2:151], z) ## this is the gain curve over site age

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
can <- raster("processed/boston/bos.can30m.tif")
isa <- raster("processed/boston/bos.isa.rereg30m.tif")
biom <- crop(biom, isa)
aoi <- crop(aoi, isa)
can <- crop(can, isa)
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can.dat <- as.data.table(as.data.frame(can))
biom.dat[, can.frac:=can.dat$bos.can30m]
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat[, isa.frac:=isa.dat$bos.isa.rereg30m]

### live MgC by area of GROUND in each pixel
biom.dat[, live.MgC.ha.ground:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[aoi>800, range(live.MgC.ha.ground, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
hist(biom.dat[aoi>800,live.MgC.ha.ground]) ## correcting for canopy cover, more mid-rage values
### live MgC/ha for "forest" fraction in each pixel
biom.dat[,live.MgC.ha.forest:=(bos.biom30m/(aoi*can.frac))*(1/2)*(1/1000)*(10^4)]
biom.dat[can.frac<=0.01, live.MgC.ha.forest:=0]
range(biom.dat[aoi>800,live.MgC.ha.forest],  na.rm=T) ## 0 - 284
hist(biom.dat[aoi>800,live.MgC.ha.forest]) ## correcting for canopy cover, more mid-rage values
## live MgC/ha in the pixel pervious cover
biom.dat[,live.MgC.ha.perv:=(bos.biom30m/(aoi*(1-isa.frac)))*(1/2)*(1/1000)*(10^4)]
biom.dat[isa.frac>0.99, live.MgC.ha.perv:=0]
range(biom.dat[aoi>800 & isa.frac<0.90,live.MgC.ha.perv],  na.rm=T) ## 0 - 6786
# hist(biom.dat[aoi>800 & isa.frac<0.98,live.MgC.ha.perv]) ## a small number of very extreme values
hist(biom.dat[aoi>800 & live.MgC.ha.perv<284, live.MgC.ha.perv])

### figure out forest "age" for the cells (using coefficients for NE total)
## age based on ground area
biom.dat[,age.ground:=log(1-(live.MgC.ha.ground/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.ground>120, age.ground:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.ground), age.ground:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.ground:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.ground:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.ground:=NA]

## age based on canopy area
biom.dat[,age.forest:=log(1-(live.MgC.ha.forest/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.forest>120, age.forest:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.forest), age.forest:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.forest:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.forest:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.forest:=NA]

## age based on pervious area
biom.dat[,age.perv:=log(1-(live.MgC.ha.perv/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.perv>120, age.perv:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.perv), age.perv:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.perv:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.perv:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.perv:=NA]

## frequency distributions of different methods
par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(biom.dat$age.ground,  main="Forest age, unadjusted", xlab="Age, yrs")
hist(biom.dat$age.forest,  main="Forest age, canopy area adj.", xlab="Age, yrs")
hist(biom.dat$age.perv,  main="Forest age, pervious area adj.", xlab="Age, yrs")

### frequency dist for biomass/ha
hist(biom.dat[aoi>800 & live.MgC.ha.ground<284, live.MgC.ha.ground])
hist(biom.dat[aoi>800 & live.MgC.ha.forest<284, live.MgC.ha.forest])
hist(biom.dat[aoi>800 & live.MgC.ha.perv<284, live.MgC.ha.perv])

a=123.67
b=0.04

## figure out next annual increment possible for the "forest" average present in each cell, based on projected "age" and corrected for area
biom.dat[,npp.ann.ground:=((a*(1-(exp(-b*(age.ground+1))))^3)-(a*(1-(exp(-b*(age.ground))))^3))*(1E-4)*aoi] ## by ground area
biom.dat[,npp.ann.forest:=((a*(1-(exp(-b*(age.forest+1))))^3)-(a*(1-(exp(-b*(age.forest))))^3))*(1E-4)*(aoi*can.frac)] ## by canopy area
biom.dat[,npp.ann.perv:=((a*(1-(exp(-b*(age.perv+1))))^3)-(a*(1-(exp(-b*(age.perv))))^3))*(1E-4)*(aoi*(1-isa.frac))] ## by pervious area
biom.dat[bos.biom30m<10, npp.ann.ground:=0] ## manually set negligible biomass cells to 0
biom.dat[bos.biom30m<10, npp.ann.forest:=0]
biom.dat[bos.biom30m<10, npp.ann.perv:=0]
biom.dat[isa.frac>0.99, npp.ann.perv:=0]
biom.dat[can.frac<0.01, npp.ann.forest:=0]

### look at some plots 
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.forest)
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.perv) ## hard to say, up to 0.2 MgC/pix
par(mfrow=c(3,1))
hist(biom.dat$npp.ann.ground, main="NPP, raw area", xlab="MgC/pix/yr")
hist(biom.dat$npp.ann.forest, main="NPP, canopy area", xlab="MgC/pix/yr")
hist(biom.dat$npp.ann.perv, main="NPP, pervious area", xlab="MgC/pix/yr")

## correct to MgC/ha/yr rather than MgC/pix/yr
biom.dat[,npp.ann.ground.ha:=(npp.ann.ground/aoi)*1E4]
biom.dat[,npp.ann.forest.ha:=(npp.ann.forest/aoi)*1E4]
biom.dat[,npp.ann.perv.ha:=(npp.ann.perv/aoi)*1E4]

hist(biom.dat[aoi>800, npp.ann.ground.ha])
hist(biom.dat[aoi>800, npp.ann.forest.ha])
hist(biom.dat[aoi>800, npp.ann.perv.ha])

## now need these as kernel density plots
par(mfrow=c(1,3))
plot(density(biom.dat[aoi>800, npp.ann.ground.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,3.4))
plot(density(biom.dat[aoi>800, npp.ann.forest.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,3.4))
plot(density(biom.dat[aoi>800, npp.ann.perv.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,3.4))

plot(density(biom.dat[aoi>800 & live.MgC.ha.ground<284, live.MgC.ha.ground], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.06))
plot(density(biom.dat[aoi>800 & live.MgC.ha.forest<284, live.MgC.ha.forest], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.06))
plot(density(biom.dat[aoi>800 & live.MgC.ha.perv<284, live.MgC.ha.perv], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.06))
## should these be groomed to remove data with no forest biomass?

par(mfrow=c(1,3))
plot(density(biom.dat[aoi>800 & bos.biom30m>10, npp.ann.ground.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,2.2))
plot(density(biom.dat[aoi>800 & bos.biom30m>10, npp.ann.forest.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,2.2))
plot(density(biom.dat[aoi>800 & bos.biom30m>10, npp.ann.perv.ha], na.rm=T, adjust=2, bw=0.05), ylim=c(0,2.2))

plot(density(biom.dat[aoi>800 & live.MgC.ha.ground<284 & bos.biom30m>10, live.MgC.ha.ground], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.032))
plot(density(biom.dat[aoi>800 & live.MgC.ha.forest<284 & bos.biom30m>10, live.MgC.ha.forest], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.032))
plot(density(biom.dat[aoi>800 & live.MgC.ha.perv<284 & bos.biom30m>10, live.MgC.ha.perv], na.rm=T, adjust=2, bw=3.5), ylim=c(0, 0.032))


## ok work on some pretty density overlays in ggplot
library(ggplot2)
x <- data.frame(v1=rnorm(100),v2=rnorm(100,1,1),v3=rnorm(100,0,2))
library(ggplot2);library(reshape2)
data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.7)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)

groom.mgc <- as.data.frame(biom.dat[aoi>800 & live.MgC.ha.perv<=300 & bos.biom30m>10, 5:7, with=F])
mgc.melt <- melt(groom.mgc)
groom.npp <- as.data.frame(biom.dat[aoi>800 & live.MgC.ha.perv<=300 & bos.biom30m>10, 14:16, with=F])
npp.melt <- melt(groom.npp)
dim(mgc.melt)
ggplot(mgc.melt, aes(x=value)) + 
  geom_density(aes(group=variable, colour=variable))
ggplot(npp.melt, aes(x=value)) + 
  geom_density(aes(group=variable, colour=variable))

## let's add the growth curve over the density plot for biomass density
a=123.67
b=0.04
AGE=seq(0,300)
y = a*(1-exp(-b*AGE))^3
plot(AGE, y)
z = diff(y)
z.rel <- z/y[2:301]
length(y)
length(z) ## this is growth after 1 year, 2 years, etc
plot(y[2:301], z)
points(y[2:151], z.rel, pch=14, cex=0.5, col="red") #hyperbolic
plot(AGE[2:151], z) ## this is the gain curve over site age
growth.curve <- as.data.frame(cbind(y[2:301],z))

ggplot(mgc.melt, aes(x=value)) + 
  geom_density(aes(group=variable, colour=variable), alpha=0.12) +
  geom_line(data=growth.curve, aes(x=V1, y=z*0.01))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_colour_discrete(name="Area basis",
                        breaks=c("live.MgC.ha.ground", "live.MgC.ha.forest", "live.MgC.ha.perv"),
                        labels=c("ground", "forest", "pervious"))+
  labs(x = "Biomass MgC/ha")

  
ggplot(npp.melt, aes(x=value, fill=variable)) + 
  geom_density(aes(group=variable, colour=variable), alpha=0.12) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
    scale_colour_discrete(name="Area basis",
                          breaks=c("npp.ann.ground.ha", "npp.ann.forest.ha", "npp.ann.perv.ha"),
                          labels=c("ground", "forest", "pervious"))+
    labs(x = "NPP MgC/ha/yr")+
    guides(fill=FALSE)



### combined biomass.growth~dbh for Street, Andy, FIA
### individual FIA tree data from sites near Boston
dat <- read.csv("data/FIA/MA_Tree_Data_ID.csv")
dat <- as.data.table(dat)
names(dat)[1] <- c("TreeID")
names(dat)[2] <- c("PlotID")
spec <- read.csv("data/FIA/REF_SPECIES.csv")
dat <- merge(x=dat, y=spec[,c("SPCD", "GENUS")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
dat$GENUS <- as.character(dat$GENUS)
dat$GENUS <- as.factor(dat$GENUS)
table(dat$GENUS)
dat$GENUS.num <- as.numeric(dat$GENUS)

b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
dat$biom0 <- biom.pred(dat$DIAM_T0)
dat$biom1 <- biom.pred(dat$DIAM_T1)
dat$biom.delt <- dat$biom1-dat$biom0
dat$biom.rel <- (dat$biom.delt/4.8)/dat$biom0 ## annualized relative growth increment

### now the andy trees
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

### lagged biomass gains, relative to previous year biomass
biom.rel <- data.frame(c(2015:1990))
for(i in 1:dim(bop)[1]){
  te <- unlist(bop[i,6:32])
  te.di <- (-1*diff(te, lag=1))
  biom.rel[,i+1] <- te.di/te[-1]
}
names(biom.rel)[-1] <- paste0("Tree", bop[,incr.ID])
names(biom.rel)[1] <- "Year" ## this has all the relative biomass gains by year (row) for each tree with increment (column)

### figure the mean relative biomass gain per year across all years in each tree
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
andy.bai[,seg.Edge:=as.factor(seg.Edge)]

### now the street trees
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
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees

par(mfrow=c(1,1))
plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean), col=as.numeric(andy.bai$seg.Edge))
plot(street[record.good==1, log(dbh.2006)], street[record.good==1, log(npp.ann.rel)]) # r2 = 0.54, slope significant
plot(dat[,log(DIAM_T0)], dat[,log(biom.rel)])

## just do a plot first
plot(dat[,log(DIAM_T0)],  ### this is FIA data
     dat[,log(biom.rel)], 
     cex=0.2, col="gray64", pch=1, ylim=c(-9, 2), xlim=c(log(5), log(110)))
points(street[record.good==1, log(dbh.2006)],  ### this is street trees
       street[record.good==1, log(npp.ann.rel)], 
       pch=15, cex=0.4, col="salmon")
andy.cols <- c("royalblue", "skyblue")
points(andy.bai[,log(avg.dbh)], andy.bai[,log(growth.mean)], 
       col=andy.cols[as.numeric(andy.bai$seg.Edge)],
       pch=15, cex=0.6)

summary(street[record.good==1,npp.ann]) ## a tiny handful of street trees show zero growth

### lets look untransformed
plot(dat[,(DIAM_T0)], 
     dat[,(biom.rel)], 
     cex=0.2, col="gray64", pch=14)
points(street[record.good==1, (dbh.2006)], 
       street[record.good==1, (npp.ann.rel)], 
       pch=15, cex=0.4, col="salmon")
andy.cols <- c("royalblue", "skyblue")
points(andy.bai[,(avg.dbh)], andy.bai[,(growth.mean)], 
       col=andy.cols[as.numeric(andy.bai$seg.Edge)],
       pch=15, cex=0.6)



