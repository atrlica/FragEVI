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
library(data.table)
library(raster)
bos.forest <- raster("processed/boston/bos.forest_only.tif")
bos.dev <- raster("processed/boston/bos.dev_only.tif")
bos.hdres <- raster("processed/boston/bos.hdres_only.tif")
bos.ldres <- raster("processed/boston/bos.ldres_only.tif")
bos.lowveg <- raster("processed/boston/bos.lowveg_only.tif")
bos.water <- raster("processed/boston/bos.water_only.tif")

for.sum <- sum(getValues(bos.forest), na.rm=T)
dev.sum <- sum(getValues(bos.dev), na.rm=T)
hdres.sum <- sum(getValues(bos.hdres), na.rm=T)
ldres.sum <- sum(getValues(bos.ldres), na.rm=T)
lowveg.sum <- sum(getValues(bos.lowveg), na.rm=T)
water.sum <- sum(getValues(bos.water), na.rm=T)
bos.aoi <- raster("processed/boston/bos.aoi.tif")
aoi.sum <- sum(getValues(bos.aoi), na.rm=T)
for.sum/aoi.sum ## 8.2%
dev.sum/aoi.sum ## 38%
hdres.sum/aoi.sum ## 38.7%
ldres.sum/aoi.sum ## 2.0%
lowveg.sum/aoi.sum # 10.7%
water.sum/aoi.sum # 2.1%

### processing script to get cumulative area by edge distance within LULC
## sum of canopy area, masked by lulc
can.sum.ma <- function(x, m) { # x is canopy 0/1 1m raster object, m is mask (LULC 1/0 map)
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
# buff.dist <- as.integer(unlist(rm_between(can.buffs, "nocan_", "mbuff", extract=TRUE)))
buff.dist <- sub("bos.nocan_", "", can.buffs); buff.dist <- as.integer(sub("mbuff.tif", "", buff.dist))
can.master <- raster("processed/boston/bos_can01_filt.tif")
# can.tot <- sum(getValues(can.master), na.rm=T) ### 31.7% of boston raster is canopy
# aoi.master <- raster("processed/boston/bos.aoi.tif")
# aoi.tot <- sum(getValues(aoi.master), na.rm=T)

## make a data table of edge pixels by lulc
can.master <- as.data.table(as.data.frame(raster("processed/boston/bos_can01_filt.tif")))
can.10mbuff <- as.data.table(as.data.frame(raster("processed/boston/bos.ed10.tif")))
lulc <- as.data.table(as.data.frame(raster("processed/boston/bos.lulc.lumped.tif")))
fat <- fwrite(cbind((can.master),
             (can.10mbuff),
             (lulc)), "processed/boston/bos.lulc_candist.csv")
write.csv(fat, "processed/boston/bos.lulc_candist.csv")
## cumulative canopy area by edge distance, for each lulc class
lu.classes <- c("forest", "dev", "hdres", "ldres", "lowveg", "water")
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


################
### combined plot, canopy edge area as fraction of total area (by LULC class)
results <- read.csv("processed/bos.can.cummdist.csv")
lu.classes <- c("forest", "dev", "hdres", "ldres", "lowveg", "water")
cols=rainbow(length(lu.classes))
cols=c("forestgreen", "gray55", "red", "yellow", "green", "lightblue")
par(mar=c(4.5, 5.5, 1, 1), oma=c(0,0,0, 0), xpd=F)

### separated cumulative area curves for each LULC
# plot(results$dist, results$frac.tot.area, pch=1, col="black", type="l", lwd=7, bty="n", lty=1, 
#      xlab="Distance from edge (m)", ylab="Cummulative area fraction",
#      ylim=c(0, 0.9), xlim=c(0, 60), yaxt="n", font.lab=2, cex.lab=2, cex.axis=2)
# axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8), labels=c("0", "20%", "40%", "60%", "80%"), cex.axis=2)
# for(d in 1:length(lu.classes)){
#   dat <- read.csv(paste("processed/bos.can.cummdist.", lu.classes[d], ".csv", sep=""))
#   lines(dat$dist, dat$frac.tot.area, col=cols[d], type="l", lwd=3)
# }
# legend("right", x=40, y=0.80, cex=1, legend=c("Boston", "Forest", "Developed", "HDRes", "LDRes", "LowVeg"), fill=c("black", cols), bty="n")


### Edge distance cummulative by LULC stack
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
contain <- as.data.table(contain)
# m <- contain[, sum(ha), by=distance]
# sum(m$V1) #3952 ha
# results[results$dist==100, "pix.less.than"]/1E4 ### 3942 ha, all LULC
# ### ok this does seem to add up to the whole area canopy cover
# plot(contain[LULC=="hdres", ha], col="salmon")
# points(contain[LULC=="forest", ha], col="forestgreen")
# points(contain[LULC=="dev", ha], col="gray55")
# points(contain[LULC=="ldres", ha], col="gold")

library(ggplot2)
contain <- as.data.table(contain)
ggplot(contain[distance!=0,], aes(x=distance, y=ha, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  xlim(1,35)+
  scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
                    name="Land Cover",
                    breaks=c("forest", "dev", "hdres", "ldres", "lowveg", "water"),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Canopy distance from edge (m)")+
  ylab("Total canopy area (ha)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = c(0.8, 0.4))

### new hotness: matching viridis color scheme
# lulc.pal <- inferno(6)
library(viridis)
lulc.pal <- viridis(6)
# lulc.pal <- inferno(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])
### fix artifact of naming "distance from edge"
contain[,distance:=distance-1]
contain <- contain[distance>=0,]

png(filename="images/Can_dist_VIRIDIS_8x8in.png", width=8, height=8, units = "in", bg="white", res=600)
par(mfrow=c(1,1), mar=c(4,4,0.5,0.5))
ggplot(contain, aes(x=distance, y=ha, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  xlim(0,35)+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c("forest", "dev", "hdres", "ldres", "lowveg", "water"),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Canopy distance from edge (m)")+
  ylab("Total canopy area (ha)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.7, 0.5),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24, face="bold"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20, face="bold"))
dev.off()

png(filename="images/Can_dist_VIRIDIS_3x3in.png", width=3, height=3, units = "in", bg="white", res=600)
par(mfrow=c(1,1), mar=c(4,4,0.5,0.5))
ggplot(contain, aes(x=distance, y=ha, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  xlim(0,35)+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c("forest", "dev", "hdres", "ldres", "lowveg", "water"),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Canopy distance from edge (m)")+
  ylab("Total canopy area (ha)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.7, 0.5),
        axis.text=element_text(size=10),
        axis.title=element_text(size=11, face="bold"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10, face="bold"))
dev.off()


### tables for example "pixel" summary cover
library(rgeos)
target <- readOGR("processed/boston/fenway_sample1.shp")
inset <- readOGR("processed/boston/fenway_inset2.shp")
cov <- raster("processed/boston/bos.cov.V2.tif")
cov <- crop(cov, inset)
plot(cov); plot(target, add=T)
biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
biom <- crop(biom, inset)
plot(biom); plot(target, add=T)
ed <- raster("processed/boston/bos.ed10.tif")
ed <- crop(ed, inset)
plot(ed); plot(target, add=T)
ndvi <- raster("processed/boston/bos.ndvi.tif")
ndvi <- crop(ndvi, inset)
plot(ndvi); plot(target, add=T)

d <- unlist(extract(cov, target))
sum(d%in%c(6,5))/length(d) ## canopy fraction
sum(d%in%c(5,4))/length(d) ## impervious fraction
e <- unlist(extract(ed, target))
sum(e)/sum(d%in%c(6,5)) ## canopy edge fraction
b <- unlist(extract(biom,target))
sum(b)/(2000) ### total biomass MgC
(sum(b)/(2000))/(length(b)/1E4) ### MgC/ha-ground
(sum(b)/(2000))/((length(b)*(sum(d%in%c(6,5))/length(d)))/1E4) ### MgC/ha-canopy
(sum(b)/(2000))/((length(b)*(1-(sum(d%in%c(5,4))/length(d))))/1E4) ### MgC/ha-perv
n <- unlist(extract(ndvi, target))
mean(n) ## mean pixel ndvi

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




#######
##
## FIGURE 2A: STEM LEVEL BIOMASS GROWTH (FIA, ANDY, STREET)
#######
library(MASS)
library(ggplot2)
library(viridis)
get_density <- function(x, y, n = 100) { ## two dimensional kernel density estimation from MASS package
  dens <- MASS::kde2d(x = x, y = y, n = n) 
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

### scale tranformation function for -1 to 1 log scaling
library(scales)
stemg.xform <- function(x){
  log(x+0.5)
}
I.stemg.xform <- function(x){
  exp(x)-0.5
}
stemg_trans <- function(){trans_new(name="stemg",
                                    transform=stemg.xform,
                                    inverse=I.stemg.xform)
}

# stemg.xform(5)
# I.stemg.xform(5)
# I.stemg.xform(stemg.xform(5)) ## it works

## set up common visual parameters
xlim.all <- c(4, 100)
ylim.all <- c(-0.1, 1)
xlim.all.log <- c(log(4), log(100))
ylim.all.log <- c(log(0.002), log(1))
title.size <- 12
axis.marks <- 9
axis.titles <- 10
pt.size <- 0.7
theme.master <-   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.title.x = element_text(face="bold", size=axis.titles),
                        axis.title.y = element_text(face="bold", size=axis.titles),
                        axis.text.x = element_text(face="plain", size=axis.marks),
                        axis.text.y = element_text(face="plain", size=axis.marks),
                        plot.title = element_text(face="bold", size=title.size))
med.lines.col <- "gray75"
med.lines.width <- 0.4
fit.col <- "gray65"
# fit.col <- "royalblue2"
fit.width <- 0.8
fit.type <- "4121"
legend.title.size=9
legend.text.size=8
alpha.master <- 0.2

####
### FIA rural trees
live <- read.csv("data/FIA/MA_Tree_Data_ID_NOMORT.csv")
live <- as.data.table(live)
names(live)[1] <- c("TreeID")
names(live)[2] <- c("PlotID")
spec <- read.csv("data/FIA/REF_SPECIES.csv")
live <- merge(x=live, y=spec[,c("SPCD", "GENUS", "SPECIES")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
live$GENUS <- as.character(live$GENUS)
live$GENUS <- as.factor(live$GENUS)
live$GENUS.num <- as.numeric(live$GENUS)
spp.allo <- read.csv("data/FIA/spp_allometrics.csv") ## manually entered selected map from spp to b0+b1
live[,spp:=paste(substr(GENUS, 1,1), ".", SPECIES, sep="")]
live <- merge(x=live, y=spp.allo[,c("spp", "b0", "b1")], by="spp", all.x=T)
live[is.na(b0), b0:=(-2.48)]
live[is.na(b1), b1:=2.4835]
biom.pred2 <- function(b0, b1, x){exp(b0+(b1*log(x)))}
## class as hard or soft wood
live[,type:="H"]
live[spp%in%c("P.strobus", "P.resinosa", "T.canadensis", "A.balsamea"), type:="S"]
live[,type:=as.factor(type)]
live[,biom0.spp:=biom.pred2(b0, b1, DIAM_T0)]
live[,biom1.spp:=biom.pred2(b0, b1, DIAM_T1)]
live[,growth.ann:=(biom1.spp-biom0.spp)/4.8]
live[,growth.ann.rel:=growth.ann/biom0.spp]
# range(live$npp.ann.rel, na.rm=T) #-12% to 52%
mod.fia.nls <- nls(growth.ann.rel ~ exp(a + b * log(DIAM_T0)), data=live, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
mm <- summary(mod.fia.nls) ## all factors significant

pred_live <- data.frame(growth.pred=predict(mod.fia.nls, newdata=data.frame(DIAM_T0=seq(live[,min(DIAM_T0, na.rm=T)],
                                                                                     100, by=0.2))),
                        DIAM_T0=seq(live[,min(DIAM_T0, na.rm=T)],
                                    100, by=0.2))

live$density <- get_density(live$DIAM_T0, live$growth.ann.rel)
live[is.na(density), density:=0]

## single color, density via alpha density
# fia.mono <- ggplot(live, aes(DIAM_T0, growth.ann.rel))+
#   geom_point(alpha=alpha.master, color=plasma(6)[1], size=pt.size)+
#   lims(x=xlim.all, y=ylim.all)+
#   # geom_vline(xintercept=live[,median(DIAM_T0, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   # geom_hline(yintercept=live[,median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_live, aes(x=DIAM_T0, y=growth.pred), color=fit.col, linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="FIA.linear")+
#   theme.master


## log+a transformed
fia.mono.log <- ggplot(live, aes(DIAM_T0, growth.ann.rel))+
  geom_point(alpha=alpha.master, color=plasma(6)[1], size=pt.size)+
  scale_y_continuous(trans="stemg", breaks=c(-0.1, 0, 0.2, 0.5, 1.0), limits=ylim.all)+
  scale_x_continuous(trans="log", limits=xlim.all, breaks=c(5, 10, 30, 60, 100))+
  # geom_vline(xintercept=live[,median(DIAM_T0, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  # geom_hline(yintercept=live[,median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  geom_line(data = pred_live, aes(x=DIAM_T0, y=growth.pred), color=fit.col, linetype=fit.type, size=fit.width)+
  labs(x = "Stem DBH (cm)", y="Stem Growth (kg/kg)", title="Rural Forest '08-'16")+
  theme.master

# ## asinh transform
# library(scales)
# asinh_trans <- function(){
#   trans_new(name = 'asinh', transform= function(x) asinh(x),
#             inverse=function(x) sinh(x))
# }
# 
# ggplot(live, aes((DIAM_T0), (growth.ann.rel)))+
#   geom_point(alpha=alpha.master, color=plasma(6)[1], size=pt.size)+
#   scale_y_continuous(trans='asinh', breaks=c(-0.2, 0, 0.2, 0.5, 1.0), limits=ylim.all)+
#   lims(x=xlim.all)+
#   # geom_vline(xintercept=live[,median(DIAM_T0, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   # geom_hline(yintercept=live[,median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   # geom_line(data = pred_live, aes(x=DIAM_T0, y=growth.pred), color=fit.col, linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="FIA.asinh")+
#   theme.master

# ### density via color gradient
# fia.dens <- ggplot(live, aes(DIAM_T0, growth.ann.rel, color=density))+
#   geom_point(alpha=alpha.master, size=pt.size)+
#   scale_colour_viridis(option="B", guide=FALSE)+
#   lims(x=xlim.all, y=ylim.all)+
#   geom_vline(xintercept=live[,median(DIAM_T0, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_hline(yintercept=live[,median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_live, aes(x=DIAM_T0, y=growth.pred), color=fit.col, linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="FIA")+
#   theme.master


####
### Andy Trees
andy.bai <- as.data.table(read.csv("processed/andy.bai.dbh.pseudo.csv"))

### log-transformed growth~dbh*edge
mod.andy.log <- lm(log(biom.rel.ann)~log(dbh.start)+seg.Edge, data=andy.bai[dbh.start>=5,]) #R2 0.31, dbh:edge not significant
pred_andy.edge <- data.frame(growth.pred=exp(predict(mod.andy.log, newdata=data.frame(dbh.start=seq(from=andy.bai[dbh.start>=5,min(dbh.start, na.rm=T)], to=100, length.out=200),
                                                                                  seg.Edge=rep("E", 200)))),
                          dbh.start=seq(from=andy.bai[dbh.start>=5,min(dbh.start, na.rm=T)],
                                       to=100, length.out=200))
pred_andy.int <- data.frame(growth.pred=exp(predict(mod.andy.log, newdata=data.frame(dbh.start=seq(from=andy.bai[dbh.start>=5,min(dbh.start, na.rm=T)], to=100, length.out=200),
                                                                                      seg.Edge=rep("I", 200)))),
                             dbh.start=seq(from=andy.bai[dbh.start>=5, min(dbh.start, na.rm=T)],
                                           to=100, length.out=200))
andy.bai[dbh.start>=5, density:=get_density(andy.bai[dbh.start>=5, dbh.start], andy.bai[dbh.start>=5, biom.rel.ann])]
andy.bai[is.na(density), density:=0.01]

## split colors (seg.Edge), density via alpha
andy.col <- plasma(6)[c(3,5)]
names(andy.col) <- levels(andy.bai$seg.Edge)
col.map <- scale_color_manual(name="Tree position", values=andy.col, labels=c("Edge <10m", "Interior"))
shape.map <- scale_shape_manual(name="Tree position", values=c(16,17), labels=c("Edge <10m", "Interior"))

# andy.mono <- ggplot(andy.bai[dbh.start>=5], aes(dbh.start, biom.rel.ann, colour=seg.Edge))+
#   geom_point(aes(shape=seg.Edge), alpha=alpha.master*1.6, size=pt.size)+
#   col.map+
#   shape.map+
#   lims(x=xlim.all, y=ylim.all)+
#   geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(dbh.start, na.rm=T)], color="gray55", size=med.lines.width)+
#   geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(biom.rel.ann, na.rm=T)], color="gray55", size=med.lines.width)+
#   geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(dbh.start, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(biom.rel.ann, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_andy.edge, aes(x=dbh.start, y=growth.pred), color="darkmagenta", linetype=fit.type, size=fit.width)+
#   geom_line(data = pred_andy.int, aes(x=dbh.start, y=growth.pred), color="chocolate3", linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="Urban Forest '90-'16")+
#   theme.master+
#   theme(legend.position = c(0.5, 0.65),
#         legend.title = element_text(size=legend.title.size, face="bold"),
#         legend.background = element_rect(fill = "white"),
#         legend.key =  element_blank(),
#         legend.key.size = unit(0.6, "lines"),
#         legend.text = element_text(size=legend.title.size-1, face="plain"),
#         legend.justification = "center")+
#   guides(shape = guide_legend(override.aes = list(size=pt.size*1.7,
#                                                   alpha=0.8,
#                                                   color=c("darkmagenta",
#                                                           "chocolate3"))))

### log+a transformed
  andy.mono.log <- ggplot(andy.bai[dbh.start>=5], aes(dbh.start, biom.rel.ann, colour=seg.Edge))+
  geom_point(aes(shape=seg.Edge), alpha=alpha.master*1.6, size=pt.size)+
  col.map+
  shape.map+
  scale_y_continuous(trans="stemg", breaks=c(-0.1, 0, 0.2, 0.5, 1.0), limits=ylim.all)+
  scale_x_continuous(trans="log", limits=xlim.all, breaks=c(5, 10, 30, 60, 100))+
  # geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(dbh.start, na.rm=T)], color="gray55", size=med.lines.width)+
  # geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(biom.rel.ann, na.rm=T)], color="gray55", size=med.lines.width)+
  # geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(dbh.start, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  # geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(biom.rel.ann, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  geom_line(data = pred_andy.edge, aes(x=dbh.start, y=growth.pred), color="darkmagenta", linetype=fit.type, size=fit.width)+
  geom_line(data = pred_andy.int, aes(x=dbh.start, y=growth.pred), color="chocolate3", linetype=fit.type, size=fit.width)+
  labs(x = "Stem DBH (cm)", y="Stem Growth (kg/kg)", title="Urban Forest '90-'16")+
  theme.master+
  theme(legend.position = c(0.55, 0.65),
        legend.title = element_text(size=legend.title.size, face="bold"),
        legend.background = element_rect(fill = "white"),
        legend.key =  element_blank(),
        legend.key.size = unit(0.6, "lines"),
        legend.text = element_text(size=legend.title.size-1, face="plain"),
        legend.justification = "center")+
  guides(shape = guide_legend(override.aes = list(size=pt.size*1.7,
                                                  alpha=0.8,
                                                  color=c("darkmagenta",
                                                          "chocolate3"))))


# ## density via color shading
# shape.map <- scale_shape_manual(name="Tree position", values=c(16,17), labels=c("Edge <10m", "Interior"))
# 
# andy.dens <- ggplot(andy.bai[dbh.start>=5], aes(dbh.start, biom.rel.ann, color=density))+
#   geom_point(aes(shape=seg.Edge), alpha=alpha.master*1.6, size=pt.size)+
#   shape.map+
#   scale_colour_viridis(option="B", guide=FALSE)+
#   lims(x=xlim.all, y=ylim.all)+
#   geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(dbh.start, na.rm=T)], color="gray55", size=med.lines.width)+
#   geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="E",median(biom.rel.ann, na.rm=T)], color="gray55", size=med.lines.width)+
#   geom_vline(xintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(dbh.start, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_hline(yintercept=andy.bai[dbh.start>=5 & seg.Edge=="I",median(biom.rel.ann, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_andy.edge, aes(x=dbh.start, y=growth.pred), color="darkmagenta", linetype=fit.type, size=fit.width)+
#   geom_line(data = pred_andy.int, aes(x=dbh.start, y=growth.pred), color="chocolate3", linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="Urban Forest '90-'15")+
#   theme.master+  
#   theme(legend.position = c(0.6, 0.5),
#         legend.title = element_text(size=legend.title.size, face="bold"),
#         legend.background = element_rect(fill = "white"),
#         legend.key =  element_blank(),
#         legend.key.size = unit(0.6, "lines"),
#         legend.text = element_text(size=legend.title.size-1, face="plain"),
#         legend.justification = "center")+
#   guides(shape = guide_legend(override.aes = list(size=pt.size*1.7,
#                                                   alpha=0.8,
#                                                   color=c("darkmagenta",
#                                                           "chocolate3"))))

####
###  street trees
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street <- as.data.table(street)
street[, record.good:=0] ## ID records with good data quality
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, growth.ann:=(biom.2014-biom.2006)/8]
street[record.good==1, growth.ann.rel:=growth.ann/biom.2006]
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees

mod.street.nls <- nls(growth.ann.rel ~ exp(a + b * log(dbh.2006)), data=street[record.good==1,], start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
mm <- summary(mod.street.nls) ## RSE is pretty high 0.33

### note here that though log-transform of dbh.2006 is not specified, predict() is generating predictions consistent with log-tranforming prior to inserting
pred_street <- data.frame(growth.pred=predict(mod.street.nls, newdata=data.frame(dbh.2006=seq(street[,min(dbh.2006, na.rm=T)],
                                                                                        100, by=0.2))),
                        dbh.2006=seq(street[,min(dbh.2006, na.rm=T)],
                                    100, by=0.2))
street[record.good==1, density:=get_density(street[record.good==1, dbh.2006], street[record.good==1, growth.ann.rel])]

# ## single color, density via alpha
# street.mono <- ggplot(street[record.good==1], aes(dbh.2006, growth.ann.rel))+
#   geom_point(alpha=alpha.master, color=plasma(6)[4], size=pt.size)+
#   lims(x=xlim.all, y=ylim.all)+
#   geom_vline(xintercept=street[record.good==1, median(dbh.2006, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_hline(yintercept=street[record.good==1, median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_street, aes(x=dbh.2006, y=growth.pred), color="gray40", linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="Street trees '06-'14")+
#   theme.master


## log x+a transform
street.mono.log <- ggplot(street[record.good==1], aes(dbh.2006, growth.ann.rel))+
  geom_point(alpha=alpha.master, color=plasma(6)[4], size=pt.size)+
  scale_y_continuous(trans="stemg", breaks=c(-0.1, 0, 0.2, 0.5, 1.0), limits=ylim.all)+
  scale_x_continuous(trans="log", limits=xlim.all, breaks=c(5, 10, 30, 60, 100))+
  # geom_vline(xintercept=street[record.good==1, median(dbh.2006, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  # geom_hline(yintercept=street[record.good==1, median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
  geom_line(data = pred_street, aes(x=dbh.2006, y=growth.pred), color="gray40", linetype=fit.type, size=fit.width)+
  labs(x = "Stem DBH (cm)", y="Stem Growth (kg/kg)", title="Street Trees '06-'14")+
  theme.master

### test -- does log transforming either scale differently fuck something up with the way the fit line is shown?
# plot(street[record.good==1, dbh.2006], street[record.good==1, growth.ann.rel])
# test <- seq(5, 100, by=1)
# 
# ## standard log transform on both axes
# plot(street[record.good==1, log(dbh.2006)], street[record.good==1, log(growth.ann.rel)])
# points(log(test), log((exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))))), pch=15, col="blue") ## linear fit
# 
# ### with two slightly different log transforms on either scale
# plot(street[record.good==1, log(dbh.2006)], street[record.good==1, log(growth.ann.rel+1)], ylim=c(-0.2, 1))
# points(log(test), log((exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))))+1), pch=15, col="blue") ## a bent linear line of best fit
# points(log(test), log((exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))))+0), pch=15, col="red") ## barely in the same range

# ## density via color shading
# street.dens <- ggplot(street[record.good==1], aes(dbh.2006, growth.ann.rel, color=density))+
#   geom_point(alpha=alpha.master, size=pt.size)+
#   scale_colour_viridis(option="B", guide=FALSE)+
#   lims(x=xlim.all, y=ylim.all)+
#   geom_vline(xintercept=street[record.good==1, median(dbh.2006, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_hline(yintercept=street[record.good==1, median(growth.ann.rel, na.rm=T)], color=med.lines.col, size=med.lines.width)+
#   geom_line(data = pred_street, aes(x=dbh.2006, y=growth.pred), color="black", linetype=fit.type, size=fit.width)+
#   labs(x = "Stem DBH (cm)", y="Relative Growth (kg/kg)", title="Street trees '06-'14")+
#   theme.master

# ### collate monocolored plots
# library(gridExtra)
# png(width=8, height=3, units="in", res=600, bg="white", filename="images/Fig2A_stem-dbh-growth_mono.png")
# grid.arrange(grobs=list(fia.mono, andy.mono, street.mono),
#              widths=c(1,1,1))
# dev.off()

### collate monocolored log-transformed plots
library(gridExtra)
png(width=8, height=3, units="in", res=600, bg="white", filename="images/Fig2A_stem-dbh-growth_mono_logxform.png")
grid.arrange(grobs=list(fia.mono.log, andy.mono.log, street.mono.log),
             widths=c(1,1,1))
dev.off()

# ### collate density colored plots
# png(width=8, height=3, units="in", res=600, bg="white", filename="images/Fig2A_stem-dbh-growth_dens.png")
# grid.arrange(grobs=list(fia.dens, andy.dens, street.dens),
#              widths=c(1,1,1))
# dev.off()

### draft plots, cumulative NPP by density
############
### cumulative npp~biomass graph w/ LULC stacked --- hybrid results
npp.dat <- as.data.table(read.csv("processed/npp.estimates.V3.csv"))
npp.dat[,biom.Mg.ha:=((bos.biom30m/1000)/(aoi*can.frac))*1E4] ### Forest basis
npp.dat[can.frac<0.005, biom.Mg.ha:=0]
# range(npp.dat[aoi>800, biom.Mg.ha], na.rm=T) ## up to 570 Mg/ha, canopy basis
# hist(npp.dat[aoi>800 & lulc==2,biom.Mg.ha])
npp.dat[,biom.bin:=cut(biom.Mg.ha, seq(0, 570, by = 10), right=FALSE, ordered_result=T)]
contain.npp <- npp.dat[!is.na(bos.biom30m) & !is.na(lulc), .(sum(hyb.npp.mod, na.rm=T), 
                                                             sum(fia.empir.npp.kg.hw.forest, na.rm=T),
                                                             sum(bos.biom30m, na.rm=T)), by=.(biom.bin, lulc)]
contain.npp <- contain.npp[order(biom.bin), ]
contain.npp[, bin.num:=as.numeric(biom.bin)]
names(contain.npp)[3:5] <- c("hyb.npp.kg", "fia.npp.kg", "tot.biom.kg")

## where is the biomass compared to where is the npp?
npp.tot <- npp.dat[, sum(hyb.npp.mod, na.rm=T)]
ar <- npp.dat[, sum(hyb.npp.mod, na.rm=T), by=lulc]
ar$npp.frac <- ar$V1/npp.tot #20% forest, 16% dev, 50% hdres

biom.tot <- npp.dat[, sum(bos.biom30m, na.rm=T)]
af <- npp.dat[, sum(bos.biom30m, na.rm=T), by=lulc]
af$biom.frac <- af$V1/biom.tot ## 13% dev, 31% forest, 43% hdres

## plots of biomass and NPP density by LULC along binned biomass density
par(mfrow=c(1,3), mar=c(3,4,4,1))
col.tmp <- c("forestgreen", "gray55", "salmon", "gold", "chartreuse3")
### summed biomass by biomass density bin and LULC
plot(contain.npp[lulc==1, bin.num], contain.npp[lulc==1, tot.biom.kg/1E6], 
     col=col.tmp[1], pch=15, type="l", lwd=3, main="Total Biomass",
     ylim=c(0,18), xaxt="n", ylab="Total Biomass (Mg x 1000)") 
axis(side = 1, at = c(0, 10, 20, 30, 40, 50), labels = seq(0, 500, by=100))
for(e in 2:5){
  lines(contain.npp[lulc==e, bin.num], contain.npp[lulc==e, tot.biom.kg/1E6], col=col.tmp[e], pch=15, type="l", lwd=3)
}
# legend(x = 40, y = 15, legend = c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other veg."), fill=col.tmp, bty="n")
# mtext(side=1, "Biomass density, canopy basis (Mg-biomass/ha)", line=2)

### summed hyb.NPP by biom.density bin and LULC
col.tmp <- c("forestgreen", "gray55", "salmon", "gold", "chartreuse3")
plot(contain.npp[lulc==1, bin.num], contain.npp[lulc==1, hyb.npp.kg/2000], 
     col=col.tmp[1], pch=15, type="l", lwd=3, main="NPP, Urban Hybrid",
     ylim=c(0, 420), xaxt="n", ylab="Total NPP (MgC/yr)") 
axis(side = 1, at = c(0, 10, 20, 30, 40, 50), labels = seq(0, 500, by=100))
for(e in 2:5){
  lines(contain.npp[lulc==e, bin.num], contain.npp[lulc==e, hyb.npp.kg/2000], col=col.tmp[e], pch=15, type="l", lwd=3)
}
# legend(x = 40, y = 400, legend = c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other veg."), fill=col.tmp, bty="n")
mtext(side=1, "Biomass density, canopy basis (Mg-biomass/ha)", line=2)

### summed fia.empir.NPP by biom.density bin and LULC
col.tmp <- c("forestgreen", "gray55", "salmon", "gold", "chartreuse3")
plot(contain.npp[lulc==1, bin.num], contain.npp[lulc==1, fia.npp.kg/2000], 
     col=col.tmp[1], pch=15, type="l", lwd=3, main="NPP, FIA",
     ylim=c(0, 420), xaxt="n", ylab="Total NPP (MgC/yr)") 
axis(side = 1, at = c(0, 10, 20, 30, 40, 50), labels = seq(0, 500, by=100))
for(e in 2:5){
  lines(contain.npp[lulc==e, bin.num], contain.npp[lulc==e, fia.npp.kg/2000], col=col.tmp[e], pch=15, type="l", lwd=3)
}
legend(x = 20, y = 400, legend = c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other veg."), fill=col.tmp, bty="n")
# mtext(side=1, "Biomass density, canopy basis (Mg-biomass/ha)", line=2)

### cumulative with ascending binned biomass density
dog <- contain.npp
yup <- data.frame()
ruler <- data.frame(bin.num=seq(1,104))
for(f in 1:6){
  a <- contain.npp[lulc==f, ]
  a <- merge(a, ruler, by="bin.num", all=T)
  a[is.na(lulc), lulc:=f]
  a[is.na(tot.npp.kg), tot.npp.kg:=0]
  a[,cum.sum:=cumsum(tot.npp.kg)]
  yup <- rbind(yup, a)
}
yup[,LULC:=as.factor(lulc)]

library(ggplot2)
ggplot(yup, aes(x=bin.num, y=cum.sum/1000, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Biomass density (Mg biomass/ha)")+
  ylab("Cumulative NPP (Mg biomass)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = c(60, 1.5E04))



### cumulative npp~biomass graph w/ LULC stacked --- FIA empirical (canopy) results
contain.npp <- npp.dat[!is.na(bos.biom30m) & !is.na(lulc), sum(fia.empir.npp.kg.hw.forest, na.rm=T), by=.(biom.bin, lulc)]
contain.npp <- contain.npp[order(biom.bin), ]
contain.npp[, bin.num:=as.numeric(biom.bin)]
names(contain.npp)[3] <- "tot.npp.kg"

### summed NPP by biom.density bin and LULC
col.tmp <- c("forestgreen", "gray55", "salmon", "gold", "chartreuse3")
plot(contain.npp[lulc==1, bin.num], contain.npp[lulc==1, tot.npp.kg/1000], 
     col=col.tmp[1], pch=15, type="l", lwd=3,
     ylim=c(0, 500), xlab="Biomass density (Mg/ha)", xaxt="n", ylab="Total NPP (Mg-biomass)") 
axis(side = 1, at = seq(0,100, by=20), labels = c(0, 10000, 20000, 30000, 40000, 50000))
for(e in 2:5){
  lines(contain.npp[lulc==e, bin.num], contain.npp[lulc==e, tot.npp.kg/1000], col=col.tmp[e], pch=15, type="l", lwd=3)
}
legend(x = 60, y = 400, legend = c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other veg."), fill=col.tmp, bty="n")

### cumulative with ascending binned biomass density
dog <- contain.npp
yup <- data.frame()
ruler <- data.frame(bin.num=seq(1,104))
for(f in 1:6){
  a <- contain.npp[lulc==f, ]
  a <- merge(a, ruler, by="bin.num", all=T)
  a[is.na(lulc), lulc:=f]
  a[is.na(tot.npp.kg), tot.npp.kg:=0]
  a[,cum.sum:=cumsum(tot.npp.kg)]
  yup <- rbind(yup, a)
}
yup[,LULC:=as.factor(lulc)]

library(ggplot2)
ggplot(yup, aes(x=bin.num, y=cum.sum/1000, fill=LULC)) + 
  geom_area(aes(fill=LULC))+
  scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Biomass density (Mg/ha)")+
  ylab("Total NPP (kg biomass)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = c(60, 1.5E04))

#####
### FIGURE 4: NPP VS BIOMASS DENSITY
########
#### SUMMED NPP, BINNED BY DENSITY AND COLORED BY LULC
library(reshape2)
npp.dat <- read.csv("processed/npp.estimates.V3.csv")
npp.dat <- as.data.table(npp.dat)
npp.dat[,biom.MgC.ha.ground:=((bos.biom30m/2000)/aoi)*1E4] ## per-pixel ground biomass density
npp.dat[,biom.MgC.ha.can:=((bos.biom30m/2000)/(aoi*can.frac))*1E4] ## per-pixel canopy biomass density
npp.dat[can.frac<0.005, biom.MgC.ha.can:=0]
npp.dat[,fia.empir.npp.MgC.ha.ground:=((fia.empir.npp.kg.hw.forest/2000)/aoi)*1E4] 
npp.dat[,fia.empir.npp.MgC.ha.can:=((fia.empir.npp.kg.hw.forest/2000)/(aoi*can.frac))*1E4]
npp.dat[can.frac<0.005, fia.empir.npp.MgC.ha.can:=0]
npp.dat[,hyb.npp.MgC.ha.ground:=((hyb.npp.mod.cap/2000)/aoi)*1E4]
npp.dat[,hyb.npp.MgC.ha.can:=((hyb.npp.mod.cap/2000)/(aoi*can.frac))*1E4]
npp.dat[can.frac<0.005, hyb.npp.MgC.ha.can:=0]
npp.dat2 <- npp.dat[aoi>800,]
# range(npp.dat[aoi>800, biom.Mg.ha], na.rm=T) ## up to 570 Mg/ha, canopy basis
hist(npp.dat2[, biom.MgC.ha.can])
hist(npp.dat2[, biom.MgC.ha.ground])
npp.dat2[,biom.bin:=cut(biom.MgC.ha.can, seq(0, 300, by = 10), right=FALSE, ordered_result=T)]
npp.dat2[,biom.bin.ground:=cut(biom.MgC.ha.ground, seq(0, 300, by = 10), right=FALSE, ordered_result=T)]

###
### plots ordered in bins along MgC/ha-CANOPY
contain.npp <- npp.dat2[!is.na(bos.biom30m) & !is.na(lulc), .(sum(hyb.npp.mod.cap, na.rm=T)/2000, 
                                                                       sum(fia.empir.npp.kg.hw.forest, na.rm=T)/2000,
                                                                       sum(bos.biom30m, na.rm=T)/2000), by=.(biom.bin, lulc)]
contain.npp <- contain.npp[order(biom.bin), ]
contain.npp[, bin.num:=as.numeric(biom.bin)]
names(contain.npp)[3:5] <- c("hyb.npp.MgC", "fia.npp.MgC", "tot.biom.MgC")
contain.npp[,LULC:=as.factor(lulc)]

### for viridis coloring
lulc.pal <- viridis(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])

hyb.npp.plot <- ggplot(contain.npp, aes(x=bin.num, y=hyb.npp.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("MgC/ha-canopy")+
  ylab("Hybrid, total NPP (MgC/yr)")+
  theme(axis.line = element_line(colour = "black"),
        axis.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
  ylim(0, 1400)
hyb.npp.plot

fia.npp.plot <- ggplot(contain.npp, aes(x=bin.num, y=fia.npp.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("MgC/ha-canopy")+
  ylab("FIA, total NPP (MgC/yr)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face="bold"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position=c(0.8, 0.6),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
  ylim(0, 1400)
fia.npp.plot

biom.plot <- ggplot(contain.npp, aes(x=bin.num, y=tot.biom.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("MgC/ha-canopy")+
  ylab("Total biomass (MgC x 1000)")+
  theme(axis.line = element_line(colour = "black"),
        axis.title = element_text(face="bold"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
  scale_y_continuous(breaks=c(0, 10000, 20000, 30000), labels=c(0, 10, 20, 30))
biom.plot

### now combine
library(gridExtra)
png(width=8, height=3, units="in", res=600, bg="white", filename="images/Fig4_biom-dens_NPP_comb.png")
grid.arrange(grobs=list(biom.plot, hyb.npp.plot, fia.npp.plot),
             widths=c(1,1,1))
dev.off()


### using the older manual color scheme
#### Biomass density in MgC/ha-canopy
### HYBRID STACKED NPP
# ggplot(contain.npp, aes(x=bin.num, y=hyb.npp.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-canopy)")+
#   ylab("Hybrid, total NPP (MgC/yr)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
#   ylim(0, 1400)
### FIA STACKED NPP
# ggplot(contain.npp, aes(x=bin.num, y=fia.npp.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-canopy)")+
#   ylab("FIA, total NPP (MgC/yr)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
#   ylim(0, 1400)
### TOTAL BIOMASS STACKED
# ggplot(contain.npp, aes(x=bin.num, y=tot.biom.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-canopy)")+
#   ylab("Total biomass (MgC)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))

### GROUND BASIS FOR CALC
### plots ordered in bins along MgC/ha-GROUND
contain.npp <- npp.dat2[!is.na(bos.biom30m) & !is.na(lulc), .(sum(hyb.npp.mod.cap, na.rm=T)/2000, 
                                                             sum(fia.empir.npp.kg.hw.forest, na.rm=T)/2000,
                                                             sum(bos.biom30m, na.rm=T)/2000), by=.(biom.bin.ground, lulc)]
contain.npp <- contain.npp[order(biom.bin.ground), ]
contain.npp[, bin.num:=as.numeric(biom.bin.ground)]
names(contain.npp)[3:5] <- c("hyb.npp.MgC", "fia.npp.MgC", "tot.biom.MgC")
contain.npp[,LULC:=as.factor(lulc)]

### with viridis coloring
lulc.pal <- viridis(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])

### HYBRID STACKED NPP
# ggplot(contain.npp, aes(x=bin.num, y=hyb.npp.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-ground)")+
#   ylab("Hybrid, total NPP (MgC/yr)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
#   ylim(0, 1400)

## viridis
ggplot(contain.npp, aes(x=bin.num, y=hyb.npp.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Biomass density (MgC/ha-ground)")+
  ylab("Hybrid, total NPP (MgC/yr)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
  ylim(0, 1400)



### FIA STACKED NPP
# ggplot(contain.npp, aes(x=bin.num, y=fia.npp.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-ground)")+
#   ylab("FIA, total NPP (MgC/yr)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
#   ylim(0, 1400)

## viridis
ggplot(contain.npp, aes(x=bin.num, y=fia.npp.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Biomass density (MgC/ha-ground)")+
  ylab("FIA, total NPP (MgC/yr)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))+
  ylim(0, 1400)


### TOTAL BIOMASS STACKED
# ggplot(contain.npp, aes(x=bin.num, y=tot.biom.MgC, fill=LULC))+
#   geom_area(position='stack')+
#   scale_fill_manual(values = c("forestgreen", "gray55", "salmon", "gold", "chartreuse3", "lightblue"),
#                     name="Land Cover",
#                     breaks=c(1,2,3,4,5,6),
#                     labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
#   xlab("Biomass density (MgC/ha-canopy)")+
#   ylab("Total biomass (MgC)")+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))

## viridis
ggplot(contain.npp, aes(x=bin.num, y=tot.biom.MgC, fill=LULC))+
  geom_area(position='stack')+
  scale_fill_manual(values = lulc.pal,
                    name="Land Cover",
                    breaks=c(1,2,3,4,5,6),
                    labels=c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water"))+
  xlab("Biomass density (MgC/ha-canopy)")+
  ylab("Total biomass (MgC)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25), labels=c(0, 50, 100, 150, 200, 250))



### experimenting with some of the viridis color pallettes
biom <- raster("processed/boston/bos.biom.tif")
lulc <- raster("processed/boston/bos.lulc.lumped.tif")
plot(biom)
image(biom)
library(viridis)
image(biom, col=inferno(256))
image(biom, col=inferno(64))
image(biom, col=viridis(64))
image(biom, col=magma(5))
image(lulc, col=magma(6))
magma(5)
image(lulc, col=inferno(6))
#order: otherveg(5)=6, forest(1)=5, water(6)=1, dev(2)=2, hdres(3)=3, ldres(4)=4
lulc.pal <- inferno(6)
lulc.pal <- viridis(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])
image(lulc, col=lulc.pal) ## very nice

lulc.pal
col2rgb(lulc.pal) ## how to enter these assholes into arc RGB
