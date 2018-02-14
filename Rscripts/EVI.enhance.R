library(raster)
library(data.table)
library(lattice)
library(rgeos)
library(rgdal)
library(sp)

### Stack 30m layers and crop to AOI
# dat.r <- stack("processed/evi.isa.30m.tif") ## just EVI and ISA
# 
# ### read in the LULC map to help constrain the pixel selection for the end members
# lulc <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/LU_polys/LU_rast_full.tif")
# lulc <- projectRaster(lulc, dat.r, method = "ngb")
# 
# ### AOI boundaries
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# AOI.r <- rasterize(AOI, dat.r)
#
# dat.r <- stack(dat.r, lulc, AOI.r)
# writeRaster(dat.r, filename="processed/EVI_enhance_stack.tif", format="GTiff", overwrite=T)
dat.r <- stack("processed/EVI_enhance_stack.tif")
names(dat.r) <-  c("evi", "isa", "lulc", "AOI")
dat <- as.data.table(as.data.frame(dat.r))

### set up framework from Zhao to find Vzi
### find values for end members

### this approach controls for which pixels to include as end-members, but not sure that makes sense and creates discontinuities at the ends of the Vzi curve vs. the binned averages
# veg <- dat[isa<0.01 & lulc%in%c(3,37), mean(evi, na.rm=T)]
# veg.err <- dat[isa<0.01 & lulc%in%c(3,37), sd(evi, na.rm=T)/sqrt(length(evi))]
# veg.n <- dat[isa<0.01 & lulc%in%c(3,37), length(evi)]
# noveg <- dat[isa>=0.99 & lulc!=20, mean(evi, na.rm=T)]
# noveg.err <- dat[isa>=0.99 & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
# noveg.n <- dat[isa>=0.99 & lulc!=20, length(evi)]

## alternate approach: just chuck the water values which will fuck everything up (v. low evi, v. low isa)
veg <- dat[isa<0.01 & AOI==1 & lulc!=20, median(evi, na.rm=T)]
veg.err <- dat[isa<0.01  & lulc!=20 & AOI==1, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<0.01  & lulc!=20 & AOI==1, length(evi)]
noveg <- dat[isa>=0.99 & lulc!=20 & AOI==1, median(evi, na.rm=T)]
noveg.err <- dat[isa>=0.99 & lulc!=20 & AOI==1, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=0.99 & lulc!=20 & AOI==1, length(evi)]

### get range of intensity values and figure a linear Vzi curve
beta.range <- seq(from=0, to=1, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)
dat[AOI==1 & lulc!=20, bin:=findInterval(isa, beta.range, all.inside=T)]
chunk <- dat[lulc!=20 & AOI==1, list(.N, m=mean(evi, na.rm=T), med=median(evi, na.rm=T), i.bin=median(isa, na.rm=T)), by=bin]

par(mar=c(4,4,2,1), mgp=c(2,1,0))
pdf(file = "images/EVI_ISA_binned30m.pdf", width = 6, height = 6)
plot(chunk[,i.bin], chunk[,med], 
     ylim=c(noveg-0.05, veg+0.05),
     main="EVI vs. %ISA (binned)", 
     xlab="Urbanization intensity", ylab="EVI",
     pch=15, col="forestgreen")
lines(beta.range, Vzi, col="red", lwd=2, lty=2)
dev.off()





### analysis of 1m data with edge classes
## look at NDVI vs edge class in 1 m data (will need to extract from sub-polys)
bos.stack <- stack("processed/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")

### load up test polygons for different areas of the city
rox <- readOGR(dsn = "processed/roxbury_test1.shp", layer="roxbury_test1")
sb <- readOGR(dsn = "processed/stonybrook_test1.shp", layer="stonybrook_test1")
allan <- readOGR(dsn = "processed/allandale_test1.shp", layer="allandale_test1")
com <- readOGR(dsn = "processed/commons_test1.shp", layer="commons_test1")

### roxbury
rox.dat <- extract(bos.stack, rox, df=T)
rox.dat <- as.data.table(rox.dat)

rox.dat[cov==2 & ed10==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(rox.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10m
rox.dat[cov==2 & is.na(ed10) & ed20==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(rox.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10-20m
rox.dat[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(rox.dat[cov==2,])[1])] ## ndvi and relative canopy % of 20-30m
rox.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), .(median(ndvi, na.rm=T), length(ndvi)/dim(rox.dat[cov==2,])[1])] ## ndvi and relative canopy % of >30m
rox.dat[, .(median(ndvi, na.rm=T), length(ndvi)/dim(rox.dat)[1]), by=cov] ### ndvi by cover and relative total fraction
rox.dat[cov==0 & isa==1, .(median(ndvi), length(ndvi)/dim(rox.dat[cov==0])[1])] ## ndvi of isa, relative fraction of barren
rox.dat[cov==0 & isa==0, .(median(ndvi), length(ndvi)/dim(rox.dat[cov==0])[1])] ## ndvi of non-impervious barren, relative fraction of barren

### allandale reserve
plot(allan, add=T)
allan.dat <- extract(bos.stack, allan, df=T)
allan.dat <- as.data.table(allan.dat)

allan.dat[cov==2 & ed10==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(allan.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10m
allan.dat[cov==2 & is.na(ed10) & ed20==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(allan.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10-20m
allan.dat[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(allan.dat[cov==2,])[1])] ## ndvi and relative canopy % of 20-30m
allan.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), .(median(ndvi, na.rm=T), length(ndvi)/dim(allan.dat[cov==2,])[1])] ## ndvi and relative canopy % of >30m
allan.dat[, .(median(ndvi, na.rm=T), length(ndvi)/dim(allan.dat)[1]), by=cov] ### ndvi by cover and relative total fraction
allan.dat[cov==0 & isa==1, .(median(ndvi), length(ndvi)/dim(allan.dat[cov==0])[1])] ## ndvi of isa, relative fraction of barren
allan.dat[cov==0 & isa==0, .(median(ndvi), length(ndvi)/dim(allan.dat[cov==0])[1])] ## ndvi of non-impervious barren, relative fraction of barren

### stonybrook
plot(sb, add=T)
sb.dat <- extract(bos.stack, sb, df=T)
sb.dat <- as.data.table(sb.dat)

sb.dat[cov==2 & ed10==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(sb.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10m
sb.dat[cov==2 & is.na(ed10) & ed20==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(sb.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10-20m
sb.dat[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(sb.dat[cov==2,])[1])] ## ndvi and relative canopy % of 20-30m
sb.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), .(median(ndvi, na.rm=T), length(ndvi)/dim(sb.dat[cov==2,])[1])] ## ndvi and relative canopy % of >30m
sb.dat[, .(median(ndvi, na.rm=T), length(ndvi)/dim(sb.dat)[1]), by=cov] ### ndvi by cover and relative total fraction
sb.dat[cov==0 & isa==1, .(median(ndvi), length(ndvi)/dim(sb.dat[cov==0])[1])] ## ndvi of isa, relative fraction of barren
sb.dat[cov==0 & isa==0, .(median(ndvi), length(ndvi)/dim(sb.dat[cov==0])[1])] ## ndvi of non-impervious barren, relative fraction of barren

### commons
plot(com, add=T)
com.dat <- extract(bos.stack, com, df=T)
com.dat <- as.data.table(com.dat)

com.dat[cov==2 & ed10==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(com.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10m
com.dat[cov==2 & is.na(ed10) & ed20==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(com.dat[cov==2,])[1])] ## ndvi and relative canopy % of 10-20m
com.dat[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(com.dat[cov==2,])[1])] ## ndvi and relative canopy % of 20-30m
com.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), .(median(ndvi, na.rm=T), length(ndvi)/dim(com.dat[cov==2,])[1])] ## ndvi and relative canopy % of >30m
com.dat[, .(median(ndvi, na.rm=T), length(ndvi)/dim(com.dat)[1]), by=cov] ### ndvi by cover and relative total fraction
com.dat[cov==0 & isa==1, .(median(ndvi), length(ndvi)/dim(com.dat[cov==0])[1])] ## ndvi of isa, relative fraction of barren
com.dat[cov==0 & isa==0, .(median(ndvi), length(ndvi)/dim(com.dat[cov==0])[1])] ## ndvi of non-impervious barren, relative fraction of barren







### extract test values from AOI over Stonybrook reserve
sb <- readOGR(dsn = "processed/stonybrook_AOI.shp", layer = "stonybrook_AOI")
sb.dat <- extract(bos.stack, sb, df=T)
sb.dat <- as.data.table(sb.dat)

sb.dat[cov==2 & ed10==1, median(ndvi)]
hist(sb.dat[cov==2 & ed10==1, ndvi])

sb.dat[cov==2 & ed20==1, median(ndvi)]
hist(sb.dat[cov==2 & ed20==1, ndvi])

sb.dat[cov==2 & ed30==1, median(ndvi)]
hist(sb.dat[cov==2 & ed30==1, ndvi])

### only interior canopy
sb.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), median(ndvi)] ### higher than 
hist(sb.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), ndvi])

sb.dat[ed30==1 & is.na(ed20) & is.na(ed10), ring:=3]
sb.dat[ed20==1 & is.na(ed10), ring:=2]
sb.dat[ed10==1, ring:=1]
sb.dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), ring:=4]
sb.dat[cov==1, ring:=5]
sb.dat[cov==0, ring:=6]

par(mar=c(1,4,2,1), mgp=c(2,1,0))
png(filename="images/SB_test_NDVI_boxpl.png", width = 6, height = 6, units = "in", res=108)
boxplot(sb.dat[, ndvi]~sb.dat[,ring], 
        col=c("pink", "orange", "yellow", "forestgreen", "lightgreen", "grey50"), xaxt="n",
        main="NDVI by edge distance and cover class", ylab="NDVI")
legend(x = 3, y = 0, legend = c("10m", "20m", "30m", ">30m", "grass", "barren"), 
       fill = c("pink", "orange", "yellow", "forestgreen", "lightgreen", "grey50"), cex=0.84)
dev.off()
