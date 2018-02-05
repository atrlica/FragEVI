library(raster)
library(data.table)
library(lattice)
library(rgeos)
library(rgdal)
library(sp)

### initial processing -- evi+isa+lulc+AOI
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
dat.r <- stack("E:/FragEVI/processed/EVI_enhance_stack.tif")
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
pdf(file = "E:/FragEVI/images/EVI_ISA_binned30m.pdf", width = 6, height = 6)
plot(chunk[,i.bin], chunk[,med], 
     ylim=c(noveg-0.05, veg+0.05),
     main="EVI vs. %ISA (binned)", 
     xlab="Urbanization intensity", ylab="EVI",
     pch=15, col="forestgreen")
lines(beta.range, Vzi, col="red", lwd=2, lty=2)
dev.off()




## look at NDVI vs edge class in 1 m data (will need to extract from sub-polys)
### read in various 1m maps
bos.cov <- raster("E:/FragEVI/processed/bos.cov.tif")
bos.ndvi <- crop(raster("E:/FragEVI/data/NDVI/NDVI_1m_res_cangrid.tif", bos.cov)) ## extent is larger than boston cover map
bos.ed1 <- raster("E:/FragEVI/processed/edge10m.tif")
bos.ed2 <- raster("E:/FragEVI/processed/edge20m.tif")
bos.ed3 <- raster("E:/FragEVI/processed/edge30m.tif")

bos.stack <- stack(bos.cov, bos.ndvi, bos.ed1, bos.ed2, bos.ed3)


