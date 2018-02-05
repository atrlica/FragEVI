library(raster)
library(data.table)
library(lattice)
library(rgeos)
library(rgdal)
library(sp)

### initial processing -- evi+isa+lulc+AOI
# dat.r <- stack("processed/evi.isa.30m.tif")
# 
# ### read in the LULC map to help constrain the pixel selection for the end members
# lulc <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/LU_polys/LU_rast_full.tif")
# lulc <- projectRaster(lulc, dat.r, method = "ngb")
# 
# ### AOI boundaries
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# AOI.r <- rasterize(AOI, dat.r)
# dat.r <- stack(dat.r, lulc, AOI.r)
# writeRaster(dat.r, filename="processed/EVI_enhance_stack.tif", format="GTiff", overwrite=T)
dat.r <- stack("E:/FragEVI/processed/EVI_enhance_stack.tif")
names(dat.r) <-  c("evi", "isa", "lulc", "AOI")

## get area of cover, snap to 30m grid in Arc
bos.cov <- raster("E:/FragEVI/processed/bos.cov.tif")

### intermediate prep step to get isolated cover values before snapping to 30 m landsat grid in Arc
# ## get labeled values for grass/canopy/barren, export to arc for aggregation to landsat grid
# grass.find <- function(x){
#   x[x!=1] <- 0
#   return(x) ## this give 900 for areas that are 100% grass
# }
# grass <- calc(bos.cov, fun=grass.find, filename="E:/FragEVI/processed/bos.grass_only.tif", format="GTiff", overwrite=T)
# 
# barr.find <- function(x){
#   x[x==0] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
# }
# barr <- calc(bos.cov, barr.find, filename="E:/FragEVI/processed/bos.barr_only.tif", format="GTiff", overwrite=T)
# 
# can.find <- function(x){
#   x[x==2] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
#   return(x)
# }
# can <- calc(bos.cov, can.find, filename="E:/FragEVI/processed/bos.can_only.tif", format="GTiff", overwrite=T)
# # can.agg <- aggregate(can, fact=30, expand=T, fun=sum, na.rm=T)

### aggregate in the stack from edge class too



### EVI_grid_agg.py to snap these to grid




### this approach controls for which pixels to include as end-members, but not sure that makes sense and creates discontinuities at the ends of the Vzi curve vs. the binned averages
# veg <- dat[isa<0.01 & lulc%in%c(3,37), mean(evi, na.rm=T)]
# veg.err <- dat[isa<0.01 & lulc%in%c(3,37), sd(evi, na.rm=T)/sqrt(length(evi))]
# veg.n <- dat[isa<0.01 & lulc%in%c(3,37), length(evi)]
# noveg <- dat[isa>=0.99 & lulc!=20, mean(evi, na.rm=T)]
# noveg.err <- dat[isa>=0.99 & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
# noveg.n <- dat[isa>=0.99 & lulc!=20, length(evi)]

## alternate approach: just chuck the water values which will fuck everything up (v. low evi, v. low isa)
veg <- dat[isa<0.01 & AOI==1, mean(evi, na.rm=T)]
veg.err <- dat[isa<0.01  & lulc!=20 & AOI==1, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<0.01  & lulc!=20 & AOI==1, length(evi)]
noveg <- dat[isa>=0.99 & lulc!=20 & AOI==1, mean(evi, na.rm=T)]
noveg.err <- dat[isa>=0.99 & lulc!=20 & AOI==1, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=0.99 & lulc!=20 & AOI==1, length(evi)]


# s <- sample(1:dim(dat)[1], size = 10000, replace=F)
# plot(dat[s,isa], dat[s,evi]) ## I see a bend in the middle of this


### set up framework from Zhao
beta.range <- seq(from=0, to=1, by=0.01) ## this has 101 bins .. need to get the bins aligned properly with the evi/isa means
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)
beta.vzi <- data.frame(cbind(beta.range, Vzi))
dat[AOI==1 & lulc!=20, bin:=findInterval(isa, beta.range)]
chunk <- dat[lulc!=20 & AOI==1, list(.N, m=median(evi, na.rm=T), I=median(isa, na.rm=T)), by=bin]

par(mar=c(3,3,2,1))
plot(chunk[,bin], chunk[,m], ylim=c(noveg-0.05, veg+0.05), main="EVI vs. %ISA(bin)", xlab="beta", ylab="EVI")
lines(beta.range*100, Vzi, col="red")
plot(chunk[,bin], chunk[,I], main="ISA vs. %ISA(bin)", ylab="frac. ISA")

## calculate EVI enhancement (either in bins or in raw data?)



## look at NDVI vs edge class in 1 m data (will need to extract from sub-polys)



