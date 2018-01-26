### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
# library(ff)

# setwd("/projectnb/buultra/atrlica/FragEVI/")
### load EVI composite for grid and process for ISA grid
evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif")
# evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
# writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)


## extent object is c(xmin, xmax, ymin, ymax)

### get ISA fraction per grid (~Urban intensity?)
### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- 
### too hard to monkey with, leave in its native state until the end

# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# # isa.dat <- ff(vmode="double", filename="processed/isa.raw.ff")
# # isa.dat <- as.vector(isa)
# # isa.dat[isa.dat==16] <- NA
# # isa <- setValues(isa, isa.dat)
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# 
# ### aggregate to 30m 
# isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.30m.tif", format="GTiff", overwrite=T)
isa.na.agg <- raster("processed/isa.30m.tif")
evi.isa <- projectRaster(from=evi.r, crs=crs(isa.na.agg), res=30, method="bilinear")
evi.isa.cr <- crop(evi.isa, extent(isa.na.agg))
writeRaster(evi.isa.cr, filename="processed/evi.isa.tif", format="GTiff", overwrite=T)

frag.stack <- stack(evi.isa.cr, isa.na.agg)

### read in AOI for crs description
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# master.crs <- crs(AOI)


### Boston test areas -- combine high-res NDVI and canopy map

#### important notes re. land cover class using Quickbird 2.4m NDVI w/ 1m canopy presence map
### not readily possible in R to align/resample the two different maps into a common grid -- severe registration error in NDVI after processing
### used Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map

## old attempted approach in R -- results in bad alignment and crummy classification maps
# ### read in 2.4m Quickbird NDVI map, filter for NA -- NOTE: working on the NDVI in R messes it up too badly to use
# bos.ndvi <- raster("E:/FragEVI/data/NDVI/NDVI.img") ## original Quickbird NDVI, manually process for 0==NA
# bos.ndvi.dat <- as.data.table(as.data.frame(bos.ndvi))
# bos.ndvi.dat[NDVI==0, NDVI:=NA]
# bos.ndvi.na <- setValues(bos.ndvi, values = bos.ndvi.dat[,NDVI]) ## good, full map, UTM 19N
# writeRaster(bos.ndvi.na, filename="E:/FragEVI/data/NDVI/NDVI_na.tif", format="GTiff", overwrite=T)

### new NDVI resampled in arc
ndvi.res <- raster("E:/FragEVI/data/NDVI/NDVI_1m_res_cangrid.tif")

## 1m Canopy presence/absence map
bos.can <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")

### test polygons imported from Arc (approx 100x100m)
t1 <- readOGR(dsn="E:/FragEVI/data/AOI", layer = "test.sm1")
t1 <- spTransform(t1, crs(bos.can))
t2 <- readOGR(dsn="E:/FragEVI/data/AOI", layer = "test.sm2")
t2 <- spTransform(t2, crs(bos.stack))
t3 <- readOGR(dsn="E:/FragEVI/data/AOI", layer = "test.sm3")
t3 <- spTransform(t3, crs(bos.stack))
t4 <- readOGR(dsn="E:/FragEVI/data/AOI", layer = "test.sm4")
t4 <- spTransform(t4, crs(bos.stack))

par(mfrow=c(1,2))

### land cover classification using the good NDVI (0 = barren, 1 = grass, 2 = canopy)
ndvi.cr <- crop(ndvi.res, bos.can)
bos <- stack(bos.can, ndvi.cr)
# bos.dat <- as.data.table(as.data.frame(bos)) ## Poland cannot into space, too large to do all at once

### test 1, good NDVI 1m, SW corner of Public Gardens -- seems to find grassy/barren areas between trees alright, sees impervious as barren
can.t1 <- crop(bos.can, t1)
ndvi.t1 <- crop(ndvi.cr, t1)
plot(can.t1); plot(t1, add=T)
plot(ndvi.t1); plot(t1, add=T)
t1.dat <- as.data.table(as.data.frame(stack(can.t1, ndvi.t1)))
names(t1.dat) <- c("can", "ndvi")
t1.dat[ndvi<0.2, cov:=0]
t1.dat[can==0 & ndvi>0.2, cov:=1]
t1.dat[can==1, cov:=2]
cov.t1 <- setValues(can.t1, t1.dat[,cov])
plot(can.t1); plot(t1, add=T)
plot(cov.t1); plot(t1, add=T)

### test 2, good NDVI 1m, Alston Suburb --some marginal confusion near canopies but also sees what look like residential lawns
can.t2 <- crop(bos.can, t2)
ndvi.t2 <- crop(ndvi.cr, t2)
plot(can.t2); plot(t2, add=T)
plot(ndvi.t2); plot(t2, add=T)
t2.dat <- as.data.table(as.data.frame(stack(can.t2, ndvi.t2)))
names(t2.dat) <- c("can", "ndvi")
t2.dat[ndvi<0.2, cov:=0]
t2.dat[can==0 & ndvi>0.2, cov:=1]
t2.dat[can==1, cov:=2]
cov.t2 <- setValues(can.t2, t2.dat[,cov])
plot(can.t2); plot(t2, add=T)
plot(cov.t2); plot(t2, add=T)

### test 3, good NDVI 1m, forest edge with ball park, Franklin Park S of track -- detects clear edge with lawn/baseball field, even detects what look like barren places at park periphery
can.t3 <- crop(bos.can, t3)
ndvi.t3 <- crop(ndvi.cr, t3)
plot(can.t3); plot(t3, add=T)
plot(ndvi.t3); plot(t3, add=T)
t3.dat <- as.data.table(as.data.frame(stack(can.t3, ndvi.t3)))
names(t3.dat) <- c("can", "ndvi")
t3.dat[ndvi<0.2, cov:=0]
t3.dat[can==0 & ndvi>0.2, cov:=1]
t3.dat[can==1, cov:=2]
cov.t3 <- setValues(can.t3, t3.dat[,cov])
plot(can.t3); plot(t3, add=T)
plot(cov.t3); plot(t3, add=T)

### test 4, good NDVI 1m, forest patch in Stonybrook Park S of Turtle Pond -- sees canopy except for a few gaps with grass
ndvi.t4 <- crop(ndvi.cr, t4)
plot(can.t4); plot(t4, add=T)
plot(ndvi.t4); plot(t4, add=T)
t4.dat <- as.data.table(as.data.frame(stack(can.t4, ndvi.t4)))
names(t4.dat) <- c("can", "ndvi")
t4.dat[ndvi<0.2, cov:=0]
t4.dat[can==0 & ndvi>=0.2, cov:=1]
t4.dat[can==1, cov:=2]
cov.t4 <- setValues(can.t4, t4.dat[,cov])
plot(can.t4); plot(t4, add=T)
plot(cov.t4); plot(t4, add=T) ## color scale goofballed -- note no barren present


### attempt to perform Boston-wide cover class, row by row
strip <- seq(from=1, to=nrow(bos.can), by = 2000)
cov <- numeric()
for(r in 1:nrow(bos.can)){
  tmp.can <- getValues(bos.can, row=r)
  tmp.ndvi <- getValues(ndvi.cr, row=r)
  d <- as.data.table(cbind(tmp.can, tmp.ndvi))
  d[tmp.ndvi<0.2, cov:=0]
  d[tmp.can==0&tmp.ndvi>=0.2, cov:=1]
  d[tmp.can==1, cov:=2]
  cov <- c(cov, d[,cov])
  print(paste("finished row", r)) ## not particularly fast!
}
bos.cov <- setValues(bos.can, cov)
writeRaster(bos.cov, filename="E:/FragEVI/processed/Boston_1mcover.tif", format="GTiff", overwrite=T)

### function for row-by-row replacement of raster values (this works but the logic is flawed for the cover repalcement)
cover.f <- function(x, y, filename) { #x is canopy, y is ndvi
  out <- raster(x)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (r in 1:nrow(out)) {
    v <- getValues(x, r)
    g <- getValues(y, r)
    
    v[g<0.2] <- 0
    v[v==0 & g>=0.2] <- 1
    v[v==1] <- 2
    out <- writeValues(out, v, r)
    print(paste("finished row", r))
  }
  out <- writeStop(out)
  return(out)
}
s <- cover.f(bos.can, ndvi.cr, filename="E:/FragEVI/processed/bos.cov.test.tif")
par(mfrow=c(1,1))
plot


### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207009_ne_19_h_20160706/m_4207009_ne_19_h_20160706.tif")
naip.dat <- as.data.table(as.data.frame(naip.test))
names(naip.dat) <- c("band1", "band2", "band3", "band4")
naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
ndvi.r <- raster(naip.test[[1]])
ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)
