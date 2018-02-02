### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(rPython)
# library(ff)


#######
### the following chunk aggregates 1m ISA to 30m landsat grid
### intermediate steps run on the cluster

### read in AOI for crs description
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# master.crs <- crs(AOI)

### load EVI composite for grid and process for ISA grid
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif") ## this is the July 2010-2012 AOI EVI composite
# evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
# writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)

### process ISA fraction per 30 m EVI gridcell (~Urban intensity?)
### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- too hard to monkey with, leave in its native state until the end
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")

# ### the data-only approach doesn't seem to work -- file is too large to read in, possibly see big.memory
# # isa.dat <- ff(vmode="double", filename="processed/isa.raw.ff")
# # isa.dat <- as.vector(isa)
# # isa.dat[isa.dat==16] <- NA
# # isa <- setValues(isa, isa.dat)

### raster-native approach -- takes a long time
### correct the 16 values to NA to get a proper aggregated mean value per cell (0 = not impervious, 1 = impervious)
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.30m.tif", format="GTiff", overwrite=T)

### read in the 30 m ISA, align to the EVI grid and stack
isa.na.agg <- raster("E:/FragEVI/processed/stack/isa.30m.tif")
# evi.isa <- projectRaster(from=evi.r, crs=crs(isa.na.agg), res=30, method="bilinear")
# evi.isa.cr <- crop(evi.isa, extent(isa.na.agg))
# writeRaster(evi.isa.cr, filename="processed/evi.isa.tif", format="GTiff", overwrite=T)
evi.r <- raster("E:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif") ## this is the July 2010-2012 AOI EVI composite
isa.stack <- projectRaster(isa.na.agg, evi.r, method="bilinear")

frag.stack <- stack(evi.r, isa.stack)
names(frag.stack) <- c("evi", "isa")
writeRaster(frag.stack, filename="E:/FragEVI/processed/EVI.ISA.30mstack.tiff", format="GTiff", overwrite=T)


### Boston -- combine high-res NDVI and canopy map to produce 1m veg classification map

#### important notes re. land cover class using Quickbird 2.4m NDVI w/ 1m canopy presence map
### not readily possible in R to align/resample the two different maps into a common grid -- severe registration error in NDVI after processing
### use Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map

### step 1: 1m Canopy presence/absence map -- set no canopy (0) values to NA
# bos.can <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")
# bos.can.dat <- as.data.table(as.data.frame(bos.can))
# bos.can.dat[bostoncanopy_1m==0,] <- NA
# bos.can.na <- raster(bos.can)
# bos.can.na <- setValues(bos.can.na, bos.can.dat$bostoncanopy_1m)
# writeRaster(bos.can.na, "E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif", format="GTiff", overwrite=T, datatype="INT1U")
bos.can.na <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif")


### step 2: read in 2.4m Quickbird NDVI map, filter for NA
bos.ndvi <- raster("E:/FragEVI/data/NDVI/NDVI.img") ## original Quickbird NDVI, manually process for 0==NA
bos.ndvi.dat <- as.data.table(as.data.frame(bos.ndvi))
bos.ndvi.dat[NDVI==0, NDVI:=NA]
bos.ndvi.na <- setValues(bos.ndvi, values = bos.ndvi.dat[,NDVI]) ## good, full map, UTM 19N
writeRaster(bos.ndvi.na, filename="E:/FragEVI/data/NDVI/NDVI_na.tif", format="GTiff", overwrite=T)

### step 3: put NDVI app through arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
### NEED PYTHON IMBED HERE

### step 4: land cover classification for Boston using the 1m resampled NDVI (0 = barren, 1 = grass, 2 = canopy)
bos.ndvi <- raster("E:/FragEVI/data/NDVI/NDVI_1m_res_cangrid.tif")
bos.ndvi <- crop(bos.ndvi, bos.can.na)

### function to process serially in chunks #### THIS WORKS BEAUTIFULLY, TAKES ~15 MIN FOR ALL OF BOSTON
cover.bl <- function(x, y, filename) { # x is canopy, y is ndvi
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## canopy map
    g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi map
    cov <- v
    cov[g>=0.2 & v==0] <- 1
    cov[v==1] <- 2
    out <- writeValues(out, cov, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
s <- cover.bl(bos.can.na, bos.ndvi, filename="processed/bos.cov.tif")
plot(s)

# ### function for row-by-row replacement of raster values (this works but is slower)
# cover.f <- function(x, y, filename) { #x is canopy, y is ndvi
#   out <- raster(x)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (r in 1:nrow(out)) {
#     v <- getValues(x, r) ## canopy map
#     g <- getValues(y, r) ## ndvi map
#     cov <- v
#     cov[g>=0.2 & v==0] <- 1
#     cov[v==1] <- 2
#     out <- writeValues(out, cov, r)
#     print(paste("finished row", r))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# s <- cover.f(bos.can, ndvi.cr, filename="E:/FragEVI/processed/bos.cov.test.tif")
# plot(s)


#### use canopy map to get buffer zones (10 20 30 m)
bos.can <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")
 
# call python script that handles scrubbing the canopy gaps for tiny gaps, then exports the 10 20 30 m buffers
pyth.path = './Rscripts/canopy_process.py'
output = system2('python.exe', args=pyth.path, stdout=TRUE)
print(paste("ArcPy working on canopy gaps: ", output))




### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207009_ne_19_h_20160706/m_4207009_ne_19_h_20160706.tif")
naip.dat <- as.data.table(as.data.frame(naip.test))
names(naip.dat) <- c("band1", "band2", "band3", "band4")
naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
ndvi.r <- raster(naip.test[[1]])
ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)
