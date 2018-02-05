### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)


#########
### Boston high-res data
#### chunk: combine 2.4m NDVI and canopy map to produce 1m veg classification map
### use Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map

### step 1: 1m Canopy presence/absence map 
bos.can <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")
### set no canopy (0) values to NA
# bos.can.dat <- as.data.table(as.data.frame(bos.can))
# bos.can.dat[bostoncanopy_1m==0,] <- NA
# bos.can.na <- raster(bos.can)
# bos.can.na <- setValues(bos.can.na, bos.can.dat$bostoncanopy_1m)
# writeRaster(bos.can.na, "E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif", format="GTiff", overwrite=T, datatype="INT1U")
# bos.can.na <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif")

### step 2: read in 2.4m Quickbird NDVI map, filter for NA
# bos.ndvi <- raster("E:/FragEVI/data/NDVI/NDVI.img") ## in UTM19N, original Quickbird NDVI
# bos.ndvi.dat <- as.data.table(as.data.frame(bos.ndvi))
# bos.ndvi.dat[NDVI==0, NDVI:=NA]
# bos.ndvi.na <- setValues(bos.ndvi, values = bos.ndvi.dat[,NDVI], filename="E:/FragEVI/data/NDVI/NDVI_na.tif", format="GTiff", overwrite=T) ## good, full map, UTM 19N

### step 3: put NDVI app through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
pyth.path = './Rscripts/NDVI_resamp.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
print(paste("ArcPy working on NDVI resample: ", output))

### step 4: land cover classification for Boston using the 1m resampled NDVI (0 = barren, 1 = grass, 2 = canopy)
bos.ndvi <- raster("E:/FragEVI/data/NDVI/NDVI_1m_res_cangrid.tif")
bos.ndvi <- crop(bos.ndvi, bos.can.na)
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
bos.cov <- cover.bl(bos.can, bos.ndvi, filename="E:/FragEVI/processed/bos.cov.tif")

#####
### Processing for 1m canopy map for edge class

### call python script for identifying canopy edge distance
pyth.path = './Rscripts/canopy_process.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
print(paste("ArcPy working on canopy edges: ", output))

### read in edge rasters and correct to produce proper edge layers
bos.cov <- raster("E:/FragEVI/processed/bos.cov.tif")
ed1 <- raster("E:/FragEVI/processed/nocan_10mbuff.tif")
ed2 <- raster("E:/FragEVI/processed/nocan_20mbuff.tif")
ed3 <- raster("E:/FragEVI/processed/nocan_30mbuff.tif")
ed1 <- extend(ed1, bos.cov) # edge and cover on same grid already
ed2 <- extend(ed2, bos.cov)
ed3 <- extend(ed3, bos.cov)
edges <- stack(ed1, ed2, ed3, bos.cov) # just to make damn sure everything lines up

### Pull in 1m edge buffers and use cover map to correct
edges.bl <- function(x, y, filename) { # x is edge class, y is cover class
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
    g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## cover map
    v[v!=0] <- NA # kill any weird values that aren't coming from the nocan==0 buffer
    v[g!=2] <- NA # cancel edge ID for non-canopy
    v[v==0] <- 1 # ID edge pixels as 1
    v[v!=1 & g%in%c(0,1,2)] <- 0 ## set non-edge values to 0 to hold space
    out <- writeValues(out, v, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
s <- edges.bl(ed1, bos.cov, filename="E:/FragEVI/processed/edge10m.tif")
t <- edges.bl(ed2, bos.cov, filename="E:/FragEVI/processed/edge20m.tif")
u <- edges.bl(ed3, bos.cov, filename="E:/FragEVI/processed/edge30m.tif")


###### get aggregated areas for 1m data at 30 m grid
#######
### the following chunk aggregates 1m ISA to 30m landsat grid
### intermediate steps run on the cluster
### read in AOI for crs description
AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
master.crs <- crs(AOI)

### load EVI composite for grid and process for ISA grid
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif") ## this is the July 2010-2012 AOI EVI composite
# evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
# writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)

#### The following should eventually be redone using the arcpy snap-to functionality -- could be artifacts with how this was made
### process ISA fraction per 30 m EVI gridcell (~Urban intensity?)
### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- too hard to monkey with, leave in its native state until the end
isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")

### raster-native approach -- takes a long time
### correct the 16 values to NA to get a proper aggregated mean value per cell (0 = not impervious, 1 = impervious)
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.30m.tif", format="GTiff", overwrite=T)

### read in 30 m ISA, align to the EVI grid and stack
isa.na.agg <- raster("processed/isa.30m.tif")
evi.r <- raster("E:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif") ## this is the July 2010-2012 AOI EVI composite
isa.pr <- projectRaster(isa.na.agg, evi.r, method="ngb", filename="processed/isa.30m.NAD83.tif", format="GTiff", overwrite=T)

frag.stack <- stack(evi.r, isa.pr)
frag.stack <- crop(frag.stack, AOI)
writeRaster(frag.stack, filename="processed/evi.isa.30m.tif", format="GTiff", overwrite=T)

bos.cov <- raster("E:/FragEVI/processed/bos.cov.tif")

#### pull out cover classes as individual layers to assist with aggregated area in 30 m pixels
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



### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207009_ne_19_h_20160706/m_4207009_ne_19_h_20160706.tif")
naip.dat <- as.data.table(as.data.frame(naip.test))
names(naip.dat) <- c("band1", "band2", "band3", "band4")
naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
ndvi.r <- raster(naip.test[[1]])
ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)
