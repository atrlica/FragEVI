### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)


#########
### Boston high-res data
### combine 2.4m NDVI and canopy map to produce 1m veg classification map
### use Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map

# ### step 1: 1m Canopy presence/absence map 
# bos.can <- raster("data/dataverse_files/bostoncanopy_1m.tif")

# ### set no canopy (0) values to NA
# bos.can.dat <- as.data.table(as.data.frame(bos.can))
# bos.can.dat[bostoncanopy_1m==0,] <- NA
# bos.can.na <- raster(bos.can)
# bos.can.na <- setValues(bos.can.na, bos.can.dat$bostoncanopy_1m)
# writeRaster(bos.can.na, "E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif", format="GTiff", overwrite=T, datatype="INT1U")
# bos.can.na <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif")

# ### step 3 (must be performed on desktop): put NDVI.img through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
# pyth.path = './Rscripts/NDVI_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(paste("ArcPy working on NDVI resample: ", output))

# ### step 4: land cover classification for Boston using the 1m resampled NDVI (0 = barren, 1 = grass, 2 = canopy)
# bos.ndvi <- raster("data/NDVI/NDVI_1m_res_cangrid.tif")
# bos.can.na <- raster("data/dataverse_files/bostoncanopy_1m_na.tif")
# bos.ndvi <- crop(bos.ndvi, bos.can.na)
# cover.bl <- function(x, y, filename) { # x is canopy, y is ndvi
#   if(file.exists("processed/bos.cov.tif")){file.remove("processed/bos.cov.tif")}
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## canopy map
#     g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi map
#     cov <- v
#     cov[g>=0.2 & v==0] <- 1
#     cov[v==1] <- 2
#     out <- writeValues(out, cov, bs$row[i])
#     print(paste("finished 1m cover block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# bos.cov <- cover.bl(bos.can, bos.ndvi, filename="processed/bos.cov.tif")

# ### Processing 1m canopy map for edge class
# ### call python script for identifying canopy edge distance (desktop only)
# pyth.path = './Rscripts/canopy_process.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)
# 
# ### read in edge rasters and correct classifications to produce raster of edge rings (>10, 10-20, 20-30, >30)
# bos.cov <- raster("processed/bos.cov.tif")
# ed1 <- raster("processed/nocan_10mbuff.tif")
# ed2 <- raster("processed/nocan_20mbuff.tif")
# ed3 <- raster("processed/nocan_30mbuff.tif")
# ed1 <- extend(ed1, bos.cov) # edge and cover on same grid already
# ed2 <- extend(ed2, bos.cov)
# ed3 <- extend(ed3, bos.cov)
# # edges <- stack(ed1, ed2, ed3, bos.cov) # just to make damn sure everything lines up
# 
# edges.bl <- function(x, y, filename) { # x is edge class, y is cover class
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
#     g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## cover map
#     v[v!=0] <- NA # kill any weird values that aren't coming from the nocan==0 buffer map
#     v[g!=2] <- NA # cancel edge ID for non-canopy
#     v[v==0] <- 1 # ID edge pixels as 1
#     v[v!=1] <- 0 ## set non-edge values to 0 to hold space
#     out <- writeValues(out, v, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# s <- edges.bl(ed1, bos.cov, filename="processed/edge10m.tif")
# t <- edges.bl(ed2, bos.cov, filename="processed/edge20m.tif")
# u <- edges.bl(ed3, bos.cov, filename="processed/edge30m.tif")

### package up all the 1m data for Boston and write to disk
ed10 <- raster("processed/edge10m.tif")
ed20 <- raster("processed/edge20m.tif")
ed30 <- raster("processed/edge30m.tif")
bos.cov <- raster("processed/bos.cov.tif")
bos.can <- raster("data/dataverse_files/bostoncanopy_1m.tif")
bos.ndvi <- raster("data/NDVI/NDVI_1m_res_cangrid.tif")
bos.ndvi <- crop(bos.ndvi, bos.can)

# ### add in 1m ISA (go get it, crop it, then reproject)
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# towns <- readOGR(dsn = "/projectnb/buultra/atrlica/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# bos.bord <- towns[towns@data$TOWN=="BOSTON",]
# bos.bord <- spTransform(bos.bord, crs(isa))
# bos.isa <- crop(isa, bos.bord)
# writeRaster(bos.isa, filename="processed/bos_isa_1m.tif", format="GTiff", overwrite=T)
# # python.load("Rscripts/bosISA_resamp.py") ### test this
bos.isa <- raster("processed/isa_cangrid.tif")
bos.isa <- crop(bos.isa, bos.can)

bos.stack <- stack(bos.can, bos.ndvi, bos.cov, bos.isa, ed10, ed20, ed30)
# writeRaster(bos.stack, "processed/boston/bos.stack.1m.tif", format="GTiff", overwrite=T)

##### work to get 1m boston data aggregated to 30m EVI grid
evi.stack <- stack("processed/EVI30m.stack.tif")
bos.stack <- stack("processed/boston/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
# plot(bos.stack)

### get Boston limits, raterize to EVI 30m grid
towns <- readOGR(dsn = "F:/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
bos.AOI <- towns[towns@data$TOWN=="BOSTON",]
bos.AOI <- bos.AOI[bos.AOI@data$SHAPE_AREA>1E07,] ## remove Harbor Islands
bos.AOI <- spTransform(bos.AOI, crs(bos.stack))
bos.AOI@data$include <- 1
# bos.AOI.r <- rasterize(bos.AOI[bos.AOI@data$include], evi.stack[[1]])
# bos.AOI.r <- crop(bos.AOI.r, bos.AOI)
# ba.dat <- as.data.table(as.data.frame(bos.AOI.r))
# ba.dat[!is.na(OBJECTID) | !is.na(OBJECTID.1), dogshit:=1] ## set all areas to same value
# bos.AOI.r <- setValues(bos.AOI.r, ba.dat$dogshit)
# 
# bos.evi <- crop(evi.stack, bos.AOI.r)
# bos.evi <- mask(bos.evi, bos.AOI.r)
# names(bos.evi) <- c("evi", "isa", "lulc", "AOI") ### the Boston city limits in 30m grid
# sum(getValues(bos.evi[["AOI"]]), na.rm=T) ## 138k pixels

##### Work boston 1m layers out into separate 0/1 classified files
### make 1m Boston AOI raster mask
bos.AOI.1m <- rasterize(bos.AOI[bos.AOI@data$include], bos.stack[[1]])
bos.AOI.1m <- crop(bos.AOI.1m, bos.AOI)
bos.AOI.1mf <- raster(bos.AOI.1m)
bos.AOI.1mf <- setValues(bos.AOI.1mf, values = rep(1, length=ncell(bos.AOI.1m)))
bos.AOI.1mf <- mask(bos.AOI.1mf, bos.AOI.1m)
writeRaster(bos.AOI.1mf, filename="processed/boston/bos.aoi_only.tif", format="GTiff", overwrite=T)
bos.aoi <- raster("processed/boston/")

# ### get separate layers for each with cover=0/1 (vs. cover=0,1,2)
# bos.stack <- stack("processed/boston/bos.stack.1m.tif")
# names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
# 
# grass.find <- function(x){
#   x[x!=1] <- 0
#   return(x)
# }
# grass.only <- calc(bos.stack[["cov"]], fun=grass.find, filename="processed/boston/bos.grass_only.tif", format="GTiff", overwrite=T)
# 
# barr.find <- function(x){
#   x[x==0] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
# }
# barr.only <- calc(bos.stack[["cov"]], barr.find, filename="processed/boston/bos.barr_only.tif", format="GTiff", overwrite=T)

# can.find <- function(x){
#   x[x==2] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
#   return(x)
# }
# can <- calc(bos.stack[["cov"]], can.find, filename="processed/bos.can_only.tif", format="GTiff", overwrite=T)

# can.only <- raster("processed/bos.can_only.tif")
grass.only <- raster("processed/boston/bos.grass_only.tif")
barr.only <- raster("processed/boston/bos.barr_only.tif")

# nonimp.only <- overlay(barr.only, bos.stack[["isa"]], fun=function(x,y){return(x-y)}, filename="processed/boston/bos.nonimp_only.tif", overwrite=T, format="GTiff")
nonimp.only <- raster("processed/boston/bos.nonimp_only.tif")

# # test if can_only == can
# can.test <- overlay(bos.stack[["can"]], can.only, fun=function(x,y){return(x-y)})
# plot(can.test) # yes they are identical

### make map of individual buffer rings
bos.stack <- stack("processed/boston/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
bos.aoi <- raster("processed/boston/bos.aoi_only.tif")

### first need to reclass the edge maps to 0/1 (are in 1/NA)
fix.buff <- function(x){ # x is buffer, y is next smallest buffer
  x[is.na(x)] <- 0
  return(x)
}
bos.stack[["ed10"]] <- calc(bos.stack[["ed10"]], fun=fix.buff)
bos.stack[["ed20"]] <- calc(bos.stack[["ed20"]], fun=fix.buff)
bos.stack[["ed30"]] <- calc(bos.stack[["ed30"]], fun=fix.buff)
bos.stack <- mask(bos.stack, bos.aoi)

### 10m ring is already done in ed10

### 20m ring
buff.20only <- overlay(bos.stack, buff.20only, fun=function(x,y){return(y-x)})

### 30m ring
buff.30only <- calc(bos.stack[["ed30"]], fix.buff)
buff.30only <- overlay(buff.20only, buff.30only, fun=function(x,y){return(y-x)})

### interior (buff.30only is already 0/1 classed)
buff.Intonly <- overlay(buff.30only, bos.stack[["can"]], fun=function(x,y){return(y-x)})

buffs.only <- stack(buff.10only, buff.20only, buff.30only, buff.Intonly)
buffs.only <- mask(buffs.only, bos.stack[["aoi"]])
writeRaster(buffs.only, filename="processed/bos.buffs_only.tif", format="GTiff", overwrite=T)

bos.stack <- stack(bos.stack, buffs.only)
writeRaster(bos.stack, filename="procesed/bos.stack.tif", format="GTiff", overwrite=T)
### boston stack is now: names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30", "aoi", "grass", "barr", "nonimp")


### fix inclusive buffer rings 0-20m and 0-30m to 0/1
bos.stack <- stack("processed/bos.stack.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30", "aoi", "grass", "barr", "nonimp", "ed10only", "ed20only", "ed30only", "Intonly")

ed20.f <- fix.buff(bos.stack[["ed20"]])
ed20.f <- mask(ed20.f, )


bos.stack <- crop(bos.stack, bos.AOI.1mf)
bos.stack <- stack(bos.stack, bos.AOI.1mf)
bos.stack <- mask(bos.stack, bos.AOI.1m) ## masks the NA/1 AOI by the 1/2 AOI
writeRaster(bos.stack, filename="processed/bos.stack.1m.tif", overwrite=T, format="GTiff")

### final 1m stack for fractional cover (0/1 classified) should be can, ndvi, isa, nonimp, barr, aoi, <10edge, <20edge, <30edge, >30edge, 
## get fractional area of each cover class per 30m aggregate cell
bos.AOI.sizes <- aggregate(bos.stack[["aoi"]], fact=30, fun=sum, na.rm=T) ## gives the size of each cell within the boston AOI
bos.stack.30m <- aggregate(bos.stack, fun=mean, na.rm=T) ## gives the mean (ndvi) or or total area (can/isa/edge/cover) per cell within the boundaries


stack.frac <- function(x, y, filename) { # x is 1m boston input, y is cell area
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    r <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## input class area
    a <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## total inside-AOI area
    v[v!=0] <- NA # kill any weird values that aren't coming from the nocan==0 buffer map
    v[g!=2] <- NA # cancel edge ID for non-canopy
    v[v==0] <- 1 # ID edge pixels as 1
    v[v!=1] <- 0 ## set non-edge values to 0 to hold space
    frac <- r/a
    out <- writeValues(out, frac, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
# s <- edges.bl(ed1, bos.cov, filename="processed/edge10m.tif")
# t <- edges.bl(ed2, bos.cov, filename="processed/edge20m.tif")
# u <- edges.bl(ed3, bos.cov, filename="processed/edge30m.tif")


######
####
### AOI aggregate to 30 m Landsat grid -- only have ISA AOI LULC and EVI as of Feb 22 2018
# ### read in AOI for crs description
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# master.crs <- crs(AOI)

### load EVI composite for grid and process for ISA grid
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif") ## this is the July 2010-2012 AOI EVI composite
# evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
# writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")

# ### aggregate raw 1m ISA to 30m
# ### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- too hard to monkey with, leave in its native state until the end
# ### raster-native approach required, file too large for data.table -- takes a long time
# ### correct the 16 values to NA to get a proper aggregated mean value per cell (0 = not impervious, 1 = impervious)
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.30m.tif", format="GTiff", overwrite=T)
# isa.na.agg <- raster("processed/isa.30m.tif")
# plot(isa.na.agg)

# ### align isa and evi grids, crop by AOI
# ### call python script for resampling 30m ISA to EVI grid
# pyth.path = './Rscripts/AOIISA_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)

isa.r <- raster("processed/isa30m_evigrd.tif")
evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83.tif") ## this is the July 2010-2012 AOI EVI composite
# AOI <- readOGR(dsn="E:/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/AOI_simple_NAD83UTM19N.shp", layer="AOI_simple_NAD83UTM19N")
# plot(isa.na.agg); plot(AOI, add=T)
# plot(evi.r); plot(AOI, add=T)
# AOI.r <- rasterize(AOI, evi.r)
# writeRaster(AOI.r, filename="processed/AOI.r.tif", format="GTiff", overwrite=T)
AOI.r <- raster("processed/AOI.r.tif")

# ### prep LULC --> EVI grid
# pyth.path = './Rscripts/LULC_EVIgrid.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)
LULC.r <- raster("processed/LULC30m.tif")

LULC.r <- extend(LULC.r, isa.r)
evi.r <- crop(evi.r, isa.r)
AOI.r <- crop(AOI.r, isa.r)
evi.r <- mask(evi.r, mask = AOI.r)

evi.stack <- stack(evi.r, isa.r, LULC.r, AOI.r)
writeRaster(evi.stack, filename="processed/EVI30m.stack.tif", format="GTiff", overwrite=T)






### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207009_ne_19_h_20160706/m_4207009_ne_19_h_20160706.tif")
naip.dat <- as.data.table(as.data.frame(naip.test))
names(naip.dat) <- c("band1", "band2", "band3", "band4")
naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
ndvi.r <- raster(naip.test[[1]])
ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)


#### Test Code: Retrieve EVI from 4-band Quickbird imagery
#### the below raw pan-sharpened data is 60cm, but is reported XXXXX integer -- assume this is corrected/orthorectified reflectance in integer form?
qb <- stack("Z:/Ultra1/Data/Remote_Sensing/Quickbird/Mosaic_pansharp/Boston_pan.tif")
names(qb) <- c("blue", "green", "red", "nir")
qb.evi <- function(x,y,z,filename){ ### enter as red blue nir
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=T, format="GTiff")
  for(i in 1:bs$n){
    red <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## red band
    red <- red/100000
    blue <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## blue band
    blue <- blue/100000
    nir <- getValues(z, row=bs$row[i], nrows=bs$nrows[i]) ## nir band
    nir <- nir/100000
    ## kill values that are all 0 (=NA)
    red.na <- red
    red.na[red==0 & blue==0 & nir==0] <- NA
    blue.na <- blue
    blue.na[red==0 & blue==0 & nir==0] <- NA
    nir.na <- nir
    nir.na[red==0 & blue==0 & nir==0] <- NA
    evi.calc <- 2.5*((nir.na-red.na)/(nir.na+6*red.na-7.5*blue.na+1))
    out <- writeValues(out, evi.calc, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
s <- qb.evi(qb[["red"]], qb[["blue"]], qb[["nir"]], filename="processed/bos.evi.1m.tif")

## this seems to give reasonable EVI values but seem pretty low even in dense canopy (~0.25 where NDVI values show 0.6)
## also there is a difference in registration (several meters) between this file and the other 1m boston rasters.
evi.bos <- raster("processed/bos.evi.1m.tif")

bos.stack <- stack("processed/bos.stack.1m.tif")
