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


### Process NAIP 1m data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207009_ne_19_h_20160706/m_4207009_ne_19_h_20160706.tif")
naip.dat <- as.data.table(as.data.frame(naip.test))
names(naip.dat) <- c("band1", "band2", "band3", "band4")
naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
ndvi.r <- raster(naip.test[[1]])
ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)
