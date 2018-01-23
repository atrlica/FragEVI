### process data layers to extract fragmentation per grid cell
library(raster)
library(rgeos)
library(rgdal)
library(sp)
library(data.table)
library(ff)

### read in AOI for crs description
AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
master.crs <- crs(AOI)

### load EVI composite for grid and process for ISA grid
evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif")
evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)

evi.AOI <- crop(evi.r, extent(AOI))
writeRaster(evi.AOI, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83_AOI.tif", format="GTiff", overwrite=T)
## need to get evi.AOI into the projection of the 1m ISA -- too hard to reproject the 1m data


### get ISA fraction per grid (~Urban intensity?)
isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
strip <- 0:15
tmp0 <- getValues(isa, row=seq(from=1, to=16), nrow=16)
tmp0[tmp0==16] <- NA
tmp0.r <- raster(matrix(data = tmp0, ncol=ncol(isa)))
for(s in 1:15){
  a <- ((strip[s]*16)+1)
  b <- ((strip[s]+1)*16)
  tmp <- getValues(isa, row=seq(from=a, to=b), nrow=16)
  tmp[tmp==16] <- NA
  assign(paste("tmp", s, sep=""), value = raster(matrix(data = tmp, ncol=ncol(isa))))
  
  isa <- setValues(isa, )
  print(paste("strip", s, sep=" "))
}
isa.dat <- as.big.matrix(isa)
isa.dat <- ff(vmode="double", dim=c(ncell(isa), 1), filename = "processed/isa.ffdata")
isa.dat[isa.dat==16, isa.dat:=NA]


