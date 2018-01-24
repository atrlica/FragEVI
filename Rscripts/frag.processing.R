### process data layers to extract fragmentation per grid cell
library(raster)
library(rgeos)
library(rgdal)
library(sp)
library(data.table)
library(ff)


### load EVI composite for grid and process for ISA grid
evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif")
evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)

evi.AOI <- crop(evi.r, extent(AOI))
writeRaster(evi.AOI, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83_AOI.tif", format="GTiff", overwrite=T)
evi.AOI <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83_AOI.tif")

## extent object is c(xmin, xmax, ymin, ymax)

### get ISA fraction per grid (~Urban intensity?)
### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- 
### too hard to monkey with, leave in its native state until the end

isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
fun <- function(x){
  x[x==16] <- NA
  return(x)
}
isa.na <- calc(isa, fun)

### aggregate to 30m 
isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
evi.isa <- projectRaster(from=evi.r, crs=crs(isa.na.agg), res=30, method="ngb")
evi.isa.cr <- crop(evi.isa, extent(isa.na.agg))


# isa.na.agg <- crop(isa.na.agg, extent(evi.isa.agg))
evi.isa <- crop(evi.isa,extent(isa.na.agg))

frag.stack <- stack(evi.isa, isa.na.agg)

AOI.isa <- spTransform(AOI, crs(isa.na.agg))
isa.na.agg <- crop(isa.na.agg, extent(AOI.isa))
isa.na.agg <- extend(isa.na.agg, extent(evi.AOI.isa))
evi.AOI.isa <- crop(evi.AOI.isa, extent(AOI.isa))
frag.stack <- stack(evi.AOI.isa, isa.na.agg)

isa.AOI <- projectRaster(from=isa.na.agg, crs=master.crs, res=30, method="ngb", filename="processed/isa_30mgrid.tif", format="GTiff", overwrite=T)


### read in AOI for crs description
AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
master.crs <- crs(AOI)


### add EVI, ISA, and LULC to a stack


plot(isa.AOI)
plot(AOI, add=T)
plot(evi.AOI)
plot(AOI, add=T)
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


