library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(data.table)


bounds <- c("UFFLA", "UFFLB", "R18", "R19", "R20")

n <- stack("data/NAIP/m_4207033_nw_19_h_20160706/m_4207033_nw_19_h_20160706.tif")
bound <- "UFFLA" ## northern plot in Lynn
ln30 <- readOGR(dsn = paste("processed/zones/", bound, "_30m.shp", sep=""),
                layer=paste(bound,"_30m", sep=""))
ln30@data$include <- 1
nn <- crop(n, extent(ln30))
plot(nn, main=paste(bound))
ln30.r <- rasterize(ln30[ln30@data$include], nn)
ln20 <- readOGR(dsn = paste("processed/zones/", bound, "_20m.shp", sep=""),
                layer=paste(bound,"_20m", sep=""))
ln20@data$include <- 1
ln20.r <- rasterize(ln20[ln20@data$include], nn)
ln10 <- readOGR(dsn = paste("processed/zones/", bound, "_10m.shp", sep=""),
                layer=paste(bound,"_10m", sep=""))
ln10@data$include <- 1
ln10.r <- rasterize(ln10[ln10@data$include], nn)



### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
nn.dat <- as.data.table(as.data.frame(nn))
names(nn.dat) <- c("band1", "band2", "band3", "band4")
nn.dat[, ndvi:=((band4-band1)/(band4+band1))]
nn.ndvi <- raster(nn[[1]])
nn.ndvi <- setValues(x = nn.ndvi, values=nn.dat[,ndvi])
plot(nn.ndvi)
plot(ln30, add=T)

gob <- stack(nn.ndvi, ln30.r, ln20.r, ln10.r)
gob.dat <- as.data.table(as.data.frame(gob))
names(gob.dat) <- c("ndvi", "buff30", "buff20", "buff10")
gob.dat$buff30 <- as.numeric(as.character(gob.dat$buff30))
gob.dat$buff20 <- as.numeric(as.character(gob.dat$buff20))
gob.dat$buff10 <- as.numeric(as.character(gob.dat$buff10))
gob.dat[is.na(buff30), buff30:=0]
gob.dat[is.na(buff20), buff20:=0]
gob.dat[is.na(buff10), buff10:=0]
gob.dat[buff30==1, median(ndvi)] ## total ndvi
gob.dat[buff30==1 & buff20==0 & buff10==0, median(ndvi)] # 30m NDVI
gob.dat[buff30==1 & buff20==1 & buff10==0, median(ndvi)] # 20m NDVI
gob.dat[buff30==1 & buff20==1 & buff10==1, median(ndvi)] # 10m NDVI
par(mfrow=c(1,3))
hist(gob.dat[buff30==1 & buff20==0 & buff10==0, ndvi]) # 30m NDVI
hist(gob.dat[buff30==1 & buff20==1 & buff10==0,ndvi]) # 20m NDVI
hist(gob.dat[buff30==1 & buff20==1 & buff10==1,ndvi]) # 10m NDVI
