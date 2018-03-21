library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(data.table)
library(tidyr)
library(reshape2)
library(ggplot2)

### Test Code: Process NAIP 1m CIR data to NDVI
### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)


## pull up guidance data
guide <- read.csv("docs/MA_Edge_Coordinates.csv", stringsAsFactors = F)

bounds <- c("UFFLA", "UFFLB", "R18", "R19", "R20") # list of plots to look at
res <- data.frame()
for(b in 1:length(bounds)){
  out <- vector()
  loc <- guide[guide$Plot_ID==bounds[b], "NAIP.scene"]
  loc1 <- substr(loc, start=1, stop=12)
  loc2 <- substr(loc, start=13, stop=14)
  loc3 <- substr(loc, start=17, stop=24)
  loc.f <- paste(loc1, loc2, "h", loc3, sep="_")
  n <- stack(paste("data/NAIP/", loc.f, "/", loc.f, ".tif", sep=""))
  if(substr(bounds[b], 1, 1)=="R"){ ## rope plots go in 50m
    whole <- readOGR(dsn = paste("processed/zones/", bounds[b], "_50m.shp", sep=""),
                     layer=paste(bounds[b],"_50m", sep=""))
    whole@data$include <- 1
    nn <- crop(n, extent(whole))
    whole.r <- rasterize(whole[whole@data$include], nn, background=0)
    buffs <- as.vector(whole.r)
    for(d in c(10,20,30,40)){
      tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_", d, "m.shp", sep=""),
                                 layer=paste(bounds[b], "_", d, "m", sep=""))
      tmp@data$include <- 1
      tmp.r <- rasterize(tmp[tmp@data$include], nn, background=0)
      tmp.dat <- as.vector(tmp.r)
      buffs <- cbind(buffs, tmp.dat)
    }
    buffs.dat <- as.data.table(as.data.frame(buffs))
    names(buffs.dat) <- c("buff50", "buff10", "buff20", "buff30", "buff40")
    ## calculate ndvi
    nn.dat <- as.data.table(as.data.frame(nn))
    names(nn.dat) <- c("band1", "band2", "band3", "band4")
    nn.dat[, ndvi:=((band4-band1)/(band4+band1))]
    nn.ndvi <- raster(nn[[1]])
    nn.ndvi <- setValues(x = nn.ndvi, values=nn.dat[,ndvi])
    t <- cbind(nn.dat[,ndvi], buffs.dat)
    names(t)[1] <- c("ndvi")
    out <- t[buff50==1, mean(ndvi)] ## total ndvi
    out <- c(out, t[buff50==1 & buff40==0, mean(ndvi)]) # 50m NDVI
    out <- c(out, t[buff40==1 & buff30==0, mean(ndvi)]) # 40m NDVI
    out <- c(out, t[buff30==1 & buff20==0, mean(ndvi)]) # 30m NDVI
    out <- c(out, t[buff20==1 & buff10==0, mean(ndvi)]) # 20m NDVI
    out <- c(out, t[buff10==1, mean(ndvi)]) # 10m NDVI
  } else{ ## 30m Andy plots
    whole <- readOGR(dsn = paste("processed/zones/", bounds[b], "_30m.shp", sep=""),
                     layer=paste(bounds[b],"_30m", sep=""))
    whole@data$include <- 1
    nn <- crop(n, extent(whole))
    whole.r <- rasterize(whole[whole@data$include], nn, background=0)
    buffs <- as.vector(whole.r)
    for(d in c(10,20)){
      tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_", d, "m.shp", sep=""),
                     layer=paste(bounds[b], "_", d, "m", sep=""))
      tmp@data$include <- 1
      tmp.r <- rasterize(tmp[tmp@data$include], nn, background=0)
      tmp.dat <- as.vector(tmp.r)
      buffs <- cbind(buffs, tmp.dat)
    }
    buffs.dat <- as.data.table(as.data.frame(buffs))
    names(buffs.dat) <- c("buff30", "buff10", "buff20")
    ## calculate ndvi
    nn.dat <- as.data.table(as.data.frame(nn))
    names(nn.dat) <- c("band1", "band2", "band3", "band4")
    nn.dat[, ndvi:=((band4-band1)/(band4+band1))]
    nn.ndvi <- raster(nn[[1]])
    nn.ndvi <- setValues(x = nn.ndvi, values=nn.dat[,ndvi])
    t <- cbind(nn.dat[,ndvi], buffs.dat)
    names(t)[1] <- c("ndvi")
    out <- t[buff30==1, mean(ndvi)] ## total ndvi
    out <- c(out, NA) # 50m NDVI
    out <- c(out, NA) # 40m NDVI
    out <- c(out, t[buff30==1 & buff20==0, mean(ndvi)]) # 30m NDVI
    out <- c(out, t[buff20==1 & buff10==0, mean(ndvi)]) # 20m NDVI
    out <- c(out, t[buff10==1, mean(ndvi)]) # 10m NDVI
  }
  res <- rbind(res, out)
}
res <- cbind(res, bounds)
names(res) <- c("NDVI.tot", "50m.NDVI", "40m.NDVI", "30m.NDVI", "20m.NDVI", "10m.NDVI", "plot")

res.long <- melt(res, id.vars=c("plot"))
res.long$dist <- rep(10, dim(res.long)[1])
res.long$dist[res.long$variable=="50m.NDVI"] <- 50
res.long$dist[res.long$variable=="40m.NDVI"] <- 40
res.long$dist[res.long$variable=="30m.NDVI"] <- 30
res.long$dist[res.long$variable=="20m.NDVI"] <- 20
res.long <- res.long[res.long$variable!="NDVI.tot",]

ggplot(data = res.long, aes(x = dist, y = value, colour = plot)) +       
  geom_line() + geom_point()+
  labs(title="NDVI vs. edge position, 1m NAIP", x="Dist from edge (m)", y="Mean NDVI")
  

