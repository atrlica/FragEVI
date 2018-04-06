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

## step 0: Prep the NAIP scenes that cover the field plots
guide <- read.csv("docs/MA_Edge_Coordinates.csv", stringsAsFactors = F)
loc <- unique(guide$NAIP.scene)
loc1 <- substr(loc, start=1, stop=12)
loc2 <- substr(loc, start=13, stop=14)
loc3 <- substr(loc, start=17, stop=24)
loc.f <- paste(loc1, loc2, "h", loc3, sep="_")
ndvi.calc <- function(x,y){
  n <- (x[]-y[])/(x[]+y[])
  return(n)
}
for(s in 1:length(loc.f)){
  tmp <- stack(paste("data/NAIP/", loc.f[s], "/", loc.f[s], ".tif", sep=""))
  tmp.n <- overlay(x = tmp[[4]], y = tmp[[1]], fun=ndvi.calc, 
                   filename=paste("data/NAIP/", loc.f[s], "/", loc.f[s], "ndvi.tif", sep=""),
                   format="GTiff", overwrite=T)
  
}

### step 1: extract data from NAIP ndvi for each plot at each distance polygon
bounds <- unique(guide$Plot_ID)
## cut out plots that are duplicated or that are to only 30m
bounds <- bounds[!bounds %in% c("UFFAA", "UFFAB", "UFFAC", "UFFAD", "UDDAE", "UFFAF", "UFFLA", "UFFLB", "R2", "R12")]
main <- data.frame()
for(b in 1:length(bounds)){
  ## load up corresponding NAIP scene
  loc <- guide$NAIP.scene[guide$Plot_ID==bounds[b]]
  loc1 <- substr(loc, start=1, stop=12)
  loc2 <- substr(loc, start=13, stop=14)
  loc3 <- substr(loc, start=17, stop=24)
  loc.f <- paste(loc1, loc2, "h", loc3, sep="_")
  r <- raster(paste("data/NAIP/", loc.f, "/", loc.f, "ndvi.tif", sep="") )

  ### make a data frame of the plot boundaries over the NDVI cropped
  tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_50m.shp", sep=""),
                 layer=paste(bounds[b], "_50m", sep=""), verbose=F)
  tmp <- spTransform(tmp, crs(r))
  r <- crop(r, tmp)
  tmp.dat <- as.vector(r)
  tmp@data$include <- 1
  buff.tmp <- rasterize(tmp[tmp@data$include], r, background=0)
  tmp.dat <- cbind(tmp.dat, as.vector(buff.tmp))
  
  dist <- c(10,20,30,40)
  for(d in 1:length(dist)){
    tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_", dist[d], "m.shp", sep=""),
                   layer=paste(bounds[b], "_", dist[d], "m", sep=""), verbose=F)
    tmp <- spTransform(tmp, crs(r))
    tmp@data$include <- 1
    buff.tmp <- rasterize(tmp[tmp@data$include], r, background=0)
    tmp.dat <- cbind(tmp.dat, as.vector(buff.tmp))
  }
  tmp.dat <- as.data.frame(tmp.dat)
  tmp.dat$plot <- rep(bounds[b], dim(tmp.dat)[1])
  names(tmp.dat) <- c("ndvi", "buff50", "buff10", "buff20", "buff30", "buff40", "plot")
  main <- rbind(main, tmp.dat)
  print(paste("processed plot", bounds[b]))
}

write.csv(main, "processed/plot.buffers.ndvi.csv")
dat <- read.csv("processed/plot.buffers.ndvi.csv")

## now do the (valid) other plots (30 m depth)
bounds <- c("UFFAA", "UFFAB",  "UFFAF", "UFFLA", "UFFLB")
main <- data.frame()
for(b in 1:length(bounds)){
  ## load up corresponding NAIP scene
  loc <- guide$NAIP.scene[guide$Plot_ID==bounds[b]]
  loc1 <- substr(loc, start=1, stop=12)
  loc2 <- substr(loc, start=13, stop=14)
  loc3 <- substr(loc, start=17, stop=24)
  loc.f <- paste(loc1, loc2, "h", loc3, sep="_")
  r <- raster(paste("data/NAIP/", loc.f, "/", loc.f, "ndvi.tif", sep="") )
  
  ### make a data frame of the plot boundaries over the NDVI cropped
  tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_30m.shp", sep=""),
                 layer=paste(bounds[b], "_30m", sep=""), verbose=F)
  tmp <- spTransform(tmp, crs(r))
  r <- crop(r, tmp)
  tmp.dat <- as.vector(r)
  tmp@data$include <- 1
  buff.tmp <- rasterize(tmp[tmp@data$include], r, background=0)
  tmp.dat <- cbind(tmp.dat, as.vector(buff.tmp))
  
  dist <- c(10,20)
  for(d in 1:length(dist)){
    tmp <- readOGR(dsn = paste("processed/zones/", bounds[b], "_", dist[d], "m.shp", sep=""),
                   layer=paste(bounds[b], "_", dist[d], "m", sep=""), verbose=F)
    tmp <- spTransform(tmp, crs(r))
    tmp@data$include <- 1
    buff.tmp <- rasterize(tmp[tmp@data$include], r, background=0)
    tmp.dat <- cbind(tmp.dat, as.vector(buff.tmp))
  }
  tmp.dat <- as.data.frame(tmp.dat)
  tmp.dat$plot <- rep(bounds[b], dim(tmp.dat)[1])
  names(tmp.dat) <- c("ndvi", "buff30", "buff10", "buff20", "plot")
  main <- rbind(main, tmp.dat)
  print(paste("processed plot", bounds[b]))
}

main$buff50 <- rep(0, dim(main)[1])
main$buff40 <- rep(0, dim(main)[1])
main <- main[,c("ndvi", "buff50", "buff10", "buff20", "buff30", "buff40", "plot")]
dat <- dat[,2:8]
dat <- rbind(dat, main)
write.csv(dat, "processed/plot.buffers.ndvi.csv")

### do some data analysis
ba1 <- read.csv("data/plots/Boston_Plots.csv", stringsAsFactors = F)
ba2 <- read.csv("data/plots/HF_Plots.csv")
ba3 <- read.csv("data/plots/Rope_Plots.csv")

### Boston Plots (0-30m)
ba1$Plot.ID <- paste(ba1$Project, ba1$Site, ba1$Plot, sep="")
sum.na <- function(x){sum(x, na.rm=T)}
ba1.t <- as.data.table(ba1)
ba1.t[,distB:=as.factor(distB)]
ba1.t[,Plot.ID:=as.factor(Plot.ID)]
ba1.t[,length(BA), by=Plot.ID] ## 40-100 trees per plot
ba1.t[,tot.BA:=sum(BA, na.rm=T), by=Plot.ID] ## total BA in each plot
ba1.t[,rel.BA:=sum(BA, na.rm=T)/(tot.BA/3), by=list(Plot.ID,distB)] ## Relative BA in each distance class

### HF plots (0-30m)
ba2.t <- as.data.table(ba2)
ba2.t$dfe <- as.factor(ba2.t$dfe)
ba2.t[,length(ba), by=Plot.ID] #50-80 trees per plot
ba2.t[,tot.BA:=sum(ba, na.rm=T), by=Plot.ID] ## total BA in each plot
ba2.t[,rel.BA:=sum(ba, na.rm=T)/(tot.BA/3), by=list(Plot.ID,dfe)] ## relative BA in each distance class

### Rope plots (0-50m)
ba3.t <- as.data.table(ba3)
ba3.t$SEGMENT <- as.factor(ba3.t$SEGMENT)
ba3.t[,length(ba.cm2), by=Plot.ID] #50-140 trees per plot
ba3.t[,tot.BA:=sum(ba.cm2, na.rm=T), by=Plot.ID] ## total BA in each plot
ba3.t[,n.bins:=length(unique(SEGMENT)), by=Plot.ID]
ba3.t[,rel.BA:=sum(ba.cm2, na.rm=T)/(tot.BA/n.bins), by=list(Plot.ID,SEGMENT)] ## relative BA in each distance class

ba1.s <- ba1.t[,list(Plot.ID, distB, tot.BA, rel.BA)]
ba1.s <- unique(ba1.s)
names(ba1.s) <- c("Plot.ID", "dist", "tot.BA", "rel.BA")
ba2.s <- ba2.t[,list(Plot.ID, dfe, tot.BA, rel.BA)]
ba2.s <- unique(ba2.s)
names(ba2.s) <- c("Plot.ID", "dist", "tot.BA", "rel.BA")
ba3.s <- ba3.t[,list(Plot.ID, SEGMENT, tot.BA, rel.BA)]
ba3.s <- unique(ba3.s)
names(ba3.s) <- c("Plot.ID", "dist", "tot.BA", "rel.BA")

ba.s <- rbind(ba1.s, ba2.s, ba3.s)

### now bind with the NDVI extracted over these plots
dat <- read.csv("processed/plot.buffers.ndvi.csv")
dat.t <- as.data.table(dat)
dat.t[buff50==1, dist:=50]
dat.t[dist==50, length(ndvi)]
dat.t[buff40==1, dist:=40]
dat.t[dist==50, length(ndvi)]
dat.t[buff30==1, dist:=30]
dat.t[dist==50, length(ndvi)]
dat.t[buff20==1, dist:=20]
dat.t[dist==50, length(ndvi)]
dat.t[buff10==1, dist:=10]
dat.t[dist==50, length(ndvi)]
dat.t[dist==40, length(ndvi)]
dat.t[dist==30, length(ndvi)]
dat.t[dist==20, length(ndvi)]
dat.t[dist==10, length(ndvi)]
dat.t[,dist:=as.factor(dist)]
dat.t <- dat.t[!is.na(dist),]
dat.t[,tot.ndvi:=mean(ndvi, na.rm=T), by=plot]
dat.t[,rel.ndvi:=mean(ndvi,na.rm=T)/tot.ndvi, by=list(plot, dist)]

dat.s <- dat.t[,list(plot, dist, tot.ndvi, rel.ndvi)]
dat.s <- unique(dat.s)
names(dat.s) <- c("Plot.ID", "dist", "tot.ndvi", "rel.ndvi")

plot.sum <- merge(ba.s, dat.s, by=c("Plot.ID", "dist"))
plot.sum$dist <- as.integer(as.character(plot.sum$dist))
class(plot.sum$dist)
ggplot(data = plot.sum, aes(x = dist, y = rel.ndvi, colour = Plot.ID)) +       
  geom_line() + geom_point()+
  labs(title="NDVI vs. edge position, 1m NAIP", x="Dist from edge (m)", y="NDVI vs. plot mean")
ggplot(data = plot.sum, aes(x = dist, y = rel.BA, colour = Plot.ID)) +       
  geom_line() + geom_point()+
  labs(title="Basal Area vs. edge position, 1m NAIP", x="Dist from edge (m)", y="BA vs. plot mean")



## For ANDY plots: process NAIP data

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
  if(substr(bounds[b], 1, 1) %in% "R"){ ## rope plots go in 50m
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
  

