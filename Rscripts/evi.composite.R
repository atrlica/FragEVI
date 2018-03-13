### process landsat scenes to produce NDVI/EVI composites
library(raster)
library(data.table)

## be sure to set correct local director
write.files <- FALSE ### force raster rewrite

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/landsat/ARDsr_030006_072005-072013/"  ### for tile 030/006 (~p12r31)
# set.dir <- "data/landsat/ARDsr_030005_072005-072013/"  ### for tile 030/005 (~p12r30)
files <- list.files(set.dir)
files <- files[!grepl(files, pattern="tar")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 16,19) %in% scn.yr]

### process each file down to a QA-filtered EVI/NDVI stack
for(f in 1:length(files)){
  dir <- paste(set.dir, files[f], sep="")
  scn.ID <- paste(substr(files[f], 1,4), substr(files[f], 9, 14), substr(files[f], 16, 23), sep="_")
  bands <- list.files(dir)
  bands.sr <- bands[grep(bands, pattern ="SRB")]
  bands.qa <- bands[grep(bands, pattern="PIXELQA")]
  
  if(write.files | !file.exists(paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""))){
    
    ## load up the band stack
    temp <- stack(paste(dir, bands.sr[1], sep="/"))
    for(g in 2:length(bands.sr)){
      temp[[g]] <- raster(paste(dir, bands.sr[g], sep="/"))
      }
    
    ### convert to data.table and QA filter
    temp.dat <- as.data.table(as.data.frame(temp))
    qa.dat <- as.data.table(as.data.frame(raster(paste(dir, bands.qa, sep="/"))))
    temp.dat <- cbind(temp.dat, qa.dat)
    names(temp.dat) <- c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7", "QA")
    temp.dat[!QA%in%c(66,130), c("SR.B1", "SR.B2","SR.B3", "SR.B4", "SR.B5", "SR.B7"):=NA] ### kill values out of QA spec
    
    ## eliminate saturated values
    temp.dat[SR.B1>10000,] <- NA
    temp.dat[SR.B2>10000,] <- NA
    temp.dat[SR.B3>10000,] <- NA
    temp.dat[SR.B4>10000,] <- NA
    temp.dat[SR.B5>10000,] <- NA
    temp.dat[SR.B7>10000,] <- NA
    
    ## scale values
    temp.dat[,c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7"):= lapply(.SD, function(x) x*0.0001)]
    
    ## caluclate NDVI and EVI
    ndvi <- temp.dat[,(SR.B4-SR.B3)/(SR.B4+SR.B3)]
    evi <- temp.dat[, 2.5*((SR.B4-SR.B3)/(SR.B4+6*SR.B3-7.5*SR.B1+1))]
    
    ### rebuild rasters
    ndvi.r <- setValues(temp[[1]], ndvi)
    evi.r <- setValues(temp[[1]], evi)
    writeRaster(ndvi.r,
                filename=paste("processed/NDVI/ARD_", scn.ID, "_NDVI.tif", sep=""),
                format="GTiff", overwrite=T)
    writeRaster(evi.r,
                filename=paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""),
                format="GTiff", overwrite=T)
    }
  print(paste("finished with ", scn.ID, sep=""))
}

## load up rasters, strip, combine, and mosaic
r.files <- list.files("~/Desktop/FragEVI/processed/EVI/")
r.files <- r.files[!grepl(r.files, pattern=".xml")]
r.files <- r.files[substr(r.files, 17, 20)%in%scn.yr]
med.fun <- function(x){median(x, na.rm=T)}

### north raster combine and rebuild (row 005)
n <- r.files[substr(r.files, 13, 15)=="005"]
north <- raster(paste("~/Desktop/FragEVI/processed/EVI/", n[1], sep=""))
n.dat <- as.data.table(as.data.frame(north))
names(n.dat) <- "north.1"
for(d in 2:length(n)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", n[d], sep=""))))
  n.dat <- cbind(n.dat, tmp)
  names(n.dat)[d] <- paste("north.", d, sep="")
}
n.dat[, n.med:=apply(n.dat, 1, med.fun)]

north <- setValues(north, n.dat$n.med)
plot(north)

### south raster combine and rebuild (row 006)
s <- r.files[substr(r.files, 13, 15)=="006"]
south <- raster(paste("~/Desktop/FragEVI/processed/EVI/", s[1], sep=""))

s.dat <- as.data.table(as.data.frame(south))
names(s.dat) <- "south.1"
for(d in 2:length(s)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", s[d], sep=""))))
  s.dat <- cbind(s.dat, tmp)
  names(s.dat)[d] <- paste("south.", d, sep="")
}
s.dat[, s.med:=apply(s.dat, 1, med.fun)]

south <- setValues(south, s.dat$s.med)

### create mosaic and save
m <- mosaic(north, south, fun=median, na.rm=T,
            filename="processed/EVI/030005-6_2010-2012_EVI.m", format="GTiff")



##### make composites of MODIS EVI

## be sure to set correct local director
write.files <- FALSE ### force raster rewrite

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/landsat/ARDsr_030006_072005-072013/"  ### for tile 030/006 (~p12r31)
# set.dir <- "data/landsat/ARDsr_030005_072005-072013/"  ### for tile 030/005 (~p12r30)
files <- list.files(set.dir)
files <- files[!grepl(files, pattern="tar")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 16,19) %in% scn.yr]

### process each file down to a QA-filtered EVI/NDVI stack
for(f in 1:length(files)){
  dir <- paste(set.dir, files[f], sep="")
  scn.ID <- paste(substr(files[f], 1,4), substr(files[f], 9, 14), substr(files[f], 16, 23), sep="_")
  bands <- list.files(dir)
  bands.sr <- bands[grep(bands, pattern ="SRB")]
  bands.qa <- bands[grep(bands, pattern="PIXELQA")]
  
  if(write.files | !file.exists(paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""))){
    
    ## load up the band stack
    temp <- stack(paste(dir, bands.sr[1], sep="/"))
    for(g in 2:length(bands.sr)){
      temp[[g]] <- raster(paste(dir, bands.sr[g], sep="/"))
    }
    
    ### convert to data.table and QA filter
    temp.dat <- as.data.table(as.data.frame(temp))
    qa.dat <- as.data.table(as.data.frame(raster(paste(dir, bands.qa, sep="/"))))
    temp.dat <- cbind(temp.dat, qa.dat)
    names(temp.dat) <- c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7", "QA")
    temp.dat[!QA%in%c(66,130), c("SR.B1", "SR.B2","SR.B3", "SR.B4", "SR.B5", "SR.B7"):=NA] ### kill values out of QA spec
    
    ## eliminate saturated values
    temp.dat[SR.B1>10000,] <- NA
    temp.dat[SR.B2>10000,] <- NA
    temp.dat[SR.B3>10000,] <- NA
    temp.dat[SR.B4>10000,] <- NA
    temp.dat[SR.B5>10000,] <- NA
    temp.dat[SR.B7>10000,] <- NA
    
    ## scale values
    temp.dat[,c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7"):= lapply(.SD, function(x) x*0.0001)]
    
    ## caluclate NDVI and EVI
    ndvi <- temp.dat[,(SR.B4-SR.B3)/(SR.B4+SR.B3)]
    evi <- temp.dat[, 2.5*((SR.B4-SR.B3)/(SR.B4+6*SR.B3-7.5*SR.B1+1))]
    
    ### rebuild rasters
    ndvi.r <- setValues(temp[[1]], ndvi)
    evi.r <- setValues(temp[[1]], evi)
    writeRaster(ndvi.r,
                filename=paste("processed/NDVI/ARD_", scn.ID, "_NDVI.tif", sep=""),
                format="GTiff", overwrite=T)
    writeRaster(evi.r,
                filename=paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""),
                format="GTiff", overwrite=T)
  }
  print(paste("finished with ", scn.ID, sep=""))
}

## load up rasters, strip, combine, and mosaic
r.files <- list.files("~/Desktop/FragEVI/processed/EVI/")
r.files <- r.files[!grepl(r.files, pattern=".xml")]
r.files <- r.files[substr(r.files, 17, 20)%in%scn.yr]
med.fun <- function(x){median(x, na.rm=T)}

### north raster combine and rebuild (row 005)
n <- r.files[substr(r.files, 13, 15)=="005"]
north <- raster(paste("~/Desktop/FragEVI/processed/EVI/", n[1], sep=""))
n.dat <- as.data.table(as.data.frame(north))
names(n.dat) <- "north.1"
for(d in 2:length(n)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", n[d], sep=""))))
  n.dat <- cbind(n.dat, tmp)
  names(n.dat)[d] <- paste("north.", d, sep="")
}
n.dat[, n.med:=apply(n.dat, 1, med.fun)]

north <- setValues(north, n.dat$n.med)
plot(north)

### south raster combine and rebuild (row 006)
s <- r.files[substr(r.files, 13, 15)=="006"]
south <- raster(paste("~/Desktop/FragEVI/processed/EVI/", s[1], sep=""))

s.dat <- as.data.table(as.data.frame(south))
names(s.dat) <- "south.1"
for(d in 2:length(s)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", s[d], sep=""))))
  s.dat <- cbind(s.dat, tmp)
  names(s.dat)[d] <- paste("south.", d, sep="")
}
s.dat[, s.med:=apply(s.dat, 1, med.fun)]

south <- setValues(south, s.dat$s.med)

### create mosaic and save
m <- mosaic(north, south, fun=median, na.rm=T,
            filename="processed/EVI/030005-6_2010-2012_EVI.m", format="GTiff")



##### make composites of MODIS EVI



## be sure to set correct local director
write.files <- FALSE ### force raster rewrite

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/landsat/ARDsr_030006_072005-072013/"  ### for tile 030/006 (~p12r31)
# set.dir <- "data/landsat/ARDsr_030005_072005-072013/"  ### for tile 030/005 (~p12r30)
files <- list.files(set.dir)
files <- files[!grepl(files, pattern="tar")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 16,19) %in% scn.yr]

### process each file down to a QA-filtered EVI/NDVI stack
for(f in 1:length(files)){
  dir <- paste(set.dir, files[f], sep="")
  scn.ID <- paste(substr(files[f], 1,4), substr(files[f], 9, 14), substr(files[f], 16, 23), sep="_")
  bands <- list.files(dir)
  bands.sr <- bands[grep(bands, pattern ="SRB")]
  bands.qa <- bands[grep(bands, pattern="PIXELQA")]
  
  if(write.files | !file.exists(paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""))){
    
    ## load up the band stack
    temp <- stack(paste(dir, bands.sr[1], sep="/"))
    for(g in 2:length(bands.sr)){
      temp[[g]] <- raster(paste(dir, bands.sr[g], sep="/"))
    }
    
    ### convert to data.table and QA filter
    temp.dat <- as.data.table(as.data.frame(temp))
    qa.dat <- as.data.table(as.data.frame(raster(paste(dir, bands.qa, sep="/"))))
    temp.dat <- cbind(temp.dat, qa.dat)
    names(temp.dat) <- c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7", "QA")
    temp.dat[!QA%in%c(66,130), c("SR.B1", "SR.B2","SR.B3", "SR.B4", "SR.B5", "SR.B7"):=NA] ### kill values out of QA spec
    
    ## eliminate saturated values
    temp.dat[SR.B1>10000,] <- NA
    temp.dat[SR.B2>10000,] <- NA
    temp.dat[SR.B3>10000,] <- NA
    temp.dat[SR.B4>10000,] <- NA
    temp.dat[SR.B5>10000,] <- NA
    temp.dat[SR.B7>10000,] <- NA
    
    ## scale values
    temp.dat[,c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7"):= lapply(.SD, function(x) x*0.0001)]
    
    ## caluclate NDVI and EVI
    ndvi <- temp.dat[,(SR.B4-SR.B3)/(SR.B4+SR.B3)]
    evi <- temp.dat[, 2.5*((SR.B4-SR.B3)/(SR.B4+6*SR.B3-7.5*SR.B1+1))]
    
    ### rebuild rasters
    ndvi.r <- setValues(temp[[1]], ndvi)
    evi.r <- setValues(temp[[1]], evi)
    writeRaster(ndvi.r,
                filename=paste("processed/NDVI/ARD_", scn.ID, "_NDVI.tif", sep=""),
                format="GTiff", overwrite=T)
    writeRaster(evi.r,
                filename=paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""),
                format="GTiff", overwrite=T)
  }
  print(paste("finished with ", scn.ID, sep=""))
}

## load up rasters, strip, combine, and mosaic
r.files <- list.files("~/Desktop/FragEVI/processed/EVI/")
r.files <- r.files[!grepl(r.files, pattern=".xml")]
r.files <- r.files[substr(r.files, 17, 20)%in%scn.yr]
med.fun <- function(x){median(x, na.rm=T)}

### north raster combine and rebuild (row 005)
n <- r.files[substr(r.files, 13, 15)=="005"]
north <- raster(paste("~/Desktop/FragEVI/processed/EVI/", n[1], sep=""))
n.dat <- as.data.table(as.data.frame(north))
names(n.dat) <- "north.1"
for(d in 2:length(n)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", n[d], sep=""))))
  n.dat <- cbind(n.dat, tmp)
  names(n.dat)[d] <- paste("north.", d, sep="")
}
n.dat[, n.med:=apply(n.dat, 1, med.fun)]

north <- setValues(north, n.dat$n.med)
plot(north)

### south raster combine and rebuild (row 006)
s <- r.files[substr(r.files, 13, 15)=="006"]
south <- raster(paste("~/Desktop/FragEVI/processed/EVI/", s[1], sep=""))

s.dat <- as.data.table(as.data.frame(south))
names(s.dat) <- "south.1"
for(d in 2:length(s)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", s[d], sep=""))))
  s.dat <- cbind(s.dat, tmp)
  names(s.dat)[d] <- paste("south.", d, sep="")
}
s.dat[, s.med:=apply(s.dat, 1, med.fun)]

south <- setValues(south, s.dat$s.med)

### create mosaic and save
m <- mosaic(north, south, fun=median, na.rm=T,
            filename="processed/EVI/030005-6_2010-2012_EVI.m", format="GTiff")



##### make composites of MODIS EVI

## be sure to set correct local director
write.files <- FALSE ### force raster rewrite

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/landsat/ARDsr_030006_072005-072013/"  ### for tile 030/006 (~p12r31)
# set.dir <- "data/landsat/ARDsr_030005_072005-072013/"  ### for tile 030/005 (~p12r30)
files <- list.files(set.dir)
files <- files[!grepl(files, pattern="tar")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 16,19) %in% scn.yr]

### process each file down to a QA-filtered EVI/NDVI stack
for(f in 1:length(files)){
  dir <- paste(set.dir, files[f], sep="")
  scn.ID <- paste(substr(files[f], 1,4), substr(files[f], 9, 14), substr(files[f], 16, 23), sep="_")
  bands <- list.files(dir)
  bands.sr <- bands[grep(bands, pattern ="SRB")]
  bands.qa <- bands[grep(bands, pattern="PIXELQA")]
  
  if(write.files | !file.exists(paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""))){
    
    ## load up the band stack
    temp <- stack(paste(dir, bands.sr[1], sep="/"))
    for(g in 2:length(bands.sr)){
      temp[[g]] <- raster(paste(dir, bands.sr[g], sep="/"))
    }
    
    ### convert to data.table and QA filter
    temp.dat <- as.data.table(as.data.frame(temp))
    qa.dat <- as.data.table(as.data.frame(raster(paste(dir, bands.qa, sep="/"))))
    temp.dat <- cbind(temp.dat, qa.dat)
    names(temp.dat) <- c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7", "QA")
    temp.dat[!QA%in%c(66,130), c("SR.B1", "SR.B2","SR.B3", "SR.B4", "SR.B5", "SR.B7"):=NA] ### kill values out of QA spec
    
    ## eliminate saturated values
    temp.dat[SR.B1>10000,] <- NA
    temp.dat[SR.B2>10000,] <- NA
    temp.dat[SR.B3>10000,] <- NA
    temp.dat[SR.B4>10000,] <- NA
    temp.dat[SR.B5>10000,] <- NA
    temp.dat[SR.B7>10000,] <- NA
    
    ## scale values
    temp.dat[,c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7"):= lapply(.SD, function(x) x*0.0001)]
    
    ## caluclate NDVI and EVI
    ndvi <- temp.dat[,(SR.B4-SR.B3)/(SR.B4+SR.B3)]
    evi <- temp.dat[, 2.5*((SR.B4-SR.B3)/(SR.B4+6*SR.B3-7.5*SR.B1+1))]
    
    ### rebuild rasters
    ndvi.r <- setValues(temp[[1]], ndvi)
    evi.r <- setValues(temp[[1]], evi)
    writeRaster(ndvi.r,
                filename=paste("processed/NDVI/ARD_", scn.ID, "_NDVI.tif", sep=""),
                format="GTiff", overwrite=T)
    writeRaster(evi.r,
                filename=paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""),
                format="GTiff", overwrite=T)
  }
  print(paste("finished with ", scn.ID, sep=""))
}

## load up rasters, strip, combine, and mosaic
r.files <- list.files("~/Desktop/FragEVI/processed/EVI/")
r.files <- r.files[!grepl(r.files, pattern=".xml")]
r.files <- r.files[substr(r.files, 17, 20)%in%scn.yr]
med.fun <- function(x){median(x, na.rm=T)}

### north raster combine and rebuild (row 005)
n <- r.files[substr(r.files, 13, 15)=="005"]
north <- raster(paste("~/Desktop/FragEVI/processed/EVI/", n[1], sep=""))
n.dat <- as.data.table(as.data.frame(north))
names(n.dat) <- "north.1"
for(d in 2:length(n)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", n[d], sep=""))))
  n.dat <- cbind(n.dat, tmp)
  names(n.dat)[d] <- paste("north.", d, sep="")
}
n.dat[, n.med:=apply(n.dat, 1, med.fun)]

north <- setValues(north, n.dat$n.med)
plot(north)

### south raster combine and rebuild (row 006)
s <- r.files[substr(r.files, 13, 15)=="006"]
south <- raster(paste("~/Desktop/FragEVI/processed/EVI/", s[1], sep=""))

s.dat <- as.data.table(as.data.frame(south))
names(s.dat) <- "south.1"
for(d in 2:length(s)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", s[d], sep=""))))
  s.dat <- cbind(s.dat, tmp)
  names(s.dat)[d] <- paste("south.", d, sep="")
}
s.dat[, s.med:=apply(s.dat, 1, med.fun)]

south <- setValues(south, s.dat$s.med)

### create mosaic and save
m <- mosaic(north, south, fun=median, na.rm=T,
            filename="processed/EVI/030005-6_2010-2012_EVI.m", format="GTiff")



##### make composites of MODIS EVI

## be sure to set correct local director
write.files <- FALSE ### force raster rewrite

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/landsat/ARDsr_030006_072005-072013/"  ### for tile 030/006 (~p12r31)
# set.dir <- "data/landsat/ARDsr_030005_072005-072013/"  ### for tile 030/005 (~p12r30)
files <- list.files(set.dir)
files <- files[!grepl(files, pattern="tar")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 16,19) %in% scn.yr]

### process each file down to a QA-filtered EVI/NDVI stack
for(f in 1:length(files)){
  dir <- paste(set.dir, files[f], sep="")
  scn.ID <- paste(substr(files[f], 1,4), substr(files[f], 9, 14), substr(files[f], 16, 23), sep="_")
  bands <- list.files(dir)
  bands.sr <- bands[grep(bands, pattern ="SRB")]
  bands.qa <- bands[grep(bands, pattern="PIXELQA")]
  
  if(write.files | !file.exists(paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""))){
    
    ## load up the band stack
    temp <- stack(paste(dir, bands.sr[1], sep="/"))
    for(g in 2:length(bands.sr)){
      temp[[g]] <- raster(paste(dir, bands.sr[g], sep="/"))
    }
    
    ### convert to data.table and QA filter
    temp.dat <- as.data.table(as.data.frame(temp))
    qa.dat <- as.data.table(as.data.frame(raster(paste(dir, bands.qa, sep="/"))))
    temp.dat <- cbind(temp.dat, qa.dat)
    names(temp.dat) <- c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7", "QA")
    temp.dat[!QA%in%c(66,130), c("SR.B1", "SR.B2","SR.B3", "SR.B4", "SR.B5", "SR.B7"):=NA] ### kill values out of QA spec
    
    ## eliminate saturated values
    temp.dat[SR.B1>10000,] <- NA
    temp.dat[SR.B2>10000,] <- NA
    temp.dat[SR.B3>10000,] <- NA
    temp.dat[SR.B4>10000,] <- NA
    temp.dat[SR.B5>10000,] <- NA
    temp.dat[SR.B7>10000,] <- NA
    
    ## scale values
    temp.dat[,c("SR.B1", "SR.B2", "SR.B3", "SR.B4", "SR.B5", "SR.B7"):= lapply(.SD, function(x) x*0.0001)]
    
    ## caluclate NDVI and EVI
    ndvi <- temp.dat[,(SR.B4-SR.B3)/(SR.B4+SR.B3)]
    evi <- temp.dat[, 2.5*((SR.B4-SR.B3)/(SR.B4+6*SR.B3-7.5*SR.B1+1))]
    
    ### rebuild rasters
    ndvi.r <- setValues(temp[[1]], ndvi)
    evi.r <- setValues(temp[[1]], evi)
    writeRaster(ndvi.r,
                filename=paste("processed/NDVI/ARD_", scn.ID, "_NDVI.tif", sep=""),
                format="GTiff", overwrite=T)
    writeRaster(evi.r,
                filename=paste("processed/EVI/ARD_", scn.ID, "_EVI.tif", sep=""),
                format="GTiff", overwrite=T)
  }
  print(paste("finished with ", scn.ID, sep=""))
}

## load up rasters, strip, combine, and mosaic
r.files <- list.files("~/Desktop/FragEVI/processed/EVI/")
r.files <- r.files[!grepl(r.files, pattern=".xml")]
r.files <- r.files[substr(r.files, 17, 20)%in%scn.yr]
med.fun <- function(x){median(x, na.rm=T)}

### north raster combine and rebuild (row 005)
n <- r.files[substr(r.files, 13, 15)=="005"]
north <- raster(paste("~/Desktop/FragEVI/processed/EVI/", n[1], sep=""))
n.dat <- as.data.table(as.data.frame(north))
names(n.dat) <- "north.1"
for(d in 2:length(n)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", n[d], sep=""))))
  n.dat <- cbind(n.dat, tmp)
  names(n.dat)[d] <- paste("north.", d, sep="")
}
n.dat[, n.med:=apply(n.dat, 1, med.fun)]

north <- setValues(north, n.dat$n.med)
plot(north)

### south raster combine and rebuild (row 006)
s <- r.files[substr(r.files, 13, 15)=="006"]
south <- raster(paste("~/Desktop/FragEVI/processed/EVI/", s[1], sep=""))

s.dat <- as.data.table(as.data.frame(south))
names(s.dat) <- "south.1"
for(d in 2:length(s)){
  tmp <- as.data.table(as.data.frame(raster(paste("~/Desktop/FragEVI/processed/EVI/", s[d], sep=""))))
  s.dat <- cbind(s.dat, tmp)
  names(s.dat)[d] <- paste("south.", d, sep="")
}
s.dat[, s.med:=apply(s.dat, 1, med.fun)]

south <- setValues(south, s.dat$s.med)

### create mosaic and save
m <- mosaic(north, south, fun=median, na.rm=T,
            filename="processed/EVI/030005-6_2010-2012_EVI.m", format="GTiff")


#########
##### make composites of MODIS EVI

## figure out which files are there, target scene years, and work each file separately
set.dir <- "data/MODIS/GTiff/"  ### from ORNL LPDAAC subset of the AOI
files <- list.files(set.dir)
files <- files[grepl(files, pattern="EVI")]
scn.yr <- c(2010, 2011, 2012)
files <- files[substr(files, 10,13) %in% scn.yr]
july.days <- seq(182, 212)
files <- files[substr(files, 14, 16) %in% july.days]
qa <- list.files(set.dir)
qa <- qa[grepl(qa, pattern="VI_Quality")]
first_k_bits <- function(int, k=16, reverse=T) {
  ## https://lpdaac.usgs.gov/products/modis_products_table/mod13q1 TABLE 2:
  ## MOD13Q1 VI Quality: "Bit 0 is the least significant (read bit words right to left)"
  integer_vector <- as.integer(intToBits(int))[1:k]
  if(reverse) integer_vector <- rev(integer_vector)
  return(paste(as.character(integer_vector), collapse=""))
}

## clean up the rasters and package as a data frame for final compositing
cont <- list()
for(f in 1:length(files)){
  scn.ID <- paste(substr(files[f], 1,3), substr(files[f], 10, 16), substr(files[f], 56,58), sep="_")
  temp.r <- raster(paste(set.dir, files[f], sep=""))
  q <- qa[grep(substr(qa, 10, 16), pattern = substr(scn.ID, 5, 11))]
  qa.r <- raster(paste(set.dir, q, sep=""))

  ### pick apart the 16 bit unsigned integer for QA
  qa.dat <- as.vector(qa.r)
  df <- data.frame(bits=sapply(qa.dat, function(x) first_k_bits(x, k=16, reverse=T)))
  df <- cbind(qa.dat, df)
  
  ### pare down the acceptable pixels based on parsing the QA strings (integers --> 16 bit binary coding)
  hit <- unique(df$qa.dat[substr(df$bits, 15,16) %in% c("00", "01")]) ## VI produced good quality, produced but check other QA
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 11, 14) %in% c("0000",
                                                                      "0001")])] ## VI quality classes, 2 highest
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 9, 10) %in% c("01",## Aerosol load low
                                                                    "10")])] ## medium load
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 8, 8) %in% c("0")])] ## no adjacent clouds
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 6, 6) %in% c("0")])] ## no mixed clouds
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 3, 5) %in% c("001", # land only (too restrictive land only), ephemeral water, 
                                                                    "100",
                                                                    "010",
                                                                    "101")])] # ocean/shoreline
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 2, 2) %in% c("0")])] ## no snow/ice
  hit <- hit[hit %in% unique(df$qa.dat[substr(df$bits, 1, 1) %in% c("0")])] ## no shadow

  temp.dat <- as.data.table(as.data.frame(temp.r))
  temp.dat <- cbind(temp.dat, qa.dat)
  names(temp.dat) <- c("evi", "qa")

  ## eliminate invalid values
  temp.dat[evi>10000 | evi<(-2000), evi:=NA]
  test <- raster(temp.r)
  test <- setValues(test, temp.dat$evi)
  plot(test)
  temp.dat[!(qa %in% hit), evi:=NA]
  test <- raster(temp.r)
  test <- setValues(test, temp.dat$evi)
  plot(test)
  cont[[f]] <- temp.dat[,evi]
}
fra <- as.data.table(do.call(cbind, cont))
median.na <- function(x){median(x, na.rm=T)}
fra[,m.evi:=apply(fra, 1, FUN=median.na)]
m <- raster(temp.r)
m <- setValues(m, fra[,m.evi])
lims <- range(scn.yr)
writeRaster(m, filename=paste("processed/EVI/MOD_", lims[1], "-", lims[2], "_EVI.tif", sep=""), format="GTiff", overwrite=T)





