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

