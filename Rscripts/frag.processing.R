### process data layers to extract fragmentation per grid cell
library(raster)

## be sure to set correct local director
setwd("~/Desktop/FragEVI/") ## for laptop
# setwd("E:/FragEVI/") ## for desktop

### load a landsat 30 m stack, use as template
rast.dir <- "LE070120312010070501T1-SC20180118171143"
band.list <- list.files(paste("data/", rast.dir, sep=""))
band.list.sr <- band.list[grep(band.list, pattern = "band")]
band.list.qa <- band.list[grep(band.list, pattern="pixel_qa")]

gim <- paste("data", rast.dir, c(band.list.sr, band.list.qa), sep="/")
temp <- stack(gim[1])
for(g in 2:length(gim)){
  temp[[g]] <- raster(gim[g])
  temp[[g]]
}


