library(raster)
library(data.table)

dat.r <- stack("processed/evi.isa.30m.tif")
names(dat.r) <-  c("evi", "isa")
dat <- as.data.table(as.data.frame(dat.r))

