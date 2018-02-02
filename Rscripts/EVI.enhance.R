library(raster)
library(data.table)

dat <- stack("E:/FragEVI/processed/EVI.ISA.30mstack.tif")
names(dat) <-  c("evi", "isa")
