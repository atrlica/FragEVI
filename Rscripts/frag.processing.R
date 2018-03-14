### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)


#########
### Boston high-res data
### combine 2.4m NDVI and canopy map to produce 1m veg classification map
### use Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map

# ### step 1: 1m Canopy presence/absence map 
# bos.can <- raster("data/dataverse_files/bostoncanopy_1m.tif")

# ### set no canopy (0) values to NA
# bos.can.dat <- as.data.table(as.data.frame(bos.can))
# bos.can.dat[bostoncanopy_1m==0,] <- NA
# bos.can.na <- raster(bos.can)
# bos.can.na <- setValues(bos.can.na, bos.can.dat$bostoncanopy_1m)
# writeRaster(bos.can.na, "E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif", format="GTiff", overwrite=T, datatype="INT1U")
# bos.can.na <- raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m_na.tif")

# ### step 3 (must be performed on desktop): put NDVI.img through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
# pyth.path = './Rscripts/NDVI_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(paste("ArcPy working on NDVI resample: ", output))

# ### step 4: land cover classification for Boston using the 1m resampled NDVI (0 = barren, 1 = grass, 2 = canopy)
# bos.ndvi <- raster("data/NDVI/NDVI_1m_res_cangrid.tif")
# bos.can.na <- raster("data/dataverse_files/bostoncanopy_1m_na.tif")
# bos.ndvi <- crop(bos.ndvi, bos.can.na)
# cover.bl <- function(x, y, filename) { # x is canopy, y is ndvi
#   if(file.exists("processed/bos.cov.tif")){file.remove("processed/bos.cov.tif")}
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## canopy map
#     g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi map
#     cov <- v
#     cov[g>=0.2 & v==0] <- 1
#     cov[v==1] <- 2
#     out <- writeValues(out, cov, bs$row[i])
#     print(paste("finished 1m cover block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# bos.cov <- cover.bl(bos.can, bos.ndvi, filename="processed/bos.cov.tif")

# ### Processing 1m canopy map for edge class
# ### call python scripbt for identifying canopy edge distance (desktop only)
# pyth.path = './Rscripts/canopy_process.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)
# 
# ### read in edge rasters and correct classifications to produce raster of edge rings (>10, 10-20, 20-30, >30)
# bos.cov <- raster("processed/bos.cov.tif")
# ed1 <- raster("processed/nocan_10mbuff.tif")
# ed2 <- raster("processed/nocan_20mbuff.tif")
# ed3 <- raster("processed/nocan_30mbuff.tif")
# ed1 <- extend(ed1, bos.cov) # edge and cover on same grid already
# ed2 <- extend(ed2, bos.cov)
# ed3 <- extend(ed3, bos.cov)
# # edges <- stack(ed1, ed2, ed3, bos.cov) # just to make damn sure everything lines up
# 
# edges.bl <- function(x, y, filename) { # x is edge class, y is cover class
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
#     g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## cover map
#     v[v!=0] <- NA # kill any weird values that aren't coming from the nocan==0 buffer map
#     v[g!=2] <- NA # cancel edge ID for non-canopy
#     v[v==0] <- 1 # ID edge pixels as 1
#     v[v!=1] <- 0 ## set non-edge values to 0 to hold space
#     out <- writeValues(out, v, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# s <- edges.bl(ed1, bos.cov, filename="processed/edge10m.tif")
# t <- edges.bl(ed2, bos.cov, filename="processed/edge20m.tif")
# u <- edges.bl(ed3, bos.cov, filename="processed/edge30m.tif")

### package up all the 1m data for Boston and write to disk
# ed10 <- raster("processed/edge10m.tif")
# ed20 <- raster("processed/edge20m.tif")
# ed30 <- raster("processed/edge30m.tif")
# bos.cov <- raster("processed/bos.cov.tif")
# bos.can <- raster("data/dataverse_files/bostoncanopy_1m.tif")
# bos.ndvi <- raster("data/NDVI/NDVI_1m_res_cangrid.tif")
# bos.ndvi <- crop(bos.ndvi, bos.can)

# ### add in 1m ISA (go get it, crop it, then reproject)
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# towns <- readOGR(dsn = "/projectnb/buultra/atrlica/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# bos.bord <- towns[towns@data$TOWN=="BOSTON",]
# bos.bord <- spTransform(bos.bord, crs(isa))
# bos.isa <- crop(isa, bos.bord)
# writeRaster(bos.isa, filename="processed/bos_isa_1m.tif", format="GTiff", overwrite=T)
# # python.load("Rscripts/bosISA_resamp.py") ### test this
# bos.isa <- raster("processed/isa_cangrid.tif")

# 
# bos.stack <- stack(bos.can, bos.ndvi, bos.cov, bos.isa, ed10, ed20, ed30)
# writeRaster(bos.stack, "processed/boston/bos.stack.1m.tif", format="GTiff", overwrite=T)


##### work to get 1m boston data aggregated to 30m EVI grid
bos.stack <- stack("processed/boston/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")

### get Boston limits, raterize to EVI 30m grid
# towns <- readOGR(dsn = "F:/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# bos.AOI <- towns[towns@data$TOWN=="BOSTON",]
# bos.AOI <- bos.AOI[bos.AOI@data$SHAPE_AREA>1E07,] ## remove Harbor Islands
# bos.AOI <- spTransform(bos.AOI, crs(bos.stack))
# bos.AOI@data$include <- 1
# bos.AOI.r <- rasterize(bos.AOI[bos.AOI@data$include], evi.stack[[1]])
# bos.AOI.r <- crop(bos.AOI.r, bos.AOI)
# ba.dat <- as.data.table(as.data.frame(bos.AOI.r))
# ba.dat[!is.na(OBJECTID) | !is.na(OBJECTID.1), dogshit:=1] ## set all areas to same value
# bos.AOI.r <- setValues(bos.AOI.r, ba.dat$dogshit)
# 
# bos.evi <- crop(evi.stack, bos.AOI.r)
# bos.evi <- mask(bos.evi, bos.AOI.r)
# names(bos.evi) <- c("evi", "isa", "lulc", "AOI") ### the Boston city limits in 30m grid
# sum(getValues(bos.evi[["AOI"]]), na.rm=T) ## 138k pixels

##### Work boston 1m layers out into separate 0/1 classified files
### make 1m Boston AOI raster mask
# bos.AOI.1m <- rasterize(bos.AOI[bos.AOI@data$include], bos.stack[[1]])
# bos.AOI.1m <- crop(bos.AOI.1m, bos.AOI)
# bos.AOI.1mf <- raster(bos.AOI.1m)
# bos.AOI.1mf <- setValues(bos.AOI.1mf, values = rep(1, length=ncell(bos.AOI.1m)))
# bos.AOI.1mf <- mask(bos.AOI.1mf, bos.AOI.1m)
# writeRaster(bos.AOI.1mf, filename="processed/boston/bos.aoi.tif", format="GTiff", overwrite=T)

### isa reprocess -- use reregistered version that fits the barren layer better
## first get the ISA back onto the cover grid
pyth.path = './Rscripts/bosISA_resamp2.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
print(output)
### this is not working for reasons that are inexplicable March 10-2018
## fuk it work with the fucked up isa
# bos.stack <- stack("processed/boston/bos.stack.1m.tif")
# names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
# 
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# bos.isa <- bos.stack[["isa"]]
# bos.isa <- crop(bos.stack, bos.aoi)
# bos.isa <- mask(bos.isa, bos.aoi)
# writeRaster(bos.isa, filename="processed/boston/bos.isa.tif")
bos.isa <- raster("processed/boston/bos.isa.tif")

### get separate layers for each with cover=0/1 (vs. cover=0,1,2)
bos.stack <- stack("processed/boston/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")

grass.find <- function(x){
  x[x!=1] <- 0
  return(x)
}
grass.only <- calc(bos.stack[["cov"]], fun=grass.find, filename="processed/boston/bos.grass.tif", format="GTiff", overwrite=T)

# barr.find <- function(x){
#   x[x==0] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
# }
# barr.only <- calc(bos.stack[["cov"]], barr.find, filename="processed/boston/bos.barr.tif", format="GTiff", overwrite=T)

# can.find <- function(x){
#   x[x==2] <- 10
#   x[x!=10] <- 0
#   x[] <- x[]/10
#   return(x)
# }
# can <- calc(bos.stack[["cov"]], can.find, filename="processed/bos.can.tif", format="GTiff", overwrite=T)
# can.only <- raster("processed/bos.can.tif")
# # test if can_only == can
# can.test <- overlay(bos.stack[["can"]], can.only, fun=function(x,y){return(x-y)})
# plot(can.test) # yes they are identical

### grass needs correction due to artifacts
### lower res on NDVI data means that some spillover is happening into probable ISA next to canopy edges
### i.e. grass is showing up in v. small isolated patches also classed as ISA
### reclass grass==1 & isa==1 as grass==0
bos.grass <- raster("processed/boston/bos.grass.tif")
bos.isa <- raster("processed/boston/bos.isa.tif")
c <- brick(bos.grass, bos.isa)
grass.fix <- function(x,y){
  x[y==1] <- 0
  return(x)
}
bos.grass <- overlay(c, fun=grass.fix, filename="processed/boston/bos.grass.tif", format="GTiff", overwrite=T)
### this raster looks better, but misalignment in some neighborhoods between the cov and isa layers means a fair amount of back yards are getting cut

## identify non-impervious fraction of barren ## Feb 24 this is not quite working right
# nonimp.only <- overlay(barr.only, bos.stack[["isa"]], fun=function(x,y){return(x-y)}, filename="processed/boston/bos.nonimp_only.tif", overwrite=T, format="GTiff")
# nonimp.only <- raster("processed/boston/bos.nonimp.tif")
# nonimp.find <- function(x, filename.nonimpbarr, filename.vegisa) { # x is first-pass nonimp
#   out1 <- raster(x)
#   out2 <- raster(x)
#   bs <- blockSize(out1)
#   out1 <- writeStart(out1, filename.nonimpbarr, overwrite=TRUE, format="GTiff")
#   out2<- writeStart(out2, filename.vegisa, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## original raster
#     d <- rep(0, length(v)) ## create blank 0 raster
#     e <- d
#     d[v==1] <- 1
#     e[v==(-1)] <- 1
#     out1 <- writeValues(out1, d, bs$row[i])
#     out2 <- writeValues(out2, e, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out1 <- writeStop(out1)
#   out2 <- writeStop(out2)
#   return(out1)
#   return(out2)
# }
# nonimp.find(nonimp.only, filename.nonimpbarr="processed/boston/bos.nonimpbarr.tif", filename.vegisa = "processed/boston/bos.vegisa.tif")
# nonimpbarr.only <- raster("processed/boston/bos.nonimpbarr.tif")
# vegisa.only <- raster("processed/boston/bos.vegisa.tif")

### new approach to finding veg-over-isa and nonimpervious barren
bos.can <- raster("processed/boston/bos.can.tif")
bos.isa <- raster("processed/boston/bos.isa.tif")
bos.barr <- raster("processed/boston/bos.barr.tif")

a <- brick(bos.can, bos.isa)
vegisa.calc <- function(x, y){
  x[x==1 & y==1] <- 2 ## mark where there is isa & canopy
  x[y==0 & x==1] <- 0 ## erase where there is canopy but no isa
  x[x==2] <- 1 ## give vegisa marker value of 1
  return(x)
}
vegisa <- overlay(a, fun=vegisa.calc, filename="processed/boston/bos.vegisa.tif", format="GTiff", overwrite=T)

b <- brick(bos.isa, bos.barr)
nib.calc <- function(x,y){
  y[x==0 & y==1] <- 2 ## mark where it's barren but not paved
  y[x==1] <- 0 ## erase any paved area
  y[y==2] <- 1 ## give nonimpbarr marker value of 1
  return(y)
}
nonimpbarr <- overlay(b, fun=nib.calc, filename="processed/boston/bos.nonimpbarr.tif", format="GTiff", overwrite=T)
### both of these, like grass, vulnerable to getting screwed up by misregistration between ISA and barren layers

# ### make map of individual buffer rings
# bos.stack <- stack("processed/boston/bos.stack.1m.tif")
# names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
# bos.aoi <- raster("processed/boston/bos.aoi_only.tif")
# 
# ### first need to reclass the edge maps to 0/1 (are in 1/NA)
# fix.buff <- function(x){ # x is buffer, y is next smallest buffer
#   x[is.na(x)] <- 0
#   return(x)
# }
# bos.stack[["ed10"]] <- calc(bos.stack[["ed10"]], fun=fix.buff)
# bos.stack[["ed20"]] <- calc(bos.stack[["ed20"]], fun=fix.buff)
# bos.stack[["ed30"]] <- calc(bos.stack[["ed30"]], fun=fix.buff)
# ### these need masking
# bos.stack <- crop(bos.stack, bos.aoi)
# bos.stack <- mask(bos.stack, bos.aoi)
# writeRaster(bos.stack, filename="processed/boston/bos.stack.1m.tif", format="GTiff", overwrite=T)
# plot(bos.stack)
# 
# bos.stack <- stack("processed/boston/bos.stack.1m.tif")
# names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")
# 
# ### 10m ring is already done in ed10
# 
# ### 20m ring
# buff.20only <- overlay(bos.stack[["ed10"]], bos.stack[["ed20"]], fun=function(x,y){return(y-x)}, filename="processed/boston/bos.buff20_only.tif", format="GTiff", overwrite=T)
# 
# ### 30m ring
# buff.30only <- overlay(bos.stack[["ed20"]], bos.stack[["ed30"]], fun=function(x,y){return(y-x)}, filename="processed/boston/bos.buff30_only.tif", format="GTiff", overwrite=T)
# 
# ### interior
# buff.Intonly <- overlay(bos.stack[["ed30"]], bos.stack[["can"]], fun=function(x,y){return(y-x)}, filename="processed/boston/bos.buffInt_only.tif", format="GTiff", overwrite=T)
# 
# buffs.only <- stack(buff.20only, buff.30only, buff.Intonly)
# writeRaster(buffs.only, filename="processed/boston/bos.buffs_only.tif", format="GTiff", overwrite=T)
# 
####
### add LULC 1m raster info
# ### call python script for rasterize LULC in Boston to 1m canopy grid
# pyth.path = './Rscripts/LULC_bos_rast.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)

# bos.lulc <- raster("processed/boston/LU_bos_r1m.tif")
# bos.aoi <- raster("processed/boston/bos.aoi_only.tif")
# bos.lulc <- crop(bos.lulc, bos.aoi)
# bos.lulc <- mask(bos.lulc, bos.aoi)
# writeRaster(bos.lulc, filename="processed/boston/bos.lulc_only.tif", format="GTiff", overwrite=T)
# 
# ## decide on a collapsed LULC scheme to flag for area fraction
# lu.classnames <- c("forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other", "water")
# lu.forest <- c(3,37) # Forest, FWet
# lu.dev <- c(15,16,7,8,18,39,24,31,19,9,29) # Comm, Ind, PartRec, SpectRec, Transp., Junk, Util, PubInst, Waste, WBRec, Marina
# lu.hdres <- c(11,10) # HDResid., MFResid., 
# lu.medres <- c(12) # MDResid.
# lu.lowres <- c(13, 38) #LDResid., VLDResid.
# lu.lowveg <- c(4,1,2,6,5,26,14,25,34,23) # NFwet, PubInst, Crop, Past, Open, Mining, Golf, SWwet, SWBeach, Cem, CBog
# lu.other <- c(17,35,40,36) #Tran, Orch, Brush, Nurs
# lu.water <- c(20)
# lulc.tot <- list(lu.for, lu.dev, lu.hdres, lu.medres, lu.lowres, lu.lowveg, lu.other, lu.water)
# 
# ## flexible blockwise function for flagging different class groups of LULC
# lulc.flag <- function(x, lu.spot, filename) { # x is 1m lulc, lu.spot is list of LULC targeted
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
#     o <- rep(0, length(v))
#     o[v%in%lu.spot] <- 1
#     out <- writeValues(out, o, bs$row[i])
#     print(paste("finished block", i, "of", bs$n, "in", lu.classnames[l]))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# 
# ### loop through LULC class groups and create separately flagged LULC rasters
# bos.aoi <- raster("processed/boston/bos.aoi_only.tif")
# for(l in 1:length(lulc.tot)){
#   tmp <- do.call(lulc.flag, 
#           args = list(bos.lulc, 
#                       lulc.tot[[l]], 
#                       paste("processed/boston/bos.", lu.classnames[l], "_only.tif", sep="")))
#   print(paste("masking to Boston AOI"))
#   tmp <- mask(tmp, bos.aoi) 
#   writeRaster(tmp, filename=paste("processed/boston/bos.", lu.classnames[l], "_only.tif", sep=""), format="GTiff", overwrite=T)
# }
# 
# bos.forest <- raster("processed/boston/bos.forest.tif")
# bos.dev <- raster("processed/boston/bos.dev.tif")
# bos.hd.res <- raster("processed/boston/bos.hd.res.tif")
# bos.med.res <- raster("processed/boston/bos.med.res.tif")
# bos.low.res <- raster("processed/boston/bos.low.res.tif")
# bos.lowveg <- raster("processed/boston/bos.lowveg.tif")
# bos.other <- raster("processed/boston/bos.other.tif")
# bos.water <- raster("processed/boston/bos.water.tif")
# plot(bos.forest, main="forest")
# plot(bos.hd.res, main="HDres")
# plot(bos.med.res, main="MDres")
# plot(bos.low.res, main="LDRes")
# plot(bos.lowveg, main="LowVeg")
# plot(bos.other, main="other")
# plot(bos.water, main="water")
# 

### prep all individual raster layers (0/1) for aggregation all at once in arc
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.stack <- stack("processed/boston/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")

ndvi.only <- bos.stack[["ndvi"]]
can.only <- bos.stack[["can"]]
isa.only <- raster("processed/boston/bos.isa.tif")
# grass.only <- raster("processed/boston/bos.grass.tif")
# grass.only <- crop(grass.only, bos.aoi)
# grass.only <- mask(grass.only, bos.aoi, filename="processed/boston/bos.grass.tif", format="GTiff", overwrite=T)
grass.only <- raster("processed/boston/bos.grass.tif")
# barr.only <- raster("processed/boston/bos.barr.tif")
# barr.only <- crop(barr.only, bos.aoi)
# barr.only <- mask(barr.only, bos.aoi, filename="processed/boston/bos.barr.tif", format="GTiff", overwrite=T)
barr.only <- raster("processed/boston/bos.barr.tif")
# nonimp.only <- raster("processed/boston/bos.nonimp.tif")
# nonimp.only <- crop(nonimp.only, bos.aoi)
# nonimp.only <- mask(nonimp.only, bos.aoi, filename="processed/boston/bos.nonimp_only.tif", format="GTiff", overwrite=T)
# nonimpbarr.only <- raster("processed/boston/bos.nonimpbarr.tif")
# nonimpbarr.only <- crop(nonimpbarr.only, bos.aoi)
# nonimpbarr.only <- mask(nonimpbarr.only, bos.aoi, filename="processed/boston/bos.nonimpbarr.tif", format="GTiff", overwrite=T)
nonimpbarr.only <- raster("processed/boston/bos.nonimpbarr.tif")
# vegisa.only <- raster("processed/boston/bos.vegisa.tif")
# vegisa.only <- crop(vegisa.only, bos.aoi)
# vegisa.only <- mask(vegisa.only, bos.aoi, filename="processed/boston/bos.vegisa.tif", format="GTiff", overwrite=T)
vegisa.only <- raster("processed/boston/bos.vegisa.tif")
ed10.only <- bos.stack[["ed10"]]
ed20.only <- bos.stack[["ed20"]]
ed30.only <- bos.stack[["ed30"]]
buffs.only <- stack("processed/boston/bos.buffs.tif")
names(buffs.only) <- c("buff.20only","buff.30only", "buff.Intonly")
buff.20only <- buffs.only[["buff.20only"]]
buff.30only <- buffs.only[["buff.30only"]]
buff.Intonly <- buffs.only[["buff.Intonly"]]
bos.forest <- raster("processed/boston/bos.forest.tif")
bos.dev <- raster("processed/boston/bos.dev.tif")
bos.hd.res <- raster("processed/boston/bos.hd.res.tif")
bos.med.res <- raster("processed/boston/bos.med.res.tif")
bos.low.res <- raster("processed/boston/bos.low.res.tif")
bos.lowveg <- raster("processed/boston/bos.lowveg.tif")
bos.other <- raster("processed/boston/bos.other.tif")
bos.water <- raster("processed/boston/bos.water.tif")


## package up the boston areamarks stack for safekeeping, then split up for aggregation
### order for area marking: 
# ndvi, can, grass, barren, isa, nonimpbarr, vegisa, (1-7)
# ed10, ed20, ed30, buff.20only, buff.30only, buff.Intonly, (8-13)
# for, dev, hd.res, med.res, low.res, lowveg, other, water (9-16)
bos.stack <- stack(ndvi.only, can.only, grass.only, barr.only, isa.only, nonimpbarr.only, vegisa.only,
                   ed10.only, ed20.only, ed30.only, buff.20only, buff.30only, buff.Intonly,
                   bos.forest, bos.dev, bos.hd.res, bos.med.res, bos.low.res, bos.lowveg, bos.other, bos.water)
names(bos.stack) <- c("ndvi", "can", "grass", "barr", "isa", "nonimpbarr", "vegisa", 
                      "ed10", "ed20", "ed30", "buff.20only", "buff.30only", "buff.Intonly", 
                      "forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other", "water")

## make sure everything has been exported as single raster files
for(l in 1:nlayers(bos.stack)){
  if(!file.exists(paste("processed/boston/bos.", names(bos.stack)[l], ".tif", sep=""))){
    print(paste("writing ", "processed/boston/bos.", names(bos.stack)[l], ".tif", sep=""))
    writeRaster(bos.stack[[l]], 
                filename=paste("processed/boston/bos.", names(bos.stack)[l], ".tif", sep=""),
                format="GTiff", overwrite=T)
  }
  print("all files exist as independent rasters")
}

# ### make up the final raster stack and export
# bos.names <- c("ndvi", "can", "grass", "barr", "isa", "nonimpbarr", "vegisa", 
#   "ed10", "ed20", "ed30", "buff.20only", "buff.30only", "buff.Intonly", 
#   "forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other", "water")
# bos.stack.30m <- raster("processed/boston/bos.ndvi.tif")
# for(n in 2:length(bos.names)){
#   bos.stack.30m <- stack(bos.stack.30m, raster(paste("processed/boston/bos", bos.names[n], "tif", sep=".")))
# }
# writeRaster(bos.stack, filename="processed/bos.stack.1m.areamarks.tif", format="GTiff", overwrite=T)
# 

## get fractional area of each cover class per 30m aggregate cell
pyth.path = './Rscripts/bos1m_agg.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)

### shovel out individual files and compile as data frame
bos.sizes <- raster("processed/boston/bos.aoi30m.tif")
bos.30m <- list.files("processed/boston/")
bos.30m <- bos.30m[grepl(pattern = "30m", x=bos.30m)]
bos.30m <- bos.30m[!grepl(pattern=".xml", x=bos.30m)]
bos.30m <- bos.30m[!grepl(pattern=".tfw", x=bos.30m)]

tmp <- stack()

for(f in 1:length(bos.30m)){
  tmp <- stack(tmp, raster(paste("processed/boston/", bos.30m[f], sep="")))
  print(f)
}

fields <- sub("bos. *(.*?) *30m.tif.*", "\\1", bos.30m)
bos.dat <- as.data.table(as.data.frame(tmp))
names(bos.dat) <- fields

## combine with 30m EVI data
evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
evi.r <- crop(evi.r, bos.sizes)
evi.r <- mask(evi.r, bos.sizes)
bos.dat[,evi:=as.vector(evi.r)]
write.csv(bos.dat, "processed/boston.30m.agg.csv")


towns <- readOGR(dsn="data/towns/TOWNS_POLY.shp", layer="TOWNS_POLY")
towns <- spTransform(towns, crs(evi.r))
plot(towns[towns@data$TOWN=="PETERSHAM",])




####
### whole-AOI data aggregated to 30 m Landsat grid -- only have ISA AOI LULC and EVI as of Feb 22 2018
### read in AOI for crs description
# AOI <- readOGR(dsn="/projectnb/buultra/atrlica/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/", layer="AOI_simple_NAD83UTM19N")
# master.crs <- crs(AOI)

# load EVI composite for grid and process for ISA grid
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI.tif") ## this is the July 2010-2012 AOI EVI composite
# evi.r <- projectRaster(from=evi.r, crs=master.crs, res = 30, method = "ngb")
# writeRaster(evi.r, filename="processed/EVI/030005-6_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")

### aggregate raw 1m ISA to 30m
### don't fart around with the ISA grid until it's stacked with a 30 m EVI and cropped -- too hard to monkey with, leave in its native state until the end
### raster-native approach required, file too large for data.table -- takes a long time
### correct the 16 values to NA to get a proper aggregated mean value per cell (0 = not impervious, 1 = impervious)
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# isa.na.agg <- aggregate(isa.na, fact=30, fun=mean, na.rm=FALSE) # get mean ISA per 30 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.30m.tif", format="GTiff", overwrite=T)
# isa.na.agg <- raster("processed/isa.30m.tif")
# plot(isa.na.agg)

# ### align isa and evi grids, crop by AOI
# ### call python script for resampling 30m ISA to EVI grid
# pyth.path = './Rscripts/AOIISA_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)

isa.r <- raster("processed/isa30m_evigrd.tif")
# evi.r <- raster("processed/EVI/030005-6_2010-2012_EVI_NAD83.tif") ## this is the July 2010-2012 AOI EVI composite
# AOI <- readOGR(dsn="E:/BosAlbedo/data/AOI/AOI_simple_NAD83UTM19N/AOI_simple_NAD83UTM19N.shp", layer="AOI_simple_NAD83UTM19N")
# plot(isa.na.agg); plot(AOI, add=T)
# plot(evi.r); plot(AOI, add=T)
# AOI.r <- rasterize(AOI, evi.r)
# writeRaster(AOI.r, filename="processed/AOI.r.tif", format="GTiff", overwrite=T)
AOI.r <- raster("processed/AOI.r.tif")

### prep LULC --> EVI grid
pyth.path = './Rscripts/LULC_EVIgrid.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)
LULC.r <- raster("processed/LULC30m.tif")

LULC.r <- extend(LULC.r, isa.r)
evi.r <- crop(evi.r, isa.r)
AOI.r <- crop(AOI.r, isa.r)
evi.r <- mask(evi.r, mask = AOI.r)

evi.stack <- stack(evi.r, isa.r, LULC.r, AOI.r)
writeRaster(evi.stack, filename="processed/EVI30m.stack.tif", format="GTiff", overwrite=T)





### aggregate AOI-wide data to 250m MODIS grid
### similar process as for aggregating to Landsat 30m, only more differenter too.

## read in AOI for crs description
AOI <- readOGR(dsn="data/AOI/AOI_simple_NAD83UTM19N", layer="AOI_simple_NAD83UTM19N")
master.crs <- crs(AOI)

## load EVI composite for grid and process to produce master grid
# evi.r <- raster("processed/EVI/MOD_2010-2012_EVI.tif") ## this is the July 2010-2012 MODIS composite for the general AOI tile
# evi.r <- projectRaster(evi.r, crs=master.crs, res=250, method="bilinear")
# writeRaster(evi.r, filename="processed/EVI/MOD_2010-2012_EVI_NAD83.tif", format="GTiff", overwrite=T)
evi.r <- raster("processed/EVI/MOD_2010-2012_EVI_NAD83.tif") # note is not cropped
# AOI.r <- rasterize(AOI, evi.r)
# AOI.r <- crop(AOI.r, AOI, filename="processed/AOI.r.250m.tif", format="GTiff", overwrite=T) ## EVI 250m grid is master grid now
AOI.r <- raster("processed/AOI.r.250m.tif")

### aggregate raw 1m ISA to 250m (NOT to EVI grid yet!)
### raster-native approach required, file too large for data.table -- takes a long time, better to run on cluster to avoid moving the 1m file around
### 1) correct the 16 values to NA to get a proper aggregated mean value per cell (0 = not impervious, 1 = impervious)
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# fun <- function(x){
#   x[x==16] <- NA
#   return(x)
# }
# isa.na <- calc(isa, fun)
# isa.na.agg <- aggregate(isa.na, fact=250, fun=mean, na.rm=T) # get mean ISA per 250 m footprint
# writeRaster(isa.na.agg, filename="processed/isa.250m.tif", format="GTiff", overwrite=T)
# isa.250m <- raster("processed/isa.250m.tif")
# isa.250m <- projectRaster(from=isa.250m, to=evi.r, filename="processed/isa.250m.tif", format="GTiff", overwrite=T)
isa.250m <- raster("processed/isa.250m.tif")


### prep LULC data into collpased categories and get fractional area for each 250m gridcell sim. to 1m-->30m
lu.classnames <- c("forest", "dev", "hd.res", "med.res", "low.res", "lowveg", "other", "water")
lu.forest <- c(3,37) # Forest, FWet
lu.dev <- c(15,16,7,8,18,39,24,31,19,9,29) # Comm, Ind, PartRec, SpectRec, Transp., Junk, Util, PubInst, Waste, WBRec, Marina
lu.hdres <- c(11,10) # HDResid., MFResid.,
lu.medres <- c(12) # MDResid.
lu.lowres <- c(13, 38) #LDResid., VLDResid.
lu.lowveg <- c(4,1,2,6,5,26,14,25,34,23) # NFwet, PubInst, Crop, Past, Open, Mining, Golf, SWwet, SWBeach, Cem, CBog
lu.other <- c(17,35,40,36) #Tran, Orch, Brush, Nurs
lu.water <- c(20)
lulc.tot <- list(lu.forest, lu.dev, lu.hdres, lu.medres, lu.lowres, lu.lowveg, lu.other, lu.water)

## flexible blockwise function for flagging different class groups of LULC
lulc.flag <- function(x, lu.spot, filename) { # x is 30m lulc, lu.spot is list of LULC targeted
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
    o <- rep(0, length(v))
    o[v%in%lu.spot] <- 1
    out <- writeValues(out, o, bs$row[i])
    print(paste("finished block", i, "of", bs$n, "in", lu.classnames[l]))
  }
  out <- writeStop(out)
  return(out)
}

### loop through LULC class groups and create separately flagged 30m LULC rasters
aoi <- raster("processed/AOI.r.tif")
aoi <- crop(aoi, AOI, filename="processed/AOI.r.30m.tif", format="GTiff", overwrite=T)
LULC.30m <- raster("processed/LULC30m.tif")
for(l in 1:length(lulc.tot)){
  tmp <- do.call(lulc.flag,
          args = list(LULC.30m,
                      lulc.tot[[l]],
                      paste("processed/LULC.30m.", lu.classnames[l], "_only.tif", sep="")))
}

test <- raster("processed/LULC.30m.forest_only.tif")
plot(test); plot(AOI, add=T)

### can use the LULC 30m like Zhao did to describe fractional cover per MODIS cell.
### because 30m can't be easily aggregated to 250m, first aggregate here to 240m (x8), then resample/snap to EVI grid in Arc
aoi.30m <- raster("processed/AOI.r.30m.tif")
aoi.250m <- raster("processed/AOI.r.250m.tif")
evi.r <- raster("processed/EVI/MOD_2010-2012_EVI_NAD83.tif")
gimme <- list.files("processed/")
gimme <- gimme[grep(gimme, pattern="LULC.30m.")]
gimme <- gimme[!grepl(gimme, pattern=".ovr")]
gimme <- gimme[!grepl(gimme, pattern=".aux")]
fields <- sub("LULC.30m. *(.*?) *_only.tif.*", "\\1", gimme)
for(g in 1:length(gimme)){
  tmp <- raster(paste("processed/", gimme[g], sep=""))
  tmp <- extend(tmp, aoi.30m)
  tmp <- mask(tmp, aoi.30m)
  tmp <- crop(tmp, aoi.30m)
  tmp <- aggregate(tmp, fact=8, fun=mean, expand=T) ## aggregate to roughly the MODIS grid size by fraction of each LULC class
  resample(tmp, evi.r, method="bilinear",   ### get the mean fractional area resampled to the slightly bigger grid
                filename=paste("processed/LULC.250m.", fields[g], ".frac.tif", sep=""),
                format="GTiff",
                overwrite=T)
  print(paste("finished aggregating", fields[g]))
}
# test <- raster("processed/LULC.250m.dev.frac.tif")
# crs(test)
# plot(test)

## this python script doesn't try to resample the LULC 240m post-aggregate stuff to the MODIS grid because of ongoing errors in Arc, no idea how to fix
## resampling to MODIS grid is done above using resample(). Don't know why I didn't do this before.
### reproject/resample the 250m ISA to the MODIS EVI grid
# pyth.path = './Rscripts/AOI250m_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)

isa.r <- raster("processed/isa.250m.res.tif")
evi.r <- raster("processed/EVI/MOD_2010-2012_EVI_NAD83.tif") ## this is the July 2010-2012 AOI EVI composite
AOI.r <- raster("processed/AOI.r.250m.tif")
LULC.l <- list()
## load up the LULC fractions
for(f in 1:length(fields)){
  r <- raster(paste("processed/LULC.250m.", fields[f], ".frac.tif", sep=""))
  LULC.l[[f]] <- r
}
LULC.r <- stack(LULC.l)
LULC.r <- crop(LULC.r, AOI.r)
isa.r <- crop(isa.r, AOI.r)
isa.r <- mask(isa.r, AOI.r)
evi.r <- crop(evi.r, AOI.r)
evi.r <- mask(evi.r, AOI.r)
AOI.250m.stack <- stack(evi.r, isa.r, LULC.r, AOI.r)
names(AOI.250m.stack) <- c("EVI.250m", "ISA.250m", "dev.frac.250m", "forest.frac.250m", 
                          "hdres.frac.250m", "lowres.frac.250m", "lowveg.frac.250m", "medres.frac.250m",
                          "other.frac.250m", "water.frac.250m", "aoi.250m")
writeRaster(AOI.250m.stack, filename="processed/AOI.250m.stack.tif", format="GTiff", overwrite=T)
AOI.250m.dat <- as.data.table(as.data.frame(AOI.250m.stack))
write.csv(AOI.250m.dat, "processed/AOI.250m.dat.csv")

AOI.1km.stack <- aggregate(AOI.250m.stack, fact=4, fun=mean, na.rm=T, exapnd=T, filename="processed/AOI.1km.stack.tif")
AOI.1km.stack <- stack("processed/AOI.1km.stack.tif")
AOI.r <- rasterize(AOI, AOI.1km.stack[[1]])
# writeRaster(AOI.r, filename="processed/AOI.1km.r.tif", format="GTiff", overwrite=T)
AOI.r <- raster("processed/AOI.1km.r.tif")
AOI.1km.stack[[11]] <- AOI.r
AOI.1km.dat <- as.data.frame(AOI.1km.stack)
names(AOI.1km.dat)<- c("EVI.1km", "ISA.1km", "dev.frac.1km", "forest.frac.1km", 
                         "hdres.frac.1km", "lowres.frac.1km", "lowveg.frac.1km", "medres.frac.1km",
                         "other.frac.1km", "water.frac.1km", "aoi.1km")
write.csv(AOI.1km.dat, "processed/AOI.1km.dat.csv")
AOI.1km.stack <- stack("processed/AOI.1km.stack.tif")
plot(AOI.1km.stack)
# ### Test Code: Process NAIP 1m CIR data to NDVI
# ### NOTE SO MUCH: BAND ORDER IN NAIP CIR DATA IS RED GREEN BLUE NIR (1,2,3,4)
# naip.test <- stack("Z:/Ultra2/Users/atrlica/FragEVI/NAIP/EssexCIR/m_4207123_ne_19_h_20160706/m_4207123_ne_19_h_20160706.tif")
# naip.dat <- as.data.table(as.data.frame(naip.test))
# names(naip.dat) <- c("band1", "band2", "band3", "band4")
# naip.dat[, ndvi:=((band4-band1)/(band4+band1))]
# ndvi.r <- raster(naip.test[[1]])
# ndvi.r <- setValues(x = ndvi.r, values=naip.dat[,ndvi])
# writeRaster(ndvi.r, filename="processed/NAIP.ndvi.test.tif", format="GTiff", overwrite=T)
# 
# 
# #### Quickbird 4 band is still in DN form -- assume not corrected for atmo or orthorectified. Possible to get some of the VI's calculated with uncorrected data but quality probably low.
# #### Produce several VI's from the Quickbird 4 band imagery
# #### the below raw pan-sharpened data is 60cm, but is reported XXXXX integer -- assume this is corrected/orthorectified reflectance in integer form?
# qb <- stack("Z:/Ultra1/Data/Remote_Sensing/Quickbird/Mosaic_pansharp/Boston_pan.tif")
# names(qb) <- c("blue", "green", "red", "nir")
# ### on closer inspection, this is likely to be the DN's stored as 16 bit integer value (range 0-65535)
# ### NIR tops out the sensor reading but other bands do not
# 
# ## out of curiosity pull a raw file from the DVD and see what the number ranges are
# raw <- raster("Z:/Ultra1/Data/Remote_Sensing/Quickbird/DVD/DVD1/052187149010_01_P001_MUL/07AUG03160142-M2AS_R1C2-052187149010_01_P001.TIF")
# ## yes, same 16 bit integer
# 
# qb.vi <- function(x, y, z, fn.evi, fn.evi2, fn.ndvi, fn.nirv){ ### enter as red blue nir, filename evi, filename evi2, filename ndvi, filename nirv
#   out.evi <- raster(x)
#   out.evi2 <- raster(x)
#   out.ndvi <- raster(x)
#   out.nirv <- raster(x)
#   bs <- blockSize(out.evi)
#   out.evi <- writeStart(out.evi, fn.evi, overwrite=T, format="GTiff")
#   out.evi2 <- writeStart(out.evi2, fn.evi2, overwrite=T, format="GTiff")
#   out.ndvi <- writeStart(out.ndvi, fn.ndvi, overwrite=T, format="GTiff")
#   out.nirv <- writeStart(out.nirv, fn.nirv, overwrite=T, format="GTiff")
#   
#   for(i in 1:bs$n){
#     red <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## red band
#     red <- red/100000
#     blue <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## blue band
#     blue <- blue/100000
#     nir <- getValues(z, row=bs$row[i], nrows=bs$nrows[i]) ## nir band
#     nir <- nir/100000
#     ## kill values that are all 0 (=NA)
#     red.na <- red
#     red.na[red==0 & blue==0 & nir==0] <- NA
#     blue.na <- blue
#     blue.na[red==0 & blue==0 & nir==0] <- NA
#     nir.na <- nir
#     nir.na[red==0 & blue==0 & nir==0] <- NA
#     ## run VI calculations
#     evi.calc <- 2.5*((nir.na-red.na)/(1+nir.na+(6*red.na-7.5*blue.na)))
#     evi.calc[evi.calc>2.5 | evi.calc<(-2.5)] <- NA
#     evi2.calc <- 2.5*((nir.na-red.na)/(1+nir.na+(2.4*red.na)))
#     evi2.calc[evi2.calc<(-2.5) | evi2.calc>2.5] <- NA
#     ndvi.calc <- (nir.na-red.na)/(nir.na+red.na)
#     ndvi.calc[ndvi.calc>1 | ndvi.calc<(-1)] <- NA
#     nirv.calc <- ndvi.calc*nir.na
#     out.evi <- writeValues(out.evi, evi.calc, bs$row[i])
#     out.evi2 <- writeValues(out.evi2, evi2.calc, bs$row[i])
#     out.ndvi <- writeValues(out.ndvi, ndvi.calc, bs$row[i])
#     out.nirv <- writeValues(out.nirv, nirv.calc, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out.evi <- writeStop(out.evi)
#   out.evi2 <- writeStop(out.evi2)
#   out.ndvi <- writeStop(out.ndvi)
#   out.nirv <- writeStop(out.nirv)
#   return(out.evi)
#   return(out.evi2)
#   return(out.ndvi)
#   return(out.nirv)
# }
# s <- qb.vi(qb[["red"]], qb[["blue"]], qb[["nir"]], 
#             fn.evi="processed/boston/bos.evi.1m.tif", 
#             fn.evi2="processed/boston/bos.evi2.1m.tif",
#             fn.ndvi="processed/boston/bos.ndvi.1m.tif",
#             fn.nirv="processed/boston/bos.nirv.1m.tif")
# 
# ## this seems to give reasonable EVI values but seem pretty low even in dense canopy (~0.25 where NDVI values show 0.6)
# ## also there is a difference in registration (several meters) between this file and the other 1m boston rasters.
# evi.bos <- raster("processed/boston/bos.evi.1m.tif")
# evi2.bos <- raster("processed/boston/bos.evi2.1m.tif")
# ndvi.bos <- raster("processed/boston/bos.ndvi.1m.tif")
# nirv.bos <- raster("processed/boston/bos.nirv.1m.tif")

