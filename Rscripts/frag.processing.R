### process data layers to extract fragmentation per grid cell
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)


###
### Boston high-res data -- where are the high-biomass "forests"?
#####
## have used a rough threshold of 20k kg per pixel (equiv) as a cutoff for where it not longer
## makes sense to treat things as a mixed urban/tree pixel and more like a forest (about 100 MgC/ha)
((20000/2000)/900)*1E4 ### 111 MgC/ha

## would like to know what sort of pixels are over 20k kg
# bos.biom <- raster("processed/boston/bos.biom30m.tif")
# #### this function reclasses everything along 20k kg threshold
# big.biom <- function(biom, filename) {
#   out <- raster(biom)
#   bs <- blockSize(biom)
#   out <- writeStart(biom, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     target <- getValues(biom, row=bs$row[i], nrows=bs$nrows[i])
#     target[target<20000] <- 0
#     target[target>=20000] <- 1
#     out <- writeValues(out, target, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# big.biom(bos.biom, "processed/boston/bos.biom20k.tif")
# big <- raster("processed/boston/bos.biom20k.tif")
# plot(big)
#####

###
### create a 1m AOI raster
#####
### get Boston limits, raterize to 1m canopy and EVI 30m grid
bos.can <- raster("processed/boston/bos.can.redux.tif")
towns <- readOGR(dsn = "F:/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
bos.AOI <- towns[towns@data$TOWN=="BOSTON",]
bos.AOI <- bos.AOI[bos.AOI@data$SHAPE_AREA>1E07,] ## remove Harbor Islands
bos.AOI <- spTransform(bos.AOI, crs(bos.can))
bos.AOI@data$include <- 1
bos.AOI.r <- rasterize(bos.AOI[bos.AOI@data$include], bos.can)
aoi.fix <- function(x, filename) { # x is lulc, lu.defs is the list of collapsed classes, filename=output
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## AOI
    v[!is.na(v)] <- 1
    out <- writeValues(out, v, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
bos.AOI.r <- aoi.fix(bos.AOI.r, filename="processed/boston/tmp.tif")
bos.AOI.r <- crop(bos.AOI.r, bos.AOI)
writeRaster(bos.AOI.r, filename="processed/boston/bos.aoi.tif", format="GTiff", overwrite=T)
#####

###
### process 1m ISA 
#####
isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
towns <- readOGR(dsn = "/projectnb/buultra/atrlica/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
bos.bord <- towns[towns@data$TOWN=="BOSTON",]
bos.bord <- spTransform(bos.bord, crs(isa))
bos.isa <- crop(isa, bos.bord)
writeRaster(bos.isa, filename="processed/bos_isa_1m.tif", format="GTiff", overwrite=T)
# python.load("Rscripts/bosISA_resamp.py") ### test this
bos.isa <- raster("processed/isa_cangrid.tif")

### handling ISA registration difficulties
### 1m ISA shows misalignment to canopy/biomass features in places, particularly Allston
### first attempt was to manually reregister to the cov==barren layer 
### some improvement was made in reregistration, but not perfect in W part of Boston
### this layer lives as processed/boston/bos.isa.rereg.tif
### this layer was manually resampled to 30m landsat EVI grid as bos.isa.rereg30m.tif
### Second attempt (much better) was to manually reregister ISA according to road features on MassGIS road centerlines
### this version is bos.isa.RR2.tif --> the 30m aggregate updated bos.isa30m.tif
#####

###
### Create lumped LULC map
#####
### call python script for rasterize LULC in Boston to 1m canopy grid
# pyth.path = './Rscripts/LULC_bos_rast.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)

### get a clean copy of the LULC 1m raster
# bos.lulc <- raster("processed/boston/LU_bos_r1m.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# # bos.lulc <- crop(bos.lulc, bos.aoi)
# # bos.lulc <- mask(bos.lulc, bos.aoi)
# # writeRaster(bos.lulc, filename="processed/boston/bos.lulc_only.tif", format="GTiff", overwrite=T)
# bos.lulc <- raster("processed/boston/bos.lulc_only.tif")

## decide on a collapsed LULC scheme to flag for area fraction
lu.classnames <- c("forest", "dev", "hdres", "ldres", "lowveg", "water")
lu.forest <- c(3,37) # Forest, FWet
lu.dev <- c(5,8,15,16,17,18,19,29,31,36,39) #Mining, Spect-rec, Comm, Ind, Transitional, Transp, Waste Disp, Marina, Urb Pub/Inst., Nursery, Junkyard
lu.hdres <- c(10,11) # HDResid., MFResid.,
lu.ldres <- c(12,13,38) # MDResid., LDResid, VLDResid
lu.lowveg <- c(1,2,4,6,7,9,14,25,26,34,40) # Crop, pasture, open, part-rec, water-rec, SWwet, SWbeach, Golf, Cemetery, Brushland
lu.water <- c(20)
lulc.tot <- list(lu.forest, lu.dev, lu.hdres, lu.ldres, lu.lowveg, lu.water)

lulc.lump <- function(x, lu.defs, filename) { # x is lulc, lu.defs is the list of collapsed classes, filename=output
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ##  LULC raster raw LUCODES
    o <- v
    o[v%in%lu.defs[[1]]] <- 1 ## forest
    o[v%in%lu.defs[[2]]] <- 2 ## dev
    o[v%in%lu.defs[[3]]] <- 3 ## hdres
    o[v%in%lu.defs[[4]]] <- 4 ## ldres
    o[v%in%lu.defs[[5]]] <- 5 ## lowveg
    o[v%in%lu.defs[[6]]] <- 6 ## water
    out <- writeValues(out, o, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
## 1m lumped LULC raster
bos.lulc <- raster("processed/boston/bos.lulc_only.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
lulc.lump(bos.lulc, lulc.tot, "processed/boston/bos.lulc.lumped.tif")

bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")
plot(bos.lulc)

## create comparable single 30m raster with collapsed LULC classes
## quick and dirty arc process on polygons poly-->raster at 30m EVI grid
bos.lulc30m <- raster("processed/boston/bos_lulc30m.tif") ## this is full 40 class LULC scheme
bos.aoi <- raster("processed/boston/bos.aoi30m.tif")
bos.lulc30m <- extend(bos.lulc30m, bos.aoi)
bos.lulc30m <- crop(bos.lulc30m, bos.aoi)
bos.lulc30m <- mask(bos.lulc30m, bos.aoi)
lulc.lump(bos.lulc30m, lulc.tot, "processed/boston/bos.lulc30m.lumped.tif")

bos.lulc30m <- raster("processed/boston/bos.lulc30m.lumped.tif")
plot(bos.lulc30m)

## flexible blockwise function for flagging different class groups of LULC
lulc.flag <- function(x, lu.spot, aoi, filename) { # x is 1m lulc, lu.spot is list of LULC targeted, aoi is aoi 1m raster
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## 1m LULC raster
    a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) # AOI mask
    o <- rep(0, length(v))
    o[v%in%lu.spot] <- 1 ## flag the in-category with 1
    o[is.na(a)] <- NA ## cancel values outside of AOI
    out <- writeValues(out, o, bs$row[i])
    print(paste("finished block", i, "of", bs$n, "in", lu.classnames[l]))
  }
  out <- writeStop(out)
  return(out)
}

### loop through LULC class groups and create separately flagged LULC rasters
for(l in 1:length(lulc.tot)){
  tmp <- do.call(lulc.flag,
                 args = list(bos.lulc,
                             lulc.tot[[l]],
                             bos.aoi,
                             paste("processed/boston/bos.", lu.classnames[l], "_only.tif", sep="")))
  # print(paste("masking to Boston AOI"))
  # tmp <- mask(tmp, bos.aoi)
  # writeRaster(tmp, filename=paste("processed/boston/bos.", lu.classnames[l], "_only.tif", sep=""), format="GTiff", overwrite=T)
}
# bos.lulc <- raster("processed/boston/bos.lulc_only.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# lulc.flag(bos.lulc, lulc.tot[[1]], bos.aoi, "processed/boston/bos.forest.tif")

bos.forest <- raster("processed/boston/bos.forest_only.tif")
bos.dev <- raster("processed/boston/bos.dev_only.tif")
bos.hdres <- raster("processed/boston/bos.hdres_only.tif")
bos.ldres <- raster("processed/boston/bos.ldres_only.tif")
bos.lowveg <- raster("processed/boston/bos.lowveg_only.tif")
bos.water <- raster("processed/boston/bos.water_only.tif")
plot(bos.forest, main="forest")
plot(bos.dev, main="Developed")
plot(bos.hdres, main="HDres")
plot(bos.ldres, main="LDRes")
plot(bos.lowveg, main="LowVeg")
plot(bos.water, main="water")
#####

###
### Create correct 1m canopy map from biomass map
#####
can.fix <- function(can, aoi, filename) { # x is biomass, a is AOI, filename=output
  out <- raster(can)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(can, row=bs$row[i], nrows=bs$nrows[i]) ##  biomass
    a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ##  biomass
    v[v>0] <- 1
    v[is.na(a)] <- NA
    out <- writeValues(out, v, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
bos.can <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.can <- crop(bos.can, bos.aoi)
s <- can.fix(bos.can, bos.aoi, filename="processed/boston/bos.can.redux.tif")
plot(s)

### fuck it do the biomass file while you'r here
bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.biom <- crop(bos.biom, bos.aoi)
writeRaster(bos.biom, filename="processed/boston/bos.biom.tif", format="GTiff", overwrite=T)
#####

### Boston high-res data -- vegetation character map
### combine 2.4m NDVI and canopy map to produce 1m veg classification map
### use Arcmap to resample + snap to 1m grid (bilinear) for NDVI 2.4m -- get good alignment with original data and features in canopy map
#####

## Jan 2019: Raciti's reported canopy coverage was based on a map that translates biomass>0 to can==1
## However, dataverse "canopy" layer has been altered via unknown smoothing process, apparently exceeds the biomass>0 coverage and creates discrepancy with Raciti et al. 2014 report
## Step 0: Create 1m canopy presence/absence from 1m biomass map
biom1m <- raster("data/dataverse_files/bostonbiomass_1m.tif")
can.from.biom <- function(biom, filename) {
  out <- raster(biom)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    r <- getValues(biom, row=bs$row[i], nrows=bs$nrows[i])
    r[r>0] <- 1
    out <- writeValues(out, r, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
cl <- can.from.biom(biom1m, "processed/boston/bos.can.redux.tif")
plot(cl)

### (must be performed on desktop): put NDVI.img through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
# pyth.path = './Rscripts/NDVI_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(paste("ArcPy working on NDVI resample: ", output))


### V3 cover mapping: Uses re-registered ISA layer
# bos.isa <- raster("processed/boston/isa.reg-res2.tif")
bos.ndvi <- raster("processed/boston/bos.ndvi.tif")
bos.can <- raster("processed/boston/bos.can.tif")
bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.ed <- raster("processed/boston/bos.ed10.tif")
# bos.isa <- crop(bos.isa, bos.aoi)
# writeRaster(bos.isa, "processed/boston/bos.isa.RR2.tif", format="GTiff", overwrite=T)
bos.isa <- raster("processed/boston/bos.isa.RR2.tif")


#### THIS VERSION spit out a 6 class raster, splits canopy by edge position
veg.class <- function(can, isa, ed, ndvi, lulc, aoi, filename) {
  out <- raster(can)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    target <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
    check <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) 
    target[check==1] <- 999 ### to flag any pixels that escape classification
    w <- getValues(can, row=bs$row[i], nrows=bs$nrows[i]) ## canopy
    x <- getValues(isa, row=bs$row[i], nrows=bs$nrows[i]) ## isa
    y <- getValues(ndvi, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi
    z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
    e <- getValues(ed, row=bs$row[i], nrows=bs$nrows[i]) ## edge
    target[w==1 & e==1] <- 6 # canopy edge
    target[w==1 & e==0] <- 5 # canopy interior
    target[w==0 & x==1] <- 4 # non-veg impervious
    target[w==0 & x==0 & y>=0.25] <- 2 # grass
    target[w==0 & x==0 & y<0.25 & z!=6] <- 3 # non-veg pervious
    target[w==0 & x==0 & y<0.25 & z==6] <- 1 # water
    target[is.na(check) | is.na(w) | is.na(x) | is.na(y) | is.na(z)] <- NA ## cancel missing values
    out <- writeValues(out, target, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
cl <- veg.class(bos.can, bos.isa, bos.ed, bos.ndvi, bos.lulc, bos.aoi, "processed/boston/bos.cov.V3-ed.tif")
plot(cl)


### THIS VERSION spits out 6 class raster, splits canopy according to position over impervious vs. pervious
### which allows for a nice see-through effect if colors are wisely chosen
# bos.isa <- raster("processed/boston/isa.reg-res2.tif")
bos.isa <- raster("processed/boston/bos.isa.RR2.tif")
bos.ndvi <- raster("processed/boston/bos.ndvi.tif")
bos.can <- raster("processed/boston/bos.can.redux.tif")
bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
crs(bos.isa); extent(bos.isa)
crs(bos.ndvi); extent(bos.ndvi)
crs(bos.can); extent(bos.can)
crs(bos.lulc); extent(bos.lulc)
crs(bos.aoi); extent(bos.aoi)
bos.can <- crop(bos.can, bos.aoi)


#### this function spits out a single 6 class raster
veg.class <- function(can, isa, ndvi, lulc, aoi, filename) {
  out <- raster(can)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    target <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
    check <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) 
    target[check==1] <- 999 ### to flag any pixels that escape classification
    w <- getValues(can, row=bs$row[i], nrows=bs$nrows[i]) ## canopy
    x <- getValues(isa, row=bs$row[i], nrows=bs$nrows[i]) ## isa
    y <- getValues(ndvi, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi
    z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
    target[w==1 & x==0] <- 6 # canopy over pervious
    target[w==1 & x==1] <- 5 # canopy over impervious
    target[w==0 & x==1] <- 4 # non-veg impervious
    target[w==0 & x==0 & y>=0.25] <- 2 # grass
    target[w==0 & x==0 & y<0.25 & z!=6] <- 3 # non-veg pervious
    target[w==0 & x==0 & y<0.25 & z==6] <- 1 # water
    target[is.na(check) | is.na(w) | is.na(x) | is.na(y) | is.na(z)] <- NA ## cancel missing values or anything out of aoi
    out <- writeValues(out, target, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
cl <- veg.class(bos.can, bos.isa, bos.ndvi, bos.lulc, bos.aoi, "processed/boston/bos.cov.V4-canisa.tif")
plot(cl)

library(viridis)
image(cl, col=magma(6))
image(cl, col=plasma(6))
image(cl, col=viridis(6))
image(cl, col=cividis(6))
cov.pal <- magma(6)
cov.pal <- plasma(6)
cov.pal <- cividis(6)
cov.pal <- c(cov.pal[1],cov.pal[6],cov.pal[3],cov.pal[2],cov.pal[4],cov.pal[5]) ## this seems like an intuitive structuring
image(cl, col=cov.pal)
col2rgb(cov.pal) ## how to enter these assholes into arc RGB
#####

##
#### This is part of process to ID plantable road buffer space
#####
### Now use cover map to locate near-road pervious pixels without canopy cover
rdbuff <- raster("E:/FragEVI/processed/boston/ROAD10m_rast.tif")
cov <- raster("E:/FragEVI/processed/boston/bos.cov.V3-canisa.tif")
rdbuff <- crop(extend(rdbuff, cov), cov)
bos.aoi <- raster("processed/boston/bos.aoi.tif")

#### find cover class 3 or 2 (non-veg pervious, veg pervious) in the buffer region
road.planting <- function(roadbuffer, cover, aoi, filename) {
  out <- raster(aoi)
  bs <- blockSize(aoi)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    target <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## dummy chunk
    road <- getValues(roadbuffer, row=bs$row[i], nrows=bs$nrows[i]) ## road buffer chunk, 1=buffer/NA
    plantme <- getValues(cover, row=bs$row[i], nrows=bs$nrows[i]) ## cover chunk
    here <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi chunk
    target[here==1] <- 0 ## blank values for anything inside AOI
    target[road==1 & plantme %in% c(2,3) & here==1] <- 1 ## mark the places inside AOI that are suitable
    out <- writeValues(out, target, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
cj <- road.planting(rdbuff, cov, bos.aoi, "processed/boston/bos.planting10m.tif")
plot(cj)

## converted to polygon in arc
### Arc is being an asshole, won't filter the polygons properly, so do it here
library(rgdal)
library(raster)
## this is the unsimplified polygon, follows raster edges perfectly
# plantP <- readOGR("processed/boston/plant10m_poly.shp")
# plantP.filt <- plantP[plantP@data$gridcode=="1000000" & plantP@data$Shape_Area>=2,]
# # writeOGR(plantP.filt, "processed/boston/plant10m_poly_filt.shp", "plant10m_poly_filt", driver="ESRI Shapefile")
plantP <- readOGR("processed/boston/plant10mPoly_Simp.shp")
plantP.filt <- plantP[plantP@data$gridcode=="1000000" & plantP@data$Shape_Area>1.3,]
writeOGR(plantP.filt, "processed/boston/plant10m_poly_simp_filt.shp", "plant10m_poly_simp_filt", driver="ESRI Shapefile")
### put this shit into arc and run a buffer on the filtered (simplified) polygons

#### read in the 4m buffered filtered polygons and filter for >50% open sky
# plantP.filt <- readOGR("processed/boston/plant10m_poly_filt.shp") ## 78618 features
# plantP.filt.buff <- readOGR("processed/boston/plant10m_poly_filt_4mBuff.shp") ##78618 features
plantP.filt.buff <- readOGR("processed/boston/plant10m_poly_simp_filt_4mBuff.shp") ##82527 features

### clean up data frame, filter for canopy sky view factor
# sum(as.numeric(as.character(plantP.filt@data$OBJECTID))-as.numeric(as.character(plantP.filt@data$Id))) ## identical
filt.d <- plantP.filt@data[,c("OBJECTID", "Shape_Leng", "Shape_Area")]
names(filt.d)[2:3] <- c("plant_Leng", "plant_Area")
buff.d <- plantP.filt.buff@data[,c("OBJECTID", "gridcode", "BUFF_DIST", "Shape_Le_1", "Shape_Area")]
names(buff.d)[4:5] <- c("buff_Leng", "buff_Area")
merg.d <- merge(x=filt.d, y=buff.d, by="OBJECTID")
plantP.filt.buff@data$OBJECTID_1 <- NULL
plantP.filt.buff@data$Id <- NULL
plantP.filt.buff@data$BUFF_DIST <- NULL
plantP.filt.buff@data$Shape_Area <- NULL
plantP.filt.buff@data$Shape_Le_1 <- NULL
plantP.filt.buff@data$Shape_Leng <- NULL
plantP.filt.buff@data$ORIG_FID <- NULL
plantP.filt.buff@data$gridcode <- NULL
plantP.filt.buff@data <- merge(x=plantP.filt.buff@data, y=merg.d, by="OBJECTID")

### filter file by extracting canopy area within the buffer, then eliminating polygons associated with buffers that have too much canopy already
keep <- integer()
ditch <- integer()
bos.can <- raster("processed/boston/bos.can.tif")
for(f in 1:dim(plantP.filt.buff@data)[1]){
  a=Sys.time()
  test <- extract(bos.can, plantP.filt.buff[f,])
  if(sum(test[[1]], na.rm=T)<0.5*plantP.filt.buff@data[f,"buff_Area"]){ ## if less than half the buffer area already in canopy
    keep <- c(keep, f)
  } else{
    ditch <- c(ditch, f)
  }
  print(paste("evaluated polygon", f, Sys.time()-a))
}

plantP.canfilt <- plantP.filt.buff[keep,]
plantP.canReject <- plantP.filt.buff[ditch,]

### this hangs forever it seems
# test <- extract(bos.can, plantP.filt.buff[,])
# # plot(crop(bos.can, extent(plantP.filt.buff[19,])))
# # plot(plantP.filt.buff[19,], add=T)
# # ((5.394/10)*dim(plantP.filt@data)[1])/60/60 ## this is an 12 hour extract
# sum.na <- function(x){sum(x, na.rm=T)}
# plantP.filt.buff@data$can_area <- sapply(test, sum.na)
# plantP.canfilt <- plantP.filt.buff[can_area<(0.5*buff_Area),]
# plantP.canReject <- plantP.filt.buff[can_area>=(0.5*buff_Area),]

writeOGR(plantP.canfilt, "processed/boston/plant10m_simp_4mbuff_canOK.shp", "plant10m_simp_4mbuff_canOK", driver="ESRI Shapefile")
writeOGR(plantP.canReject, "processed/boston/plant10m_simp_4mbuff_canBAD.shp", "plant10m_simp_4mbuff_canBAD", driver="ESRI Shapefile")
#####

###
### Processing 1m canopy map for edge class
#####
# ### call python scripbt for identifying canopy edge distance (desktop only)
# pyth.path = './Rscripts/canopy_process.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)
### notes: Feb 6 2019 -- have reprocessed this with the un-smoothed canopy cover 1m raster (had to simplify polygons)

### Canopy area by cumulative distance from edge
can.sum <- function(c, b, l) { # c is full canopy, b is buffer, l is lulc
  bs <- blockSize(c)
  tot <- integer()
  forest <- integer()
  dev <- integer()
  hdres <- integer()
  ldres <- integer()
  lowveg <- integer()
  water <- integer()
  for (i in 1:bs$n) {
    v <- getValues(c, row=bs$row[i], nrows=bs$nrows[i])
    buff <- getValues(b, row=bs$row[i], nrows=bs$nrows[i])
    lulc <- getValues(l, row=bs$row[i], nrows=bs$nrows[i])
    v[v==1 & buff==1] <- 0 ## anything canopy that gets labeled as interior to buff, 0 out
    tot <- c(tot, sum(v, na.rm=T))
    forest <- c(forest, sum(v[lulc==1], na.rm=T))
    dev <- c(dev, sum(v[lulc==2], na.rm=T))
    hdres <- c(hdres, sum(v[lulc==3], na.rm=T))
    ldres <- c(ldres, sum(v[lulc==4], na.rm=T)) 
    lowveg <- c(lowveg, sum(v[lulc==5], na.rm=T))
    water <- c(water, sum(v[lulc==6], na.rm=T))
    print(paste("finished block", i, "of", bs$n))
  }
  return(c(sum(tot, na.rm=T), sum(forest, na.rm=T), sum(dev, na.rm=T), 
           sum(hdres, na.rm=T), sum(ldres, na.rm=T), sum(lowveg, na.rm=T),
           sum(water, na.rm=T)))
}
# test <- raster("processed/boston/bos.nocan_10mbuff.tif")
# test <- (extend(test, bos.aoi))
# test <- crop(test, bos.aoi)
# d <- can.sum(bos.can, test, bos.lulc) ## this is exactly = all can area (by LULC) within 10m of edge (see Table S5)

### load up the list of canopy buffer files to process
can.buffs <- list.files("processed/boston/")
can.buffs <- can.buffs[grep(pattern = "bos.nocan_", x = can.buffs)]
can.buffs <- can.buffs[grep(pattern=".tif", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".vat", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".aux", x=can.buffs)]
can.buffs <- can.buffs[!grepl(pattern = ".ovr", x=can.buffs)]
buff.dist <- as.integer(unlist(rm_between(can.buffs, "nocan_", "mbuff", extract=TRUE)))

bos.can <- raster("processed/boston/bos.can.redux.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.can <- crop(bos.can, bos.aoi)
bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")
bos.lulc <- crop(bos.lulc, bos.aoi)

### set up the row 0 of the buffers area file -- row "0" is the whole area...
dist <- integer() ## keep track of the buffer distances as you process the rasters in whatever random order they come in
buff.area.inclusive <- data.frame()

## loop through successive buffer files and and get inclusive canopy buffer areas by buffer distance
for(g in 1:length(can.buffs)){
  print(paste("initializing buffer raster", can.buffs[g]))
  r <- raster(paste("processed/boston/", can.buffs[g], sep=""))
  ## make sure every goddamn one is the same ncell as aoi/can/lulc
  r <- extend(r, bos.aoi)
  r <- crop(r, bos.aoi)
  buff.area.inclusive <- rbind(buff.area.inclusive, can.sum(bos.can, r, bos.lulc)/1E4) ## calc total and LULC canopy area inside this buffer distance
  dist <- c(dist, buff.dist[g])
}
buff.area.inclusive <- cbind(dist, buff.area.inclusive)
names(buff.area.inclusive) <- c("buff.dist", "total.can.ha", "forest.can.ha", "dev.can.ha", "hdres.can.ha", "ldres.can.ha",
                        "lowveg.can.ha", "water.can.ha")


## get the single-distance ring areas for each buffer 
can.tot <- integer()
can.forest <- integer()
can.dev <- integer()
can.hdres <- integer()
can.ldres <- integer()
can.lowveg <- integer()
can.water <- integer()
buff.track <- integer()
for(d in 2:100){ ## working from buffer distance 2 up (will add in 1m buffer and whole canopy areas below)
  can.tot <- c(can.tot,
               buff.area.inclusive$total.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$total.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.forest <- c(can.forest,
                  buff.area.inclusive$forest.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$forest.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.dev <- c(can.dev,
               buff.area.inclusive$dev.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$dev.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.hdres <- c(can.hdres, 
                 buff.area.inclusive$hdres.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$hdres.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.ldres <- c(can.ldres,
                 buff.area.inclusive$ldres.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$ldres.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.lowveg <- c(can.lowveg,
                  buff.area.inclusive$lowveg.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$lowveg.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  can.water <- c(can.water,
                 buff.area.inclusive$water.can.ha[buff.area.inclusive$buff.dist==d]-buff.area.inclusive$water.can.ha[buff.area.inclusive$buff.dist==(d-1)])
  buff.track <- c(buff.track, d) ## we are calculating areas up to buffer distance 2, which requires buffer distance 1 to get isolated
}
## add in buffer distance 1 area
d=1
can.tot <- c(can.tot, buff.area.inclusive$total.can.ha[buff.area.inclusive$buff.dist==d])
can.forest <- c(can.forest, buff.area.inclusive$forest.can.ha[buff.area.inclusive$buff.dist==d])
can.dev <- c(can.dev, buff.area.inclusive$dev.can.ha[buff.area.inclusive$buff.dist==d])
can.hdres <- c(can.hdres, buff.area.inclusive$hdres.can.ha[buff.area.inclusive$buff.dist==d])
can.ldres <- c(can.ldres, buff.area.inclusive$ldres.can.ha[buff.area.inclusive$buff.dist==d])
can.lowveg <- c(can.lowveg, buff.area.inclusive$lowveg.can.ha[buff.area.inclusive$buff.dist==d])
can.water <- c(can.water, buff.area.inclusive$water.can.ha[buff.area.inclusive$buff.dist==d])
buff.track <- c(buff.track, d)

##package up processed ring areas
report <- data.frame(
  buff.dist=buff.track,
  total.can.ha=can.tot,
  forest.can.ha=can.forest,
  dev.can.ha=can.dev,
  hdres.can.ha=can.hdres,
  ldres.can.ha=can.ldres,
  lowveg.can.ha=can.lowveg,
  water.can.ha=can.water
)
### add in total canopy area by LULC
all.can <- read.csv("processed/boston/results/TableS5_area_totals.csv")
report <- rbind(report,
                     c(0, 
                       all.can$CAN[7],
                       all.can$CAN[1:6]))
## export report
report <- report[order(report$buff.dist),]
write.csv(report, "processed/boston/results/Canopy_area_by_edge_buffer.csv")

# ## have a look
# plot(report$buff.dist[2:101], report$total.can.ha[2:101], xlim=c(0, 50))
# points(report$buff.dist[2:101], report$forest.can.ha[2:101], xlim=c(0, 50), pch=15, cex=0.6, col="darkgreen")
# points(report$buff.dist[2:101], report$hdres.can.ha[2:101], xlim=c(0, 50), pch=15, cex=0.6, col="salmon")
# points(report$buff.dist[2:101], report$lowveg.can.ha[2:101], xlim=c(0, 50), pch=15, cex=0.6, col="goldenrod")
# points(report$buff.dist[2:101], report$dev.can.ha[2:101], xlim=c(0, 50), pch=15, cex=0.6, col="gray40")
### the distance 0 figures are damn close to the LULC totals in Table S5

###
### Create canopy edge class ring rasters (>10, 10-20, 20-30, >30)
### note: With simplified polygon step during the edge buffering process, the resulting
### canopy edge rasters are a half-pixel off the canopy grid; had to be manually resampled in Arc
ed1 <- raster("processed/boston/nocan_10mbuff_res.tif")
ed2 <- raster("processed/boston/nocan_20mbuff_res.tif")
ed3 <- raster("processed/boston/nocan_30mbuff_res.tif")
bos.can <- raster("processed/boston/bos.can.redux.tif")
bos.aoi <- raster("processed/boston/bos.aoi.tif")
## slice and dice every-fucking-thing
bos.can <- crop(bos.can, bos.aoi)
ed1 <- extend(ed1, bos.aoi) # edge and cover on same grid already
ed2 <- extend(ed2, bos.aoi)
ed3 <- extend(ed3, bos.aoi)
ed1 <- crop(ed1, bos.aoi)
ed2 <- crop(ed2, bos.aoi)
ed3 <- crop(ed3, bos.aoi)
# edges <- stack(ed1, ed2, ed3, bos.can) # just to make damn sure everything lines up

edges.bl <- function(x, y, z, filename) { # x is edge class, y is canopy flag, z is aoi
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge flag at designated distance
    g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## canopy cover
    a <- getValues(z, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
    t <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## dummy layer to modify
    v[!(v%in%c(0,1))] <- NA # kill any weird values that aren't coming from the nocan==0 buffer map
    t[g==1 & v==0] <- 1 ### wherever there's canopy that is inside the edge buffer (edge==0), flag it
    t[g==1 & v==1] <- 0 ### wherever there's canopy that did not get edge buffer flagged (ie edge==1), cancel it
    t[a!=1] <- NA ### cancel anything outside AOI
    out <- writeValues(out, t, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
s <- edges.bl(ed1, bos.can, bos.aoi, filename="processed/boston/bos.ed10m.redux.tif")
t <- edges.bl(ed2, bos.can, bos.aoi, filename="processed/boston/bos.ed20m.redux.tif")
u <- edges.bl(ed3, bos.can, bos.aoi, filename="processed/boston/bos.ed30m.redux.tif")

## make a corresponding biomass map that is just the <10m edge biomass
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.ed <- raster("processed/boston/bos.ed10m.redux.tif")
bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.biom <- crop(bos.biom, bos.aoi)
bos.ed <- crop(bos.ed, bos.aoi)
edge.biom <- function(x, y, z, filename) { # x is edge10 canopy, y is biomass, z is aoi
  out <- raster(y)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge canopy
    g <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## biomass
    a <- getValues(z, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
    g[v!=1] <- 0
    g[a!=1] <- NA
    out <- writeValues(out, g, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
ff <- edge.biom(bos.ed, bos.biom, bos.aoi, filename="processed/boston/bos.ed10m.biom.redux.tif")
#####

###
### Get area/biomass metrics for landscape structure summary
#####
bos.aoi <- raster("processed/boston/bos.aoi.tif")
bos.can <- raster("processed/boston/bos.can.redux.tif")
bos.ed <- raster("processed/boston/bos.ed10m.redux.tif")
bos.isa <- raster("processed/boston/bos.isa.RR2.tif")
bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")

# bos.aoi.tot <- sum(getValues(bos.aoi), na.rm=T) ## 12455 ha, lower than below
bos.can <- crop(bos.can, bos.aoi)
# bos.can.tot <- sum(getValues(bos.can), na.rm=T) ### 3144 ha, this is not the same as the number that comes out of the bottom here
bos.biom <- crop(bos.biom, bos.aoi)
# bos.biom.tot <- sum(getValues(bos.biom), na.rm=T) ## 357 GgC
bos.ed <- crop(bos.ed, bos.aoi)
# bos.ed.tot <- sum(getValues(bos.ed), na.rm=T) ## 2663 ha
bos.isa <- crop(bos.isa, bos.aoi)
# bos.isa.tot <- sum(getValues(bos.isa), na.rm=T) ## 7138 ha
bos.lulc <- crop(bos.lulc, bos.aoi)

sum.na <- function(x){sum(x, na.rm=T)}

bs <- blockSize(bos.aoi)
tmp.aoi <- matrix(ncol=6)
tmp.can <- matrix(ncol=6)
tmp.edge <- matrix(ncol=6)
tmp.isa <- matrix(ncol=6)
tmp.biom <- matrix(ncol=6)
forest <- data.table()
dev <- data.table()
hdres <- data.table()
ldres <- data.table()
lowveg <- data.table()
water <- data.table()
aoi.tab <- data.table()
for (i in 1:bs$n) {
  df <- data.frame(
    aoi=getValues(bos.aoi, row=bs$row[i], nrows=bs$nrows[i]),
    can=getValues(bos.can, row=bs$row[i], nrows=bs$nrows[i]),
    edge=getValues(bos.ed, row=bs$row[i], nrows=bs$nrows[i]),
    isa=getValues(bos.isa, row=bs$row[i], nrows=bs$nrows[i]),
    biom=getValues(bos.biom, row=bs$row[i], nrows=bs$nrows[i]),
    lulc=getValues(bos.lulc, row=bs$row[i], nrows=bs$nrows[i]))
  df <- as.data.table(df)
  forest <- rbind(forest, df[lulc==1, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  dev <- rbind(dev, df[lulc==2, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  hdres <- rbind(hdres, df[lulc==3, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  ldres <- rbind(ldres, df[lulc==4, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  lowveg <- rbind(lowveg, df[lulc==5, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  water <- rbind(water, df[lulc==6, .(sum.na(aoi),
                          sum.na(can),
                          sum.na(edge),
                          sum.na(isa),
                          sum.na(biom))])
  aoi.tab <- rbind(aoi.tab, df[, .(sum.na(aoi),
                                    sum.na(can),
                                    sum.na(edge),
                                    sum.na(isa),
                                    sum.na(biom))])
  
  print(paste("finished block", i, "of", bs$n))
}
area.tots <- rbind(apply(forest, MARGIN=2, FUN=sum.na),
        apply(dev, MARGIN=2, FUN=sum.na),
        apply(hdres, MARGIN=2, FUN=sum.na),
        apply(ldres, MARGIN=2, FUN=sum.na),
        apply(lowveg, MARGIN=2, FUN=sum.na),
        apply(water, MARGIN=2, FUN=sum.na),
        apply(aoi.tab, MARGIN=2, FUN=sum.na))
area.tots[,1:4] <- area.tots[,1:4]/1E4
area.tots[,5] <- area.tots[,5]/(2*1E6)
area.tots <- data.frame(area.tots)
area.tots <- cbind(c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total"),
                   area.tots)
names(area.tots) <- c("LULC", "AOI", "CAN", "EDGE", "ISA", "BIOM")
write.csv(cbind(c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total"),
                area.tots),
          "processed/boston/results/TableS5_area_totals.csv")
#####


###
### Compare canopy map from 2006 to map from 2014
#####
bos.can06 <- raster("data/dataverse_files/bostoncanopy_1m.tif")
bos.can14 <- raster("data/BosTrees2014/LC_boston2014_NAD83_1m.tif")
bos.aoi <- raster("processed/boston/bos.AOI.1m.tif")
bos.can06 <- crop(bos.can06, bos.aoi)
bos.can06 <- mask(bos.can06, bos.aoi)
bos.can14 <- crop(bos.can14, bos.aoi)
bos.can14 <- mask(bos.can14, bos.aoi)

## function to ID changes in canopy from 06-14
can.comp <- function(can1, can2, aoi, filename) { # can1 is can at t=1, can2 is can at t=2
  out <- raster(can1)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    a <- getValues(can1, row=bs$row[i], nrows=bs$nrows[i]) ## can1
    b <- getValues(can2, row=bs$row[i], nrows=bs$nrows[i]) # can2
    m <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) # aoi
    
    ## clean up data
    b[b==0] <- NA
    a[!a%in%c(0,1)] <- NA ## kill out of range values
    b[!b%in%c(1,2,3,4,5,6,7,8,9,10,11)] <- NA
    m[m!=1] <- NA
    
    ## make new vector that records differences
    A <- rep(0, length(a))
    # a[a==0 & b!=1] <- 0
    # a[a==1 & b==1] <- 0 # stable
    A[a==1 & b!=1] <- 2 # loss
    A[a==0 & b==1] <- 1 # gain
    A[is.na(m)] <- NA # mask the rest
    out <- writeValues(out, A, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
can.comp(bos.can06, bos.can14, bos.aoi, "processed/boston/bos.canChange0614.tif")

## plots for comparision
par(mar=c(1,1,4,1))
plot(bos.can06, main="Canopy 2006")
plot(bos.can14, main="Canopy 2014")
COMP <- raster("processed/boston/bos.canChange0614.tif")
plot(COMP, main="1=gain, 2=loss")


addme <- function(can1, value) { # can1 is can target, value is targeted value(s)
  out <- raster(can1)
  bs <- blockSize(out)
  g=0
  for (i in 1:bs$n) {
    a <- getValues(can1, row=bs$row[i], nrows=bs$nrows[i]) ## can1
    g=g+sum(a%in%value, na.rm=T)
    print(paste("finished block", i, "of", bs$n))
  }
  return(g)
}
aoi.tot <- addme(bos.aoi, 1) ## total city area (06 study boundary)
can06.tot <- addme(bos.can06, 1) ## total canopy area '06, 3941.64 ha
aoi06.tot <- addme(bos.can06, c(0,1)) ##12455.09 ha
can14.tot <- addme(bos.can14, 1) ## total canopy area '14, 3418.30 ha
aoi14.tot <- addme(bos.can14, c(1,2,3,4,5,6,7,8,9,10,11)) #12228.03 ha ## tiny registration (?) error means some edge trimming with aoi mask
gain <- addme(COMP, 1) ## area gained by 2014, 617.38 ha
loss <- addme(COMP, 2) ## area lost by 2014, 1125.14 ha
(gain-loss)/1E4 ## 507.75 ha loss
(can14.tot-can06.tot)/1E4 ### not quite the same -523.33 ha
can06.tot/aoi06.tot # 31.6% canopy in 06
can14.tot/aoi14.tot # 30.0% canopy in 14
sum(getValues(bos.can06), na.rm=T)/aoi.tot ## same 31.6% canopy
#####


###
### Aggregate 1m into 30m raster w Landsat grid
#####
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# bos.isa <- raster("processed/boston/bos.isa.RR2.tif")
# bos.can <- raster("processed/boston/bos.can.redux.tif")
# # bos.lulc <- raster("processed/boston/bos.lulc.lumped.tif")
# bos.ed10 <- raster("processed/boston/bos.ed10m.redux.tif")
# bos.biom <- raster("processed/boston/bos.biom.tif")

### get fractional area of each cover class per 30m aggregate cell
pyth.path = './Rscripts/bos1m_agg.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)

### prep LULC --> EVI grid
pyth.path = './Rscripts/LULC_EVIgrid.py'
output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE); print(output)
LULC.r <- raster("processed/LULC30m.tif")
#####


###
### aggregate AOI-wide data to 250m MODIS grid
#####
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
  tmp <- aggregate(tmp, fact=8, fun=mean, expand=T, filename="processed/LULC.240m.test.tif", format="GTiff") ## aggregate to roughly the MODIS grid size by fraction of each LULC class
  resample(tmp, evi.r, method="bilinear",   ### get the mean fractional area resampled to the slightly bigger grid
                filename=paste("processed/LULC.250m.", fields[g], ".frac.tif", sep=""),
                format="GTiff",
                overwrite=T)
  print(paste("finished aggregating", fields[g]))
}


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
#####


# ###Process NAIP 1m CIR data to NDVI
#####
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
#####

