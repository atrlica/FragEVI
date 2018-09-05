library(data.table)
library(raster)
library(rgdal)
library(rgeos)


#### pull in ACES FF CO2 data
## CentFL ACRU: lin: a+b*age-or-dbh
a=1.516
b=0.414
a+(b*25) ## 11.87 --> this equation is age~dbh, so this is a 12 year old tree

biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- readOGR("data/towns/TOWNS_POLY.shp")
aoi <- spTransform(aoi, crs(biom))
aoi.bos <- aoi[aoi@data$TOWN%in%c("BOSTON", 
                                  "BROOKLINE", 
                                  "NEWTON", 
                                  "CAMBRIDGE", 
                                  "EVERETT", 
                                  "CHELSEA", 
                                  "SOMERVILLE",
                                  "WATERTOWN", 
                                  "MILTON", 
                                  "QUINCY",
                                  "DEDHAM",
                                  "NEEDHAM"),]
plot(aoi.bos)
plot(biom)
plot(aoi.bos, add=T)

### this is the geotif of annual (2011?) FF emissions in kgC/km2/yr
ff <- raster("data/ACES/sdat_1501_12_20180720_100404481.tif")
ff <- projectRaster(ff, crs=crs(biom))
MA.tot <- extract(ff, aoi)
(sum(unlist(MA.tot))/1E9)*44/12 ### kgC/km2/yr-->MMTCO2/yr, 44 for MA, checks out

ff.bos <- crop(ff, aoi.bos)
plot(ff.bos); plot(aoi.bos, add=T)
bud <- unlist(extract(ff.bos, aoi.bos[aoi.bos@data$TOWN=="BOSTON",]))
hist(bud)
sum(bud)/1E3
# a <- values(ff.bos)
# a[a>1E08] <- NA
# ff.cor <- raster(ff.bos)
# ff.cor <- setValues(ff.cor, a)
# plot(ff.cor); plot(aoi.bos, add=T)

### need netcdf files of HOURLY (2013-2014)to get a July month extract
### see https://thredds.daac.ornl.gov/thredds/catalog/ornldaac/1501/catalog.html, 7 GB for each year

