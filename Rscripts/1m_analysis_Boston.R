library(raster)
library(data.table)
library(rgeos)
library(rgdal)
library(sp)
library(ggplot2)
library(knitr)


################
### analysis of 1m data with edge classes, select zones of Boston city limits
bos.stack <- stack("processed/bos.stack.1m.tif")
names(bos.stack) <- c("can", "ndvi", "cov", "isa", "ed10", "ed20", "ed30")

### load up test polygons for different areas of the city
rox <- readOGR(dsn = "processed/zones/roxbury_test1.shp", layer="roxbury_test1")
sb <- readOGR(dsn = "processed/zones/stonybrook_test1.shp", layer="stonybrook_test1")
allan <- readOGR(dsn = "processed/zones/allandale_test1.shp", layer="allandale_test1")
com <- readOGR(dsn = "processed/zones/commons_test1.shp", layer="commons_test1")
matt <- readOGR(dsn = "processed/zones/mattapan_test1.shp", layer="mattapan_test1")
golf <- readOGR(dsn = "processed/zones/gwrightgolf_test1.shp", layer="gwrightgolf_test1")
ncc <- readOGR(dsn = "processed/zones/newcalvcem_test1.shp", layer="newcalvcem_test1")
wrox <- readOGR(dsn = "processed/zones/wroxbury_test1.shp", layer="wroxbury_test1")
dor <- readOGR(dsn = "processed/zones/dorchester_test1.shp", layer="dorchester_test1")
beach <- readOGR(dsn = "processed/zones/beachmont_test1.shp", layer="beachmont_test1")
send <- readOGR(dsn = "processed/zones/southend_test1.shp", layer="southend_test1")
allst <- readOGR(dsn = "processed/zones/allston_test1.shp", layer="allston_test1")

objs <- list(rox, sb, allan, com, matt, golf, ncc, wrox, dor, beach, send, allst)
tests <- c("Roxbury", "Stonybrook", "Allandale", "Commons", "Mattapan", "GWGolf", "NewCalvCem", "WRoxbury", "Dorchester", "Beachmont", "Southend", "Allston")
types <- c("residential", "forest", "forest", "park", "residential", "golf", "cemetery", "residential", "residential", "non-forest", "residential", "residential")
areas <- numeric()
g <- data.frame()
y <- data.frame()
print("starting loop through test zones")
### loop the test polygons, extract and produce summary NDVI statistics
for(m in 3:length(objs)){
  dat <- extract(bos.stack, objs[[m]], df=T)
  dat <- as.data.table(dat)
  a <- round(dim(dat)[1]/10000, 1) ## area in ha
  e10.sum <- round(dat[cov==2 & ed10==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(dat[cov==2])[1])], 4)
  if(dim(e10.sum)[1]==0){e10.sum <- cbind(0, 0)}
  e20.sum <- round(dat[cov==2 & is.na(ed10) & ed20==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(dat[cov==2])[1])],4)
  if(dim(e20.sum)[1]==0){e20.sum <- cbind(0, 0)}
  e30.sum <- round(dat[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, .(median(ndvi, na.rm=T), length(ndvi)/dim(dat[cov==2])[1])], 4)
  if(dim(e30.sum)[1]==0){e30.sum <- cbind(0, 0)}
  e40.sum <- round(dat[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), .(median(ndvi, na.rm=T), length(ndvi)/dim(dat[cov==2])[1])], 4)
  if(dim(e40.sum)[1]==0){e40.sum <- cbind(0, 0)}
  can.a <- round(dim(dat[cov==2,])[1]/dim(dat)[1], 2)
  grass.sum <- round(dat[cov==1, .(median(ndvi, na.rm=T), dim(dat[cov==1,])[1]/dim(dat[,])[1])], 4)
  if(dim(grass.sum)[1]==0){grass.sum <- cbind(0, 0)}
  imp.sum <- round(dat[cov==0 & isa==1, .(median(ndvi), length(ndvi)/dim(dat[cov==0,])[1])], 4) ## ndvi of isa, relative fraction of barren
  if(dim(imp.sum)[1]==0){imp.sum <- cbind(0, 0)}
  barr.sum <- round(dat[cov==0 & isa==0, .(median(ndvi), length(ndvi)/dim(dat[cov==0,])[1])], 4) ## ndvi of non-impervious barren, relative fraction of barren
  if(dim(barr.sum)[1]==0){barr.sum <- cbind(0, 0)}
  isa.a <- round(dim(dat[cov==0,])[1]/dim(dat)[1],2)
  
  d <- unname(as.vector(c(cbind(a, e10.sum, e20.sum, e30.sum, e40.sum, can.a, grass.sum, imp.sum, barr.sum, isa.a))))
  g <- rbind(g, d)
  print(paste("finished summary of zone", m))
  
  dat$zone <- as.factor(m)
  y <- rbind(y, dat)
  print(paste("integrated data from zone ", m))
}
colnames(g) <- c("Area(ha)", "NDVI.10m", "10m.frac", "NDVI.20m", "20m.frac", "NDVI.30m", "30m.frac", "NDVI.Cint", "Cint.frac", "Can.Tfrac", "NDVI.grass", "grass.Tfrac", "NDVI.imp", "imp.frac", "NDVI.nonimp", "nonimp.frac", "barren.Tfrac")
b <- cbind(tests, types, g)
colnames(b)[1:2] <- c("Zone", "type")
write.csv(b, "processed/bos.zones.summ.csv")

####
### clustered boxplots
## make new class scheme including forest and barren positions
### this combined data frame may be too large to work with (~1M pixels per test area)
y <- as.data.table(y)
y[,cov.ed:=0] ## barren, nonimpervious
y[cov==1, cov.ed:=1] ## grass
y[cov==2 & ed10==1, cov.ed:=2] ## 10m canopy
y[cov==2 & is.na(ed10) & ed20==1, cov.ed:=3] ## 20m canopy
y[cov==2 & is.na(ed10) & is.na(ed20) & ed30==1, cov.ed:=4] ## 30m canopy
y[cov==2 & is.na(ed10) & is.na(ed20) & is.na(ed30), cov.ed:=5] ## >30m canopy
y[cov==0 & isa==1, cov.ed:=6] ## barren, impervious
fill.pal <- c("gray80", "lightgreen", "darksalmon", "orange", "yellow3", "darkgreen", "gray40")
# y$cov.ed <- as.factor(y$cov.ed)
# y[zone==1, zone.name:="Roxbury"]
# y[zone==2, zone.name:="Stonybrook"]
# y[zone==3, zone.name:="Allandale"]
# y[zone==4, zone.name:="Commons"]
# y[zone==5, zone.name:="Mattapan"]
# y[zone==6, zone.name:="GWGolf"]
# y[zone==7, zone.name:="NewCalvCem"]
# y[zone==8, zone.name:="WRoxbury"]
# y[zone==9, zone.name:="Dorchester"]
# y[zone==10, zone.name:="Beachmont"]
# y[zone==11, zone.name:="Southend"]
# y[zone==12, zone.name:="Allston"]
# y$zone.name <- as.factor(y$zone.name)
# write.csv(y, "processed/Bos.tests.extract.data.csv")
p <- ggplot(aes(x=zone, y=ndvi, fill=cov.ed), data=y)+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width=.8))+
  scale_fill_manual(values=fill.pal, labels=c("barren", "grass", "10m", "20m", "30m", ">30m", "imperv."))+
  labs(title="NDVI by edge+cover, 1m", y="NDVI 1m")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_discrete(name="Test Zone", labels=tests)
p




