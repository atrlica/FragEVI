library(raster)
library(data.table)
library(rgeos)
library(rgdal)
library(sp)
library(ggplot2)
library(knitr)

### Stack 30m layers and crop to AOI
dat.r <- stack("processed/EVI30m.stack.tif")
names(dat.r) <-  c("evi", "isa", "lulc", "AOI")
dat <- as.data.table(as.data.frame(dat.r))
dat <- dat[AOI==1,]
dat <- dat[is.finite(lulc),] ## 734k pixels with valid EVI/ISA/LULC

### Zhao framework to find Vzi

## alternate approach: just chuck the water values which will fuck everything up (v. low evi, v. low isa)
veg <- dat[isa<0.01  & lulc!=20, median(evi, na.rm=T)]
veg.err <- dat[isa<0.01  & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<0.01  & lulc!=20, length(evi)]
noveg <- dat[isa>=0.99 & lulc!=20, median(evi, na.rm=T)]
noveg.err <- dat[isa>=0.99 & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=0.99 & lulc!=20, length(evi)]

### get range of intensity values and figure a linear Vzi curve
beta.range <- seq(from=0, to=1, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)
dat[lulc!=20, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[lulc!=20, .(m.isa=median(isa, na.rm=T), m.evi=median(evi, na.rm=T)), by=bin]
evi.bin <- evi.bin[order(m.isa),]

### binned EVI vs. ISA, all LULC
par(mar=c(4,4,2,1), mgp=c(2,1,0), mfrow=c(1,2))
pdf(file = "images/EVI_ISA_binned30m.pdf", width = 6, height = 6)
plot(evi.bin[,m.isa]*100, evi.bin[,m.evi],
     main="EVI vs. %ISA, all LULC",
     xlab="%ISA", ylab="EVI",
     pch=15, col="forestgreen")
lines(beta.range*100, Vzi, col="red", lwd=2, lty=2)
plot(evi.bin[,m.isa]*100, evi.bin[,m.evi]-Vzi,
     main="EVI enhancement, all LULC",
     xlab="%ISA", ylab="EVI enhancement",
     pch=16, col="royalblue")
abline(h=0, lty=1, lwd=1.6)
dev.off()

#### look at EVI vs. ISA plot by subtype
dat[,lulc:=(as.factor(lulc))]
vtypes <- dat[,unique(lulc)]
vtypes <- vtypes[!is.na(vtypes)]
vtypes <- vtypes[!vtypes==20]
vnum <- c(1:20, 23:26, 29, 31, 34:40)
vnames <- c("Crop", "Pasture", "Forest", "NFwet", "Mining", "Open", "PartRec", 
            "SpectRec", "WBRec", "MFResid", "HDResid", "MDResid", "LDResid", "SWwet",
            "Comm", "Ind", "Tran", "Transp", "Waste", "Water", "CBog", "Util", 
            "SWBeach", "Golf", "Marina", "PubInst", "Cem", "Orch", "Nurs", "FWet", "VLDResid", 
            "Junk", "Brush")
vlook <- data.frame(cbind(vnum, vnames))
pal <- rainbow(length(vtypes))
for(v in 1:length(vtypes)){
  veg <- dat[isa<0.01 & lulc==vtypes[v], median(evi, na.rm=T)]
  noveg <- dat[isa>=0.99 & lulc==vtypes[v], median(evi, na.rm=T)]
  if(is.na(noveg)){noveg <- dat[isa>=0.90 & lulc==vtypes[v], median(evi, na.rm=T)]} ## relax assumptions for lulc that doesn't cover the whole range
  if(is.na(veg)){veg <- dat[isa<=0.10 & lulc==vtypes[v], median(evi, na.rm=T)]}
  Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)
  evi.tmp <- dat[lulc==vtypes[v],]
  evi.tmp[, bin:=findInterval(isa, beta.range, all.inside=T)]
  plotme <- evi.tmp[, .(m.isa=median(isa, na.rm=T), m.evi=median(evi, na.rm=T)), by=bin]
  plot(plotme[,m.isa]*100, plotme[,m.evi],
       main=paste("EVI vs. ISA,", vlook[vnum==vtypes[v], "vnames"]),
       xlab="%ISA", ylab="EVI",
       pch=15, col=pal[v],
       ylim=c(0, 0.65), xlim=c(0,100))
  lines(beta.range*100, Vzi, col="red", lwd=2, lty=2)
}

### fractional cover by LULC
vfracs <- dat[,length(evi), by=lulc]
vfracs.t <- vfracs
vfracs.t$lulc.name <- NA
for(h in 1:dim(vfracs.t)[1]){
  hunt <- vfracs.t$lulc[h]
  vfracs.t$lulc.name[h] <- as.character(vlook$vnames[vlook$vnum==hunt])
}
vfracs.t$Tfrac <- round(vfracs$V1/dat[!is.na(AOI),length(evi)], digits=3)
keep <- vfracs.t[order(Tfrac, decreasing = T),]
keep <- keep[Tfrac>=0.01,]
sum(keep$Tfrac)



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





#########
#### Analysis of 30m Boston data with 1m cover fractions
dat.r <- stack("processed/boston/")

