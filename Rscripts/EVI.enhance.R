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
dat <- read.csv("processed/boston.30m.agg.csv")
dat <- as.data.table(dat)
dat <- dat[!is.na(aoi) & water<=0.01,] # 137 k pixels inside AOI and nearly water free
dat <- dat[aoi>675,] # only take pixels that have >75% sub-pixel data == 134k pixels
# hist(dat[,evi])
# hist(dat[,ndvi])

### look at values binned by evi range
bin.range <- seq(from=dat[,min(evi, na.rm=T)], to=dat[,max(evi, na.rm=T)], length.out = 100)
dat[,evi.bin:=findInterval(evi, bin.range, all.inside=T)]
evi.b <- dat[,.(m.barr=mean(barr),
                m.grass=mean(grass),
                m.can=mean(can),
                m.isa=mean(isa),
                m.nonimpbarr=mean(nonimpbarr),
                m.vegisa=mean(vegisa),
                m.ed10=mean(ed10),
                m.ed20=mean(ed20),
                m.ed30=mean(ed30),
                m.buff.20only=mean(buff.20only),
                m.buff.30only=mean(buff.30only),
                m.buff.Intonly=mean(buff.Intonly),
                m.forest=mean(forest),
                m.lowveg=mean(lowveg),
                m.low.res=mean(low.res),
                m.med.res=mean(med.res),
                m.hd.res=mean(hd.res),
                m.dev=mean(dev),
                m.other=mean(other),
                m.ndvi=mean(ndvi), 
                m.evi=mean(evi, na.rm=T),
                n.evi=length(evi, na.rm=T)), by=evi.bin] 
### basic plots -- might be good to display bin size 
### should look at evi binned according to each sub-pixel area factor (rather than binned to evi)
plot(evi.b$m.evi, evi.b$m.ndvi) ### ndvi looks non-informative below 0 EVI and above 0.7 EVI
plot(evi.b$m.evi, evi.b$m.barr) ## above ~95% barren all EVI is <0.05, nice line down to ~3% baren then a large range in EVI (0.6-0.8)
plot(evi.b$m.evi, evi.b$m.grass) ## weird, response between 0-0.3 grass, then 0.5-0.8 EVI covers whole range from 0.2-0.8 grass
plot(evi.b$m.evi, evi.b$m.can) ## 0-0.4 EVI ~ 0-0.5 can, then above 0.6 EVI can be anything from 0.2-0.8 can
plot(evi.b$m.evi, evi.b$m.isa) ## relatively diagnostic: anything below 0 EVI is 100% ISA, anything above 0.7 EVI is 0-10% ISA -- greenest parts of the city contain some pavement
plot(evi.b$m.evi, evi.b$m.nonimpbarr) ## bell curve: low EVI~low NIB, peak NIB at 0.3 EVI, EVI>0.7 ~ 0 NIB
plot(evi.b$m.evi, evi.b$m.vegisa) # similar to NIB, peak at VISA 0.2 ~ EVI 0.4, highest but sloppy EVI at VISA 0-7% -- the greenest parts of the city have some ISA below the canopy
plot(evi.b$m.evi, evi.b$m.ed10) ## up to 0.6 ed10 area, plateau ed10 (0.4) EVI 0.4-0.6, then decline in ed10 up till EVI>0.7, which covers 0.1 to 0.6 ed10 (greenest pixels can contain relatively little or a lot of edge canopy)
plot(evi.b$m.evi, evi.b$m.ed20) ## similar to ed10 but shifted right, EVI>0.7 can be 0.2-0.7 ed10
plot(evi.b$m.evi, evi.b$m.ed30) # linearish up to EVI0.6, then range of 0.2-0.7 ed30
plot(evi.b$m.evi, evi.b$m.buff.20only)
plot(evi.b$m.evi, evi.b$m.buff.30only)
plot(evi.b$m.evi, evi.b$m.buff.Intonly) ## interior forest peaks 0.6-0.8EVI, but again greenest pixels are 0% interior
plot(evi.b$m.evi, evi.b$m.forest) ## forest 0.4-0.8 for EVI >0.6, but EVI~0.8 is forest 0.4-1.0
plot(evi.b$m.evi, evi.b$m.lowveg) ## very sloppy, generall EVI>0.7 can be 0-0.5 lowveg
plot(evi.b$m.evi, evi.b$m.low.res) ## noise, too little range covered
plot(evi.b$m.evi, evi.b$m.med.res) ## peak 0.2-0.6EVI, but low coverage (<0.05 med.res)
plot(evi.b$m.evi, evi.b$m.hd.res) ## smooth bell curve 0-0.6EVI, but above 0.6 very low hd.res
plot(evi.b$m.evi, evi.b$m.dev) ### low control: 100% dev can be EVI0-0.2, EVI>0.4 can be 0-20% dev
plot(evi.b$m.evi, evi.b$m.other) ## chaos, too little coverage

plot(evi.b$m.grass, evi.b$m.can)

diag <- evi.b[m.evi>0.6, ]
diag <- diag[order(m.evi),]
View(diag)

d <- lm(m.evi~m.grass*m.can, data=evi.b); summary(d) ## R2 0.96, no interaction effect

### evi as response in intervals binned by the cover area
library(ggplot2)

## canopy
bin.range <- seq(from=dat[,min(can, na.rm=T)], to=dat[,max(can, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(can, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.can=mean(can, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.can, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="Canopy")
# uncomplicated except tails -- can>80% shows steeper curve up to 0.6EVI, and can<0.03 EVI range 0.1-0.2

## isa
bin.range <- seq(from=dat[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(isa, bin.range, all.inside=T)]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.isa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="Canopy")
# classic decline, falls steeply above first step at 0, bellies out at higher ISA, crashes again at end (looks like full AOI)


## grass
bin.range <- seq(from=dat[,min(grass, na.rm=T)], to=dat[,max(grass, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(grass, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.grass=mean(grass, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.grass, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="Grass")
# grass <10% can be EVI 0.3-0.2, somewhat noisy, tops out at grass>0.8 but EVI<0.5 (ok so even very grassy pixels don't get as high as canopy)

## barren
bin.range <- seq(from=dat[,min(barr, na.rm=T)], to=dat[,max(barr, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(barr, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.barr=mean(barr, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.barr, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="Barren")
# fully barren is 0.1EVI, fully not is EVI 0.6, steady curve in between -- looks a lot like ISA vs. EVI curve, discontinuities at either end (must be very different members in those groups, lots more pixels)


## ed10
bin.range <- seq(from=dat[,min(ed10, na.rm=T)], to=dat[,max(ed10, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(ed10, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.ed10=mean(ed10, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.ed10, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="ed10")
# # anemic curve, mainly 0 ed10, saturating at EVI0.4~ed10 0.6, a few high ed10 with EVI close .45 and up

## ed20
bin.range <- seq(from=dat[,min(ed20, na.rm=T)], to=dat[,max(ed20, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(ed20, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.ed20=mean(ed20, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.ed20, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="ed20")
# linear increase up to ed20 0.9 ~ EVI0.5


## ed30
bin.range <- seq(from=dat[,min(ed30, na.rm=T)], to=dat[,max(ed30, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(ed30, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.ed30=mean(ed30, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.ed30, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="ed30")
# linear increase up to ed20 1.0 ~ EVI0.6

## buff.Intonly
bin.range <- seq(from=dat[,min(buff.Intonly, na.rm=T)], to=dat[,max(buff.Intonly, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(buff.Intonly, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.buff.Intonly=mean(buff.Intonly, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.buff.Intonly, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="buff.Intonly")
# not diagnostic: 75% are buff.Intonly = 0, all EVI with any interior canopy are >0.6, trending slightly up with greater interior canopy

## nonimpbarr (NIB)
bin.range <- seq(from=dat[,min(nonimpbarr, na.rm=T)], to=dat[,max(nonimpbarr, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(nonimpbarr, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.nonimpbarr=mean(nonimpbarr, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.nonimpbarr, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="nonimpbarr")
# very sparse above NIB>0.3, negative trend but overall little control on EVI


## vegisa (VISA)
bin.range <- seq(from=dat[,min(vegisa, na.rm=T)], to=dat[,max(vegisa, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(vegisa, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.vegisa=mean(vegisa, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.vegisa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="VISA")
# sparse above VISA 0.5, stalls at EVI 0.35 but scattered very high EVI (0.6+ with very high VISA)

## forest
bin.range <- seq(from=dat[,min(forest, na.rm=T)], to=dat[,max(forest, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(forest, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.forest=mean(forest, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.forest, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="forest")
# very sparse, 88%% have almost 0 forest, range EVI 0.45-0.6+, trending slightly upward -- i.e. any forest exposure tends to be quite high EVI


## hd.res
bin.range <- seq(from=dat[,min(hd.res, na.rm=T)], to=dat[,max(hd.res, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(hd.res, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.hd.res=mean(hd.res, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.hd.res, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="hd.res")
# essentially bimodal: 50% no hd.res = 0.285EVI, 29% full hd.res = 0.31EVI, no signal on EVI -- close to AOI mean of 0.29 EVI

## lowveg
bin.range <- seq(from=dat[,min(lowveg, na.rm=T)], to=dat[,max(lowveg, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(lowveg, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.lowveg=mean(lowveg, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.lowveg, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="lowveg")
# too rare -- 91% no lowbveg, basically mean EVI -- tiny upward trend with more lowveg

## dev
bin.range <- seq(from=dat[,min(dev, na.rm=T)], to=dat[,max(dev, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(dev, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.dev=mean(dev, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.dev, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="dev")
# bimodal: 47% no dev, EVI 0.38; 100% dev 0.18, slight downward trend


### pure pixels analysis
forest.p <- dat[forest>0.95, ]
dev.p <- dat[dev>0.95, ]
hdres.p <- dat[hd.res>0.95, ]
lowveg.p <- dat[lowveg>0.95, ]

forest.p[,median(evi, na.rm=T)]#0.64
dev.p[,median(evi, na.rm=T)]#0.13     
hdres.p[,median(evi, na.rm=T)]#0.31
lowveg.p[,median(evi, na.rm=T)]#0.48

## these pure classes account for 75% of Boston (about 25% is mixed pixels)
(forest.p[,length(evi)]+lowveg.p[,length(evi)]+hdres.p[,length(evi)]+dev.p[,length(evi)])/dim(dat)[1] 
(forest.p[,length(evi)])/dim(dat)[1] # pure forest is 5.4%
(lowveg.p[,length(evi)])/dim(dat)[1] # pure lowveg is 5.0%
(hdres.p[,length(evi)])/dim(dat)[1] # pure hd.res is 30.8%
(dev.p[,length(evi)])/dim(dat)[1] # pure dev is 34.0%

### what is relationship with isa etc. inside pure pixel classes?

### pure forest 7k
bin.range <- seq(from=forest.p[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
a <- forest.p
a[,bin:=findInterval(isa, bin.range, all.inside=T)]
a[,visafrac:=0]
a[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
a <- a[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            count=.N,
            m.visa=mean(vegisa, na.rm=T),
            m.can=mean(can, na.rm=T),
            m.visafrac=mean(visafrac, na.rm=T)), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.isa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="isa in Forest")
# 76% of pure forest is 0 ISA, stable slightly lower from 0-40% ISA (EVI0.55), drops chaotically to ~0.25 above 40% ISA -- so lightly paved forests look a lot like other forest, but the weird highly paved look more like developed
ggplot(dog, aes(x=m.isa, y=m.visa))+geom_point(aes(size=count.frac))+labs(title="VISA in Forest")
## steady increase in VISA up till 0.4 ISA 20%VISA, then gets chaotic -- up to 40-60% VISA in above 0.5 ISA
ggplot(dog, aes(x=m.isa, y=m.can))+geom_point(aes(size=count.frac))+labs(title="Canopy in Forest")
## nearly 100% canopy, sub-linear drop in canopy with ISA (i.e. not trading canopy for increasing ISA at 1:1)
### canopy in forest stays intact with mild paving then drops, but almost no "forest" is paved in general (i.e. "forest" really does mean nearly not developed)
ggplot(dog, aes(x=m.isa, y=m.visafrac))+geom_point(aes(size=count.frac))+labs(title="Fraction of canopy over ISA, Forest")
### nearly linear 1:1 (slope slightly <1) -- for every increment of ISA you make the remaining forest canopy more dominated by VISA
ggplot(dog, aes(x=m.can, y=m.evi))+geom_point(aes(size=count.frac, color=m.isa))+labs(title="EVI vs. canopy, FOREST")
### basically forest is all 100% unpaved, paving it at all causes instant decline in EVI


## pure dev, 45k
bin.range <- seq(from=dev.p[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
a <- dev.p
a[,bin:=findInterval(isa, bin.range, all.inside=T)]
a[,visafrac:=0]
a[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
a <- a[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            count=.N,
            m.visa=mean(vegisa, na.rm=T),
            m.can=mean(can, na.rm=T),
            m.visafrac=mean(visafrac, na.rm=T)), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.isa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="isa in Dev")
# simple linear decline in EVI with ISA, except EVI in full ISA drops off a cliff; dev is 46% 1.0ISA EVI0.08, 7% 0 ISA EVI0.41, sparse in between but thicker above ISA 0.5 -- so looking like grass when unpaved but linear decline to basically pavement
ggplot(dog, aes(x=m.isa, y=m.visa))+geom_point(aes(size=count.frac))+labs(title="VISA in Dev")
## VISA 0.65 ISA 12%VISA, then increasing ISA results in lower VISA
ggplot(dog, aes(x=m.isa, y=m.can))+geom_point(aes(size=count.frac))+labs(title="Canopy in Dev")
## canopy varies 0.15-0.3 from ISA 0.0-0.5, then declines
ggplot(dog, aes(x=m.isa, y=m.visafrac))+geom_point(aes(size=count.frac))+labs(title="Fraction of canopy over ISA, DEV")
### fraction of canopy over ISA increases up to peak about VISAfrac0.45 isa 0.9, crashes at highest ISA where very little veg left.
ggplot(dog, aes(x=m.can, y=m.evi))+geom_point(aes(size=count.frac, color=m.isa))+labs(title="EVI vs. canopy, DEV")
## bulk of observations are eighter low can/low EVI, peak EVI at 0.3/CAN0.3, but highest EVI at CAN 0.14 -- it's where ISA is 0
### so... you add canopy from ISA zero and the EVI AND ISA DECLINES. -- could be a gremlin, not very many zero-to-low ISA pixels
ggplot(dog, aes(x=m.can, y=m.visafrac))+geom_point(aes(size=count.frac, color=m.isa))+labs(title="EVI vs. VISAFRAC, DEV")


## pure hd.res ### not very many, 41k
bin.range <- seq(from=hdres.p[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
a <- hdres.p
a[,bin:=findInterval(isa, bin.range, all.inside=T)]
a[,visafrac:=0]
a[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
a <- a[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            count=.N,
            m.visa=mean(vegisa, na.rm=T),
            m.can=mean(can, na.rm=T),
            m.visafrac=mean(visafrac, na.rm=T)), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.isa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="isa in HDResid")
# nice even distribution of ISA -- sigmoidal, steep decline with increasing ISA but about 0.4 ISA it levels out in its decline some
ggplot(dog, aes(x=m.isa, y=m.visa))+geom_point(aes(size=count.frac))+labs(title="VISA in HDResid")
## vegetation over ISA peaks about 0.65 ISA 22%VISA, then increasing IS results in lower VISA
ggplot(dog, aes(x=m.isa, y=m.can))+geom_point(aes(size=count.frac))+labs(title="Canopy in HDResid")
## classic decline from 0.8 can/0ISA to 1.0ISA/0 canopy, sublinear with steep drop only at the very end (above ISA 90%)
## HD Resid has nearly as much canopy as a forest at low ISA, maintains some canopy cover almost up to total ISA
ggplot(dog, aes(x=m.isa, y=m.visafrac))+geom_point(aes(size=count.frac))+labs(title="Fraction of canopy over ISA, HDResid")
### nearly perfectly linear increase in VISAfrac up to 80%/90%ISA, then crashes in once ISA gets too high

ggplot(dog, aes(x=m.can, y=m.evi))+geom_point(aes(size=count.frac, color=m.isa))+labs(title="EVI vs. canopy, HDRESID")
### so peak EVI is not as high as forest (forest is nearly all above 0.55), but peak canopy is nearly as  high. Back off from full
## canopy and EVI drops initially but then it slows down in the middle -- VISAfrac is increasing up till you fall to 0.2 canopy
## after 0.2-0.3 canopy, you're playing for scraps and EVI declines quickly again
ggplot(dog, aes(x=m.can, y=m.visafrac))+geom_point(aes(size=count.frac, color=m.isa))+labs(title="CAN VS VISAFRAC, HDRESID")
## as canopy declines a greater fraction of it is in VISA until very low canopy, after which the remaining trees must be in a yard or something


## pure lowveg 6.6k
bin.range <- seq(from=lowveg.p[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
a <- lowveg.p
a[,bin:=findInterval(isa, bin.range, all.inside=T)]
a[,visafrac:=0]
a[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
a <- a[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            count=.N,
            m.visa=mean(vegisa, na.rm=T),
            m.can=mean(can, na.rm=T),
            m.visafrac=mean(visafrac, na.rm=T)), by=bin]
dog[,count.frac:=count/sum(count)]
ggplot(dog, aes(x=m.isa, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="isa in lowveg")
# 51% are 0 ISA, get very rare above ISA 0.35, but almost no effect on EVI from 0-0.4 ISA, then steep decline -- wtf it's the opposite of hd.res, initial ISA additions don't really degrade EVI
ggplot(dog, aes(x=m.isa, y=m.visa))+geom_point(aes(size=count.frac))+labs(title="VISA in lowveg")
## VISA rises to about 0.25/ISA40%, then gets chaotic
ggplot(dog, aes(x=m.isa, y=m.can))+geom_point(aes(size=count.frac))+labs(title="Canopy in lowveg")
## low canopy at 0ISA (15%can), weak peak at~ISA40%
ggplot(dog, aes(x=m.isa, y=m.visafrac))+geom_point(aes(size=count.frac))+labs(title="Fraction of canopy over ISA, VISA")
### VISA fraction increases up to ISA 25%, then steady, chaotic above 50% ISA

ggplot(dog, aes(x=m.can, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="EVI vs. canopy, lowveg")
## bulk of observations are can 0.15, .46 EVI -- other high EVI are in higher bins, disorganized



