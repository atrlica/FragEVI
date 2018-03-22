library(raster)
library(data.table)
library(rgeos)
library(rgdal)
library(sp)
library(ggplot2)
library(knitr)

#### 30m EVI enhancement assessment
### Stack 30m layers and crop to AOI
dat.r <- stack("processed/EVI30m.stack.tif")
names(dat.r) <-  c("evi", "isa", "lulc", "AOI")
dat <- as.data.table(as.data.frame(dat.r))
dat <- dat[AOI==1,]
dat <- dat[is.finite(lulc),] ## 7.3M pixels with valid EVI/ISA/LULC

## get values for vegetation and non-veg endmembers
veg <- dat[isa<0.01  & lulc!=20, median(evi, na.rm=T)]
veg.err <- dat[isa<0.01  & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<0.01  & lulc!=20, length(evi)]
noveg <- dat[isa>=0.99 & lulc!=20, median(evi, na.rm=T)]
noveg.err <- dat[isa>=0.99 & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=0.99 & lulc!=20, length(evi)]
beta.range <- seq(from=0, to=.99, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)

## bin the values by isa fraction
dat[lulc!=20, bin:=findInterval(isa, beta.range, all.inside=F)]
evi.bin <- dat[lulc!=20, .(m.isa=median(isa, na.rm=T), 
                           m.evi=median(evi, na.rm=T),
                           bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]

### binned EVI vs. ISA, all LULC
# pdf(file = "images/EVI_ISA_binned30m.pdf", width = 6, height = 6)
# plot(evi.bin[,m.isa]*100, evi.bin[,m.evi],
#      main="EVI vs. %ISA, all LULC",
#      xlab="%ISA", ylab="EVI",
#      pch=15, col="forestgreen")
# lines(beta.range*100, Vzi, col="red", lwd=2, lty=2)
# plot(evi.bin[,m.isa]*100, evi.bin[,m.evi]-Vzi,
#      main="EVI enhancement, all LULC",
#      xlab="%ISA", ylab="EVI enhancement",
#      pch=16, col="royalblue")
# abline(h=0, lty=1, lwd=1.6)
# dev.off()

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="blue")+
  scale_size_continuous(range=c(2, 8), breaks=c(0.01, 0.02, 0.40), name="% total")+
  labs(title="EVI vs. ISA by bin, 30m", x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)

ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="blue")+
  labs(title="EVI enhancement", x="Frac. ISA", y="EVI enhancement")+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)
## huge concentration of data in 0% ISA category (~44%)

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

## do EVI vs ISA coloring the points according to majority area
lulc.look <- data.frame(class=c(1,2,3,4,5,6,7,0))
lulc.look$name <- c("low.res",
                    "med.res",
                    "hi.res",
                    "dev",
                    "forest",
                    "lowveg",
                    "other",
                    "none")
lulc.look$col <- c("yellow", "orange", "salmon", "grey45", "forestgreen", "lightgreen", "red", "black")
lu.lowres <- c(13, 38) #LDResid., VLDResid.
lu.medres <- c(12) # MDResid.
lu.hdres <- c(11,10) # HDResid., MFResid.,
lu.dev <- c(15,16,7,8,18,39,24,31,19,9,29) # Comm, Ind, PartRec, SpectRec, Transp., Junk, Util, PubInst, Waste, WBRec, Marina
lu.forest <- c(3,37) # Forest, FWet
lu.lowveg <- c(4,1,2,6,5,26,14,25,34,23) # NFwet, PubInst, Crop, Past, Open, Mining, Golf, SWwet, SWBeach, Cem, CBog
lu.other <- c(17,35,40,36) #Tran, Orch, Brush, Nurs

dat[lulc%in%lu.lowres, lulc.class:=1]
dat[lulc%in%lu.medres, lulc.class:=2]
dat[lulc%in%lu.hdres, lulc.class:=3]
dat[lulc%in%lu.dev, lulc.class:=4]
dat[lulc%in%lu.forest, lulc.class:=5]
dat[lulc%in%lu.lowveg, lulc.class:=6]
dat[lulc%in%lu.other, lulc.class:=7]

for(m in 1:100){
  tot.m <- dat[bin==m, .N]
  uu <- dat[bin==m, .N/tot.m, by=lulc.class]
  hit <- uu[V1==max(V1), lulc.class] ## most common pixel class, not necessarily majority
  # if(uu[lulc.class==hit,V1]<0.5){hit <- 0} ## mixed pixels,  no class >0.5
  evi.bin[bin==m, lulc.maj:=hit]
}
evi.bin <- merge(evi.bin, lulc.look, by.x="lulc.maj", by.y="class") ## add lulc.class names and colors
evi.bin <- evi.bin[order(bin),]
evi.bin$lulc.maj <- as.factor(evi.bin$lulc.maj)

## get values for vegetation and non-veg endmembers
veg <- dat[isa<0.01  & lulc!=20, median(evi, na.rm=T)]
veg.err <- dat[isa<0.01  & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<0.01  & lulc!=20, length(evi)]
noveg <- dat[isa>=0.99 & lulc!=20, median(evi, na.rm=T)]
noveg.err <- dat[isa>=0.99 & lulc!=20, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=0.99 & lulc!=20, length(evi)]
beta.range <- seq(from=0, to=1, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)

#EVI vs. ISA plot, color by majority area fraction
ggplot(aes(x=m.isa, y=m.evi, color=lulc.maj), data=evi.bin)+
  geom_point(aes(size=bin.frac))+
  scale_color_manual(breaks=c("1", "2", "3", "4", "5"),
                     values=c("yellow", "orange", "salmon", "grey45", "forestgreen"),
                     guide=guide_legend(title="largest LULC"),
                     labels=c("VLD+LD Resid.", "MD Resid.", "HD+MF Resid.", "Developed", "Forest"))+
  labs(title="EVI vs. binned ISA, 30m AOI", x="Frac. ISA", y="Median EVI")+
  scale_size(range=c(2,8),breaks=c(0.01, 0.02, 0.40),
             labels=c(">1%","1-2%", ">40%"),
             guide=guide_legend(title="% total"))+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100], color=lulc.maj), data=evi.bin)+
  geom_point(aes(size=bin.frac))+
  scale_color_manual(breaks=c("1", "2", "3", "4", "5"),
                     values=c("yellow", "orange", "salmon", "grey45", "forestgreen"),
                     guide=guide_legend(title="largest LULC"),
                     labels=c("VL+L Resid", "M Resid", "H+MF Resid", "Dev", "For"))+
  labs(title="EVI enhancement, 30m AOI", x="Frac. ISA", y="EVI enhancement")+
  scale_size(range=c(0.7,10),
             breaks=c(0.01, 0.02, 0.40),
             labels=c(">1%","1-2%", ">40%"),
             guide=guide_legend(title="Bin size"))+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)



#########
#### Analysis of 30m Boston data with 1m cover fractions
dat <- read.csv("processed/boston.30m.agg.csv")
dat <- as.data.table(dat)
dat <- dat[!is.na(aoi) & water<=0.01,] # 137 k pixels inside AOI and nearly water free
dat <- dat[aoi>675,] # only take pixels that have >75% sub-pixel data == 134k pixels
# hist(dat[,evi])
# hist(dat[,ndvi])


## do all canopy edge fractions match to total canopy fraction?
frog <- dat[,(ed10+buff.20only+buff.30only+buff.Intonly)-can] 
hist(frog); range(frog, na.rm=T)
frog.n <- rep(0, length(frog))
frog.n[frog>0.005] <- 1
sum(frog.n) #only 73 pixels differ on the %cover canopy vs. canopy edge classes


## do we get sensible figures for vegisa?
dat[can>0, range(vegisa)] # 0 to 98% vegisa (this is not likely)
hist(dat[can>0, vegisa]) # a long tail of >0.6 vegisa wtf
duck <- dat[vegisa>can,] ## no instances of vegisa greater than can
# dat[,vegisa.fucked:=0]
# dat[vegisa>can, vegisa.fucked:=1]
# r.duck <- raster("processed/boston/bos.aoi30m.tif")
# r.duck <- setValues(r.duck, values = dat[,vegisa.fucked])
# r.duck <- writeRaster(r.duck, "processed/boston/vegisa.fucked.tif", format="GTiff", overwrite=T)

# 
# ### look at values binned by evi range
# bin.range <- seq(from=dat[,min(evi, na.rm=T)], to=dat[,max(evi, na.rm=T)], length.out = 100)
# dat[,evi.bin:=findInterval(evi, bin.range, all.inside=T)]
# evi.b <- dat[,.(m.barr=mean(barr),
#                 m.grass=mean(grass),
#                 m.can=mean(can),
#                 m.isa=mean(isa),
#                 m.nonimpbarr=mean(nonimpbarr),
#                 m.vegisa=mean(vegisa),
#                 m.ed10=mean(ed10),
#                 m.ed20=mean(ed20),
#                 m.ed30=mean(ed30),
#                 m.buff.20only=mean(buff.20only),
#                 m.buff.30only=mean(buff.30only),
#                 m.buff.Intonly=mean(buff.Intonly),
#                 m.forest=mean(forest),
#                 m.lowveg=mean(lowveg),
#                 m.low.res=mean(low.res),
#                 m.med.res=mean(med.res),
#                 m.hd.res=mean(hd.res),
#                 m.dev=mean(dev),
#                 m.other=mean(other),
#                 m.ndvi=mean(ndvi), 
#                 m.evi=mean(evi, na.rm=T),
#                 n.evi=length(evi, na.rm=T)), by=evi.bin] 
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
## part of the "EVI enhancement" is that the line (Vzi) is determined by 0 and 100% canopy
## 100% canopy is VERY green; 0% is VERY not-green
## EVI responds less to increasing/decreasing canopy at low range, and takes off above ~90%
## below 25% canopy, loss of canopy is associated with greater EVI decline than between 25-75%
### this curve looks a lot like an inverted ISA/EVI curve

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
# grass does the same thing as canopy down low, but never gets to the hyper-green effect at the end of curve
## i.e. a 100% grass pixel tops out lower than a 100% canopy pixel
## and upeer reaches of curve are less constrained -- grass fraction doesn't exert the same control on EVI as canopy
## and the bulk of boston pixels are below 25% grass

## barren -- this is basically what Zhao was using to get urban area for his analysis
## barren is the residual of grass+canopy
bin.range <- seq(from=dat[,min(barr, na.rm=T)], to=dat[,max(barr, na.rm=T)], length.out = 100)
a <- dat
a[,bin:=findInterval(barr, bin.range, all.inside=T)]
a[,count:=.N, by=bin]
a[, ed10.frac:=0]
a[, vegisa.frac:=0]
a[can>0,ed10.frac:=ed10/can]
a[can>0,vegisa.frac:=vegisa/can]
a[can==0, ed10.frac:=NA]
a[ed10.frac>1, ed10.frac:=NA]
a[can==0, vegisa.frac:=NA]
dog <- a[,.(m.evi=mean(evi, na.rm=T), 
            m.ndvi=mean(ndvi, na.rm=T),
            m.barr=mean(barr, na.rm=T),
            m.can=mean(can, na.rm=T),
            m.isa=mean(isa, na.rm=T),
            m.grass=mean(grass, na.rm=T),
            m.vegisa=mean(vegisa, na.rm=T),
            m.nonimpbarr=mean(nonimpbarr, na.rm=T),
            m.ed10=mean(ed10, na.rm=T),
            m.ed10frac=mean(ed10.frac, na.rm=T),
            m.vegisafrac=mean(vegisa.frac, na.rm=T),
            count=.N), by=bin]
dog[,count.frac:=count/sum(count)]

### estimate Vzi
veg <- dog[bin==1, m.evi]
noveg <- dog[bin==99, m.evi]
beta.range <- seq(from=dog[,min(m.barr)], to=dog[,max(m.barr)], length.out = dog[,length(unique(bin))])
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)
dog <- dog[order(bin),]
dog[,Vzi:=Vzi]

## plots
ggplot(dog, aes(x=m.barr, y=m.evi))+
  geom_point(aes(size=count.frac))+
  labs(title="Barren")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
ggplot(dog, aes(x=m.barr, y=m.evi-Vzi))+
  geom_point(aes(size=count.frac))+
  labs(title="EVI enhancement")+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)
## basically dominated by the two end members -- lots of pixels are either 100% barren or 100% vegetated
## fast decrease in EVI up to ~12%, then linear up to about 75%, then fast loss again above ~75%
## looks like the inverse of canopy, but the breaks are sharper on either end and the distribution isn't as lopsided in canopy

### barrenness predicted quite well by canopy cover
ggplot(dog, aes(x=m.barr, y=m.can))+
  geom_point(aes(size=count.frac))+
  labs(title="Canopy")
## 0% barren are ~70% canopy, then nearly linear decline to 100% barren -- i.e. canopy loss is mainly a 1:1 trade with barrenness
# Basically every gain in barrenness is at a loss of canopy, especially above 0.25 barren -- so there really is a trade between canopy area and barrenness that does NOT drag EVI down in linear fashion
summary(lm(dog$m.can~dog$m.barr)) ### slope -0.73 -- a unit increase in barrenness gives a 0.73 decrease in canopy

ggplot(dog, aes(x=m.barr, y=m.grass))+geom_point(aes(size=count.frac))+labs(title="Grass")
# sloppier, sigmoidal, maxes out at 0.35
## grass makes up a component of non-barren that is always less than 30%, also lost fairly regularly with greater increase in barrenness

## grass is mostly below 25%
dat[grass<0.25, length(evi)]/dim(dat)[1] ## 84% of pixels below 25% grass
dat[grass>0.75, length(evi)]/dim(dat)[1] ## 3% of pixels above 75% grass
dat[can<0.25, length(evi)]/dim(dat)[1] ## 51% of pixels below 25% canopy
dat[can>0.75, length(evi)]/dim(dat)[1] ## 14%% of pixels avove 75% canopy

ggplot(dog, aes(x=m.barr, y=m.ed10frac))+geom_point(aes(size=count.frac))+labs(title="fraction canopy <10m edge, Barren")
### below 25% barren, canopy is above 50% and you start to see some >10m edge classes
# 25-75% barren is linear canopy decline, everything is within 10m
## basically by 25% barren all canopy is within 10m of edge
### not sure I trust this still -- why have low fraction at high range of barren?

### fraction of canopy that is within 10m of edge is very high
ggplot(dog, aes(x=m.can, y=m.ed10frac))+geom_point(aes(size=count.frac))+labs(title="fraction canopy <10m edge, Canopy")
### same idea: only when canopy cover is >0.4 do you get interior canopies
### another thing: Even 100% canopy is still 50% edge (???)




ggplot(dog, aes(x=m.ed10, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="area canopy <10m edge")
ggplot(dog, aes(x=m.ed10frac, y=m.evi))+geom_point(aes(size=count.frac))+labs(title="fraction edge canopy vs EVI")

### plot of EVI vs. Barren with Vzi Line and colored by edge canopy fraction
dat <- read.csv("processed/boston.30m.agg.csv")
dat <- as.data.table(dat)
dat <- dat[!is.na(aoi) & water<=0.01,] # 137 k pixels inside AOI and nearly water free
dat <- dat[aoi>675,] # only take pixels that have >75% sub-pixel data == 134k pixels

## get values for vegetation and non-veg endmembers
veg <- dat[barr<0.01  & water<0.02, median(evi, na.rm=T)]
veg.err <- dat[barr<0.01  & water<0.02, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[barr<0.01 & water<0.02, length(evi)]
noveg <- dat[barr>=0.99 & water<0.02, median(evi, na.rm=T)]
noveg.err <- dat[barr>=0.99 & water<0.02, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[barr>=0.99 & water<0.02, length(evi)]
beta.range <- seq(from=0, to=1, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)

dat[,ed10.frac:=ed10/can]
dat[ed10.frac>1, ed10.frac:=NA]
dat[water<0.02, bin:=findInterval(barr, beta.range, all.inside=F)]
evi.bin <- dat[water<0.02, .(m.barr=median(barr, na.rm=T), 
                           m.evi=median(evi, na.rm=T),
                           m.can=median(can, na.rm=T),
                           m.ed10.frac=median(ed10.frac, na.rm=T),
                           bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.barr),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.barr),]


ggplot(evi.bin, aes(x=m.can, y=m.evi, colour=m.ed10.frac))+
  geom_point(aes(size=bin.frac))+
  labs(title="Canopy vs EVI", x="Frac. Barren", y="Median Evi")+
  scale_color_continuous(low="darkgreen", high="salmon", limits=c(0.45, 1))
ggplot(evi.bin, aes(x=m.barr, y=m.evi, colour=m.ed10.frac))+
  geom_point(aes(size=bin.frac))+
  labs(title="Barren vs EVI, 30m, Boston", x="Barren Frac.", y="Median EVI")+
  scale_color_continuous(low="darkgreen", high="salmon", limits=c(0.45, 1), name="% <10m")+
  scale_size_continuous(name="% total", breaks=c(0.01, 0.02, 0.10), range=c(2, 8))+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)


plot(dog$m.evi, dog$m.ndvi)





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


#######
### pure pixels analysis
#### Analysis of 30m Boston data with 1m cover fractions
dat <- read.csv("processed/boston.30m.agg.csv")
dat <- as.data.table(dat)
dat <- dat[!is.na(aoi) & water<=0.01,] # 137 k pixels inside AOI and nearly water free
dat <- dat[aoi>675,] # only take pixels that have >75% sub-pixel data == 134k pixels
bin.range <- seq(from=dat[,min(isa, na.rm=T)], to=dat[,max(isa, na.rm=T)], length.out = 100)
dat[, isa.bin:=findInterval(isa, bin.range, all.inside=T)] ## set up markers for isa bin ranges
bin.range <- seq(from=dat[,min(barr, na.rm=T)], to=dat[,max(barr, na.rm=T)], length.out = 100)
dat[, isa.bin:=findInterval(barr, bin.range, all.inside=T)] ## set up markers for isa bin ranges


### set up median tables for the different LULC classes and combine for plotting
forest.p <- dat[forest>0.95, ]
forest.p[,visafrac:=0]
forest.p[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
forest.p[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
forest.m <- forest.p[,.(m.evi=mean(evi, na.rm=T), 
                        m.ndvi=mean(ndvi, na.rm=T),
                        m.isa=mean(isa, na.rm=T),
                        m.barr=mean(barr, na.rm=T),
                        m.grass=mean(grass, na.rm=T),
                        count=.N,
                        m.visa=mean(vegisa, na.rm=T),
                        m.can=mean(can, na.rm=T),
                        m.visafrac=mean(visafrac, na.rm=T),
                        m.ed10=mean(ed10, na.rm=T),
                        m.ed20=mean(ed20, na.rm=T),
                        m.ed30=mean(ed30, na.rm=T),
                        m.buff.Intonly=mean(buff.Intonly, na.rm=T)), by=barr.bin]
forest.m[,count.frac:=count/dim(forest.p)[1]]
forest.m[,pure.class:="forest"]
# plot(forest.m$m.isa, forest.m$m.barr)
# plot(forest.m$isa, forest.m$m.evi)
# plot(forest.m$m.barr, forest.m$m.evi)

dev.p <- dat[dev>0.95, ]
dev.p[,visafrac:=0]
dev.p[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
dev.p[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
dev.m <- dev.p[,.(m.evi=mean(evi, na.rm=T), 
                        m.ndvi=mean(ndvi, na.rm=T),
                        m.isa=mean(isa, na.rm=T),
                        m.barr=mean(barr, na.rm=T),
                        m.grass=mean(grass, na.rm=T),
                        count=.N,
                        m.visa=mean(vegisa, na.rm=T),
                        m.can=mean(can, na.rm=T),
                        m.visafrac=mean(visafrac, na.rm=T),
                        m.ed10=mean(ed10, na.rm=T),
                        m.ed20=mean(ed20, na.rm=T),
                        m.ed30=mean(ed30, na.rm=T),
                        m.buff.Intonly=mean(buff.Intonly, na.rm=T)), by=barr.bin]
dev.m[,count.frac:=count/dim(dev.p)[1]]
dev.m[,pure.class:="dev"]
plot(dev.m$m.isa, dev.m$m.evi)
plot(dev.m$m.isa, dev.m$m.barr)
plot(dev.m$m.barr, dev.m$m.evi)

hdres.p <- dat[hd.res>0.95, ]
hdres.p[,visafrac:=0]
hdres.p[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
hdres.p[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
hdres.m <- hdres.p[,.(m.evi=mean(evi, na.rm=T), 
                  m.ndvi=mean(ndvi, na.rm=T),
                  m.isa=mean(isa, na.rm=T),
                  m.barr=mean(barr, na.rm=T),
                  m.grass=mean(grass, na.rm=T),
                  count=.N,
                  m.visa=mean(vegisa, na.rm=T),
                  m.can=mean(can, na.rm=T),
                  m.visafrac=mean(visafrac, na.rm=T),
                  m.ed10=mean(ed10, na.rm=T),
                  m.ed20=mean(ed20, na.rm=T),
                  m.ed30=mean(ed30, na.rm=T),
                  m.buff.Intonly=mean(buff.Intonly, na.rm=T)), by=barr.bin]
hdres.m[,count.frac:=count/dim(hdres.p)[1]]
hdres.m[,pure.class:="hdres"]
plot(hdres.m$m.isa, hdres.m$m.evi)
plot(hdres.m$m.isa, hdres.m$m.barr)
plot(hdres.m$m.barr, hdres.m$m.evi)

lowveg.p <- dat[lowveg>0.95, ]
lowveg.p[,visafrac:=0]
lowveg.p[can>0,visafrac:=vegisa/can] ## fraction of canopy area that is over isa
lowveg.p[visafrac<=1,] ### there are weird and not that rare artifacts where e.g. very high grass pixels show higher VISA than canopy -- wtf? No visual anomaly on the rasters -- complate description of cover with ISA/NIB/CAN/VISA/GRASS
lowveg.m <- lowveg.p[,.(m.evi=mean(evi, na.rm=T), 
                      m.ndvi=mean(ndvi, na.rm=T),
                      m.isa=mean(isa, na.rm=T),
                      m.barr=mean(barr, na.rm=T),
                      m.grass=mean(grass, na.rm=T),
                      count=.N,
                      m.visa=mean(vegisa, na.rm=T),
                      m.can=mean(can, na.rm=T),
                      m.visafrac=mean(visafrac, na.rm=T),
                      m.ed10=mean(ed10, na.rm=T),
                      m.ed20=mean(ed20, na.rm=T),
                      m.ed30=mean(ed30, na.rm=T),
                      m.buff.Intonly=mean(buff.Intonly, na.rm=T)), by=barr.bin]
lowveg.m[,count.frac:=count/dim(lowveg.p)[1]]
lowveg.m[,pure.class:="lowveg"]
plot(lowveg.m$m.isa, lowveg.m$m.evi)
plot(lowveg.m$m.isa, lowveg.m$m.barr)
plot(lowveg.m$m.barr, lowveg.m$m.evi)

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

### combine into common data frame
puro <- rbind(forest.m, dev.m, hdres.m, lowveg.m)
puro[,pure.class:=(as.factor(pure.class))]
puro[m.visafrac>1, m.visafrac:=NA]

plot(puro[pure.class=="forest", m.isa], puro[pure.class=="forest", count.frac]) ### dominated by 0% ISA
plot(puro[pure.class=="lowveg", m.isa], puro[pure.class=="lowveg", count.frac]) ## dominated by 0% ISA
plot(puro[pure.class=="dev", m.isa], puro[pure.class=="dev", count.frac]) ### dominated by 100% ISA
plot(puro[pure.class=="hdres", m.isa], puro[pure.class=="hdres", count.frac]) ## some range of ISA, but still overrepresented in 100% ISA

### combined plots with all four pure-pixel classes
ggplot(puro, aes(x=m.barr, y=m.evi, color=pure.class, size=count.frac))+
  geom_point()+labs(title="EVI vs. Barren", colour="LULC class", size="Relative %", x="Barren", y="EVI")+
  # scale_size("% of LULC", range=c(0,8))+
  scale_size(range=c(1,10),breaks=c(0, 0.01, 0.02, 0.05, 0.10, 0.40),labels=c(">=0",">=1%",">=2%",">=5%",">=10%",">=40%"),guide="legend")## fine but point sizes are shitballs
ggplot(puro, aes(x=m.barr, y=m.can, color=pure.class, size=count.frac))+
  geom_point()+labs(title="%Canopy vs. Barrn", colour="LULC class", size="Relative %", x="Barren", y="Canopy")+
  # scale_size("% of LULC", range=c(0,8))+
  scale_size(range=c(1,10),breaks=c(0, 0.01, 0.02, 0.05, 0.10, 0.40),labels=c(">=0",">=1%",">=2%",">=5%",">=10%",">=40%"),guide="legend")## fine but point sizes are shitballs

#### I don't trust the veg-over-isa category, think the processing must be screwed up somehow.
ggplot(puro, aes(x=m.barr, y=m.visa, color=pure.class, size=count.frac))+
  geom_point()+labs(title="Canopy over ISA vs. ISA", colour="LULC class", size="Relative %", x="ISA", y="Vegetation over ISA")+
  # scale_size("% of LULC", range=c(0,8))+
  scale_size(range=c(1,10),breaks=c(0, 0.01, 0.02, 0.05, 0.10, 0.40),labels=c(">=0",">=1%",">=2%",">=5%",">=10%",">=40%"),guide="legend")## fine but point sizes are shitballs
ggplot(puro, aes(x=m.barr, y=m.visafrac, color=pure.class, size=count.frac))+
  geom_point()+labs(title="%Canopy over ISA vs. ISA", colour="LULC class", size="Relative %", x="ISA", y="Fraction of Vegetation over ISA")+
  # scale_size("% of LULC", range=c(0,8))+
  scale_size(range=c(1,10),breaks=c(0, 0.01, 0.02, 0.05, 0.10, 0.40),labels=c(">=0",">=1%",">=2%",">=5%",">=10%",">=40%"),guide="legend")## fine but point sizes are shitballs


### same analysis viz. barren fraction





### invividual pure-pixels plots
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
#### multiple plots
library(reshape2)
df <- data.frame(a = rnorm(10), c = rnorm(10), g = rnorm(10), class = sample(letters[20:23], 10, TRUE))
df.m <- melt(df)


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


