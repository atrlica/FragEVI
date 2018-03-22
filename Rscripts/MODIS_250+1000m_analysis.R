library(raster)
library(data.table)
library(rgeos)
library(rgdal)
library(sp)
library(ggplot2)
library(knitr)


######
###### MODIS 250m + 1 km analysis
### 250m analysis
resolution <- 1000
if(resolution==250){dat <- read.csv("processed/AOI.250m.dat.csv")}
if(resolution==1000){dat <- read.csv("processed/AOI.1km.dat.csv")}
dat <- dat[,2:12]
names(dat) <-  c("evi", "isa", "dev", "forest", "hdres", "lowres", "lowveg", "medres", "other", "water", "AOI")
dat <- as.data.table(dat)
dat <- dat[AOI==1,] 
dim(dat)[1] ## number of pixels available (106k for 250m, 7k for 1 km)
dat[,frac.tot:=apply(dat[,3:10], 1, sum)] ## most pixels are within +/- 5% of complete area with the LULC category fractions
dat[frac.tot<0.95, length(frac.tot)] ## only 500 pixels are missing more than 5% area fraction
dat[,evi:=evi/10000]

## get values for vegetation and non-veg endmembers
crithi <- 0.99
critlow <- 0.01
if(resolution==1000){crithi <- 0.92} ## correct for 1km where no 100% ISA available
veg <- dat[isa<critlow  & water<0.02, median(evi, na.rm=T)] ## 10k pixels are 0% ISA
veg.err <- dat[isa<critlow  & water<0.02, sd(evi, na.rm=T)/sqrt(length(evi))]
veg.n <- dat[isa<critlow  & water<0.02, length(evi)]
noveg <- dat[isa>=crithi & water<0.02, median(evi, na.rm=T)]
noveg <- dat[isa>=crithi & water<0.02, median(evi, na.rm=T)]
noveg.err <- dat[isa>=crithi & water<0.02, sd(evi, na.rm=T)/sqrt(length(evi))]
noveg.n <- dat[isa>=crithi & water<0.02, length(evi)] ### only 20 pixels are 100% ISA
beta.range <- seq(from=0, to=1, by=0.01) 
Vzi <- ((1-beta.range)*veg)+(beta.range*noveg)

## bin the values by isa fraction
dat[water<0.02, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02, .(m.isa=median(isa, na.rm=T), 
                             m.evi=median(evi, na.rm=T),
                             m.forest=median(forest, na.rm=T),
                             m.dev=median(dev, na.rm=T),
                             m.hdres=median(hdres, na.rm=T),
                             m.medres=median(medres, na.rm=T),
                             m.lowres=median(lowres, na.rm=T),
                             m.lowveg=median(lowveg, na.rm=T),
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
  scale_size_continuous(range=c(1, 6), breaks=c(0.001, 0.01, 0.02, 0.03, 0.04), name="% total")+
  labs(title=paste("EVI vs. binned ISA, ", resolution, "m, AOI", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="red")+
  labs(title=paste("EVI vs. binned ISA, ", resolution, "m, AOI", sep=""), x="Frac. ISA", y="EVI enhancement")+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)
## @250m 82k pixels below 50% ISA, 6k above 50%
dim(dat[isa<0.25,])[1]/dim(dat)[1] ## at 250m, 94% <50% ISA, at 1km, 95% ie. 6.7k pixels are <50% ISA

### stack plot of fractional area by ISA bin
dat[water<0.02, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin2 <- dat[water<0.02, .(m.isa=median(isa, na.rm=T), 
                              m.evi=median(evi, na.rm=T),
                              s.forest=sum(forest, na.rm=T),
                              s.dev=sum(dev, na.rm=T),
                              s.hdres=sum(hdres, na.rm=T),
                              s.medres=sum(medres, na.rm=T),
                              s.lowres=sum(lowres, na.rm=T),
                              s.lowveg=sum(lowveg, na.rm=T),
                              bin.count=.N), by=bin]
evi.bin2 <- evi.bin2[order(m.isa),]
tot <- evi.bin2[,sum(bin.count)]
evi.bin2[,bin.frac:=bin.count/tot, by=bin]
evi.bin2 <- evi.bin2[!is.na(m.isa),]

### total area present across ISA range
d <- reshape(evi.bin2, varying=c("s.forest", "s.lowveg", "s.lowres", "s.medres", "s.hdres", "s.dev"),
             idvar="bin", direction="long")
d$time <- as.factor(d$time)
ggplot(d, aes(x=m.isa, y=s/16))+
  geom_area(aes(fill=time))+
  scale_fill_manual(values=c("gray45", "forestgreen", "salmon", "gold", "darkolivegreen3", "darkorange"),
                    breaks=c("forest", "lowveg", "lowres", "medres", "hdres", "dev"), name="LULC",
                    labels=c("Forest", "Low Veg.", "LD+VLD Resid.", "MD Resid.", "HD+MF Resid.", "Developed"))+
  labs(x="Frac. Impervious", y="Total area (km2)", title=paste("LULC composition vs ISA, ", resolution, "m", sep=""))

#### relative area present across ISA range
evi.bin2[,tot.area:=.(s.forest+s.dev+s.hdres+s.medres+s.lowres+s.lowveg), by=bin]
e <- reshape(evi.bin2, varying=c("s.forest", "s.lowveg", "s.lowres", "s.medres", "s.hdres", "s.dev"),
             idvar="bin", direction="long")
e$time <- as.factor(e$time)
e[,rel.area:=s/tot.area]
ggplot(e, aes(x=m.isa, y=rel.area))+
  geom_area(aes(fill=time))+
  scale_fill_manual(values=c("gray45", "forestgreen", "salmon", "gold", "darkolivegreen3", "darkorange"),
                    breaks=c("forest", "lowveg", "lowres", "medres", "hdres", "dev"), name="LULC",
                    labels=c("Forest", "Low Veg.", "LD+VLD Resid.", "MD Resid.", "HD+MF Resid.", "Developed"))+
  labs(x="Frac. Impervious", y="Frac. of AOI", title=paste("LULC composition vs ISA, ", resolution, "m", sep=""))


### plots of EVI vs. fractional LULC
ggplot(aes(x=m.forest, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="forestgreen")+
  labs(title=paste("EVI vs. forest,", resolution, "m, AOI", sep=""), x="Frac. forest", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## saturates EVI above 75% forest, 0 forest ranges from 0.1-0.3

ggplot(aes(x=m.dev, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="gray45")+
  labs(title=paste("EVI vs. Dev,", resolution, "m, AOI", sep=""), x="Frac. Dev", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## sloppy, everything above 50% dev is <0.2

ggplot(aes(x=m.hdres, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="salmon")+
  labs(title=paste("EVI vs. Hdres,", resolution, "m, AOI", sep=""), x="Frac. HD Resid.", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## even at max HD resid (>60%) EVI is still ~0.3

ggplot(aes(x=m.medres, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="orange")+
  labs(title=paste("EVI vs. MDres,", resolution, "m, AOI", sep=""), x="Frac. MD Resid.", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## medres never gets very common; the 0% includes very urban and very not urban at either end of scale, but only a slight decrease with increasing MDres

ggplot(aes(x=m.lowres, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="goldenrod")+
  labs(title=paste("EVI vs. LDres,", resolution, "m, AOI", sep=""), x="Frac. LD Resid.", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## similar to MD resid -- even a bit rarer, but little effect on dragging down EVI

ggplot(aes(x=m.lowveg, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="lightgreen")+
  labs(title=paste("EVI vs. lowveg,", resolution, "m, AOI", sep=""), x="Frac. low veg", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
## stays rare (<12%), but weakly positive effect



## do EVI vs ISA coloring the points according to majority area
lulc.list <- c("dev", "forest", "hdres", "lowres", "medres", "other")
evi.bin[,lulc.maj:=NA]
for(m in evi.bin[order(bin), bin]){
  tot.m <- dat[bin==m, .N]
  uu <- unlist(dat[bin==m, .(median(dev), median(forest),median(hdres),median(lowres),median(medres),median(other))])
  # if(max(uu)>0.5) ## this just gets the plurality -- not majority necessarily
  evi.bin[bin==m, lulc.maj:=lulc.list[which(uu==max(uu))]]
}
evi.bin <- evi.bin[order(bin),]
evi.bin$lulc.maj <- as.factor(evi.bin$lulc.maj)


#EVI vs. ISA plot, color by majority area fraction
# not valid at 1km -- no pixels emerge
ggplot(aes(x=m.isa, y=m.evi, color=lulc.maj), data=evi.bin)+
  geom_point(aes(size=bin.frac))+
  scale_color_manual(breaks=c("1", "2", "3"),
                     values=c("forestgreen", "gray45", "salmon"),
                     guide=guide_legend(title="largest LULC"),
                     labels=c("Forest", "Dev", "HD Resid"))+
  labs(title=paste("EVI vs. binned ISA, ", resolution, "m, AOI", sep=""), x="Frac. ISA", y="Median EVI")+
  # scale_size(range=c(0.7,10),breaks=c(0.01, 0.02, 0.40),
  #            labels=c(">1%","1-2%", ">40%"),
  #            guide=guide_legend(title="Bin size"))+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100], color=lulc.maj), data=evi.bin)+
  geom_point(aes(size=bin.frac))+
  scale_color_manual(breaks=c("1", "2", "3"),
                     values=c("forestgreen", "gray45", "salmon"),
                     guide=guide_legend(title="largest LULC"),
                     labels=c("Forest", "Dev", "HD Resid"))+
  labs(title=paste("EVI enhancement, ", resolution, "m, AOI", sep=""), x="Frac. ISA", y="EVI enhancement")+
  # scale_size(range=c(0.7,10),
  #            breaks=c(0.01, 0.02, 0.40),
  #            labels=c(">1%","1-2%", ">40%"),
  #            guide=guide_legend(title="Bin size"))+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)

### EVI enhancement in pure pixel samples
### developed
dat[water<0.02 & dev>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & dev>0.5, .(m.isa=median(isa, na.rm=T), 
                                       m.evi=median(evi, na.rm=T),
                                       bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #5k pixels (143 pix at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="gray45")+
  scale_size(name="% total")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% Developed", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="blue")+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% Developed", sep=""), x="Frac. ISA", y="EVI enhancement")+
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)

## hdres
dat[water<0.02 & hdres>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & hdres>0.5, .(m.isa=median(isa, na.rm=T), 
                                         m.evi=median(evi, na.rm=T),
                                         bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #6k pixels (243 pix at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="salmon")+
  scale_size(name="% total")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% HD+MF Resid.", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
# not run: not all ISA bins represented
# ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
#   geom_point(aes(size=bin.frac), colour="blue")+
#   labs(title="EVI enhancement, >50% HD Resid", x="Frac. ISA", y="EVI enhancement")+
#   geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)

## medres
dat[water<0.02 & medres>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & medres>0.5, .(m.isa=median(isa, na.rm=T), 
                                          m.evi=median(evi, na.rm=T),
                                          bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #5k pixels (141 pix at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="orange")+
  scale_size(name="% total")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% MD Resid.", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
# ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
#   geom_point(aes(size=bin.frac), colour="blue")+
#   labs(title="EVI enhancement, >50% HD Resid", x="Frac. ISA", y="EVI enhancement")+
#   geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)

## lowres
dat[water<0.02 & lowres>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & lowres>0.5, .(m.isa=median(isa, na.rm=T), 
                                          m.evi=median(evi, na.rm=T),
                                          bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #4k pixels (47 pix at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="gold")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% VLD+LD Resid.", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)
# ggplot(aes(x=m.isa, y=m.evi-Vzi[1:100]), data=evi.bin)+
#   geom_point(aes(size=bin.frac), colour="blue")+
#   labs(title="EVI enhancement, >50% HD Resid", x="Frac. ISA", y="EVI enhancement")+
#   geom_hline(yintercept=0, linetype="solid", color="black", size=0.7)

## other 
dat[water<0.02 & other>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & other>0.5, .(m.isa=median(isa, na.rm=T), 
                                         m.evi=median(evi, na.rm=T),
                                         bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #199 pixels (5 pix at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="red")+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% other", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)


## forest
dat[water<0.02 & forest>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & forest>0.5, .(m.isa=median(isa, na.rm=T), 
                                          m.evi=median(evi, na.rm=T),
                                          bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #44k pixels (2.5k at 1km)

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="forestgreen")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% Forest", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)


## lowveg
dat[water<0.02 & lowveg>0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & lowveg>0.5, .(m.isa=median(isa, na.rm=T), 
                                          m.evi=median(evi, na.rm=T),
                                          bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #6k pixels, 229 pix at 1km

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="darkolivegreen3")+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(noveg, veg))+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, >50% low veg", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2)


## mixed
dat[water<0.02 & dev<0.5 & forest<0.5 & hdres<0.5 & lowres<0.5 & lowveg<0.5 & medres<0.5 & other<0.5, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[water<0.02 & dev<0.5 & forest<0.5 & hdres<0.5 & lowres<0.5 & lowveg<0.5 & medres<0.5 & other<0.5, .(m.isa=median(isa, na.rm=T), 
                                                                                                                   m.evi=median(evi, na.rm=T),
                                                                                                                   bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]
evi.bin[,sum(bin.count)] #17k pixels, 1.6k at 1km

ggplot(aes(x=m.isa, y=m.evi), data=evi.bin)+
  geom_point(aes(size=bin.frac), colour="orchid")+
  labs(title=paste("EVI vs. ISA, ", resolution, "m, mixed pixels", sep=""), x="Frac. ISA", y="Median EVI")+
  geom_abline(intercept=veg, slope=(noveg-veg), color="red", linetype="dashed", size=1.2) ## mixed pixels look like they tend to be less paved

