library(data.table)
library(reticulate)
library(raster)

hf.dbh <- read.csv("docs/ian/Harvard_Forest_biomass.csv") ## edge plto data for HF
length(unique(hf.dbh$Spp)) ## 18 species ID'd

andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv") ##edge plot data for Andy's suburban sites
length(unique(andy.dbh$Species_Code)) ## 17 species ID'd

andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv")
dim(andy.bai) ## 195 individuals, dbh by year -- convert to biomass by year? Edge plot data from surburb sites

rope <- read.csv("docs/ian/Rope_Plots_Data_Entry.csv") ## Rope plots, probably not useful, just a description of 

### should we also see if Britt's biomass increment~biomass is worth trying in here?

street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
length(unique(street$Species)) ## 77 species!
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
mean.na <- function(x){mean(x, na.rm=T)}
tapply(street$ann.npp, INDEX = street$Species, FUN = mean.na)
tapply(street$ann.npp, INDEX = street$Species, FUN=length)
# genspec <- strsplit(as.character(street$Species), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# spec <- unlist(lapply(genspec, "[[", 2))
# a <- regexpr("^[[:alpha:]]+", street$Species)
# b <- a+attr(a, "match.length")-1
# unique(gen)
# street$genus <- gen

### clean up taxonomies
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
a <- street[,length(ann.npp)/dim(street)[1], by=genus]
a <- a[order(a$V1, decreasing = T),]
sum(a$V1[1:12])
focus <- a$genus[1:12]
## dominated by Acer, Gleditsia, Tilia, Zlekova, Ulmus, Fraxinus
plotme <- street[standing.dead.2014==0 & genus %in% focus,]
plotme <- plotme[ann.npp>0,]
plot(street$dbh.2014, street$ann.npp, ylim=c(0, 110))
plot(plotme$dbh.2014, plotme$ann.npp)
plot(plotme$biomass.2006, plotme$ann.npp)
plot(plotme$biomass.2014, plotme$ann.npp)
biom.top <- plotme[,mean(biomass.2014)+1.96*(sd(biomass.2014))]
biom.bot <- plotme[,mean(biomass.2014)-1.96*(sd(biomass.2014))]
plot(plotme[biomass.2014<biom.top, biomass.2014], plotme[biomass.2014<biom.top, ann.npp])
plot(plotme[biomass.2014<biom.top, dbh.2014], plotme[biomass.2014<biom.top, ann.npp])

lin.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~plotme[biomass.2014<biom.top, biomass.2014])
summary(lin.mod)
pol.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~poly(plotme[biomass.2014<biom.top, biomass.2014], 2))
summary(pol.mod)


### biomass in 2014 doesn't really predict annual growth much in the preceding years, R2 0.43 for poly2
plot(plotme[biomass.2014<biom.top, biomass.2014], plotme[biomass.2014<biom.top, ann.npp])
plot(fitted(pol.mod), residuals(pol.mod))
q <- seq(0, 1500)
y2 <- pol.mod$coefficients[1]+(pol.mod$coefficients[2]*q)+(pol.mod$coefficients[3]*(q^2))
y <- lin.mod$coefficients[1]+(lin.mod$coefficients[2]*q)
lines(q, y)
plot(q, y2)

for(g in 1:length(focus)){
  # plot(plotme[genus==focus[g], dbh.2014], plotme[genus==focus[g], ann.npp], 
  #      main=focus[g])
  print(summary(lm(plotme[genus==focus[g], ann.npp]~plotme[genus==focus[g], dbh.2014])))
  
}
## everything grows faster at higher dbh, some are sloppy (R2 0.3), some are tighter (0.6-0.8)
## combined top genera have 0.65*dbh2014 for ann.prod, R2 0.46
plot(plotme$biomass.2014, plotme$ann.npp)

## mixed hardwood allometric
test <- 1:100
b0 <- -2.48
b1 <- 2.48
plot(test, exp(b0+(b1*log(test))))

### process biomass to 30m aggregate
### reran modified bos1m_agg to get summed 30m biomass
### need to think about what the raster values mean when summed to 30m grid from 1m

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m 
biom1m <- raster("data/dataverse_files/bostonbiomass_1m.tif") ## this is in MgC/ha
can1m <- raster("data/dataverse_files/bostoncanopy_1m.tif")
## conversion factors: the summed biomass/30m (900m2) cell 

(51240/1E04) ## a total of about 5 MgC/cell at max (convert MgC/ha to /m2 -- grid cell value is sum of the 1m2 MgC/ha values)
(51240/1E04)*2 ## a total of about 10 Mg.biomass per cell at max
51240*(1E-04)*(1E04/900) ## or about 57 MgC/ha max, in 900m2 cell footprint

## compare to 1m map of biomass (MgC/ha) ## max of 84 MgC/ha at 1m
84*(1E-04)*1000 ### in 1m grid, max of 8.4 kgC/m2
84*(1E-04)*1000*900 ### in 900m2 30m grid cell, 7560 kgC/cell = 7.6 MgC/cell
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[, npp.ss:=5.104+(bos.biom30m/1E04)*2*1000*0.029] ### based on linear model of street tree growth~
biom.dat[bos.biom30m==0, npp.ss:=0] ## manually correct for places w no trees
npp.ss.r <- raster(biom)
npp.ss.r <- setValues(npp.ss.r, biom.dat$npp.ss) ## ok, just a linear transform, 0-300 kg biomass/cell/yr, equiv to 1.66 MgC/ha/yr

## check: can we get figures same as reported in Raciti et al?
plot(biom1m) ## seem to get the same figure as Raciti, in MgC/ha
biom1m.dat <- as.data.table(as.data.frame(biom1m))
names(biom1m.dat) <- c("biom.MgC.ha")
biom1m.dat[,MgC.cell:=biom.MgC.ha/1E4]
biom1m.totC <- biom1m.dat[,sum(MgC.cell, na.rm=T)] # 71k MgC total storage, contrast 355k MgC in Raciti
biom1m.totar <- biom1m.dat[!is.na(MgC.cell), length(MgC.cell)]
biom1m.totar/ncell(biom1m) ### 37% of study area has valid pixels for biomass 
biom1m.totar.km <- biom1m.totar/1E6 ## 124 km2 (12 kha), contrast 123.3 km2 reported total
biom1m.totC/(biom1m.totar.km*100) ## this represents avg of 5.75 MgC/ha for the study area, 
28.8*biom1m.totar/10^4 ## would need 357k MgC total storage to make the figures reported in Raciti
#### Summary: Raciti reports greater overall C storage than is represented in his data, not clear where/how this deviates from other reported results


plot(can1m)
can1m.dat <- as.data.table(as.data.frame(can1m))
names(can1m.dat) <- c("can")
can.totcan <- can1m.dat[!is.na(can), sum(can)]
can.totar <- can1m.dat[!is.na(can), length(can)]
can.totcan/can.totar ## avg cover 25.8%, contrast 25.5% reported 
### summary: canopy data looks about like what was reported in Raciti
