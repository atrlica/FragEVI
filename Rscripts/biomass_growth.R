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
plot(plotme$biomass.2014, plotme$ann.npp) ### this is biomass kg per tree, related to biomass gain per year
biom.top <- plotme[,mean(biomass.2014)+1.96*(sd(biomass.2014))]
biom.bot <- plotme[,mean(biomass.2014)-1.96*(sd(biomass.2014))]
plot(plotme[biomass.2014<biom.top, biomass.2014], plotme[biomass.2014<biom.top, ann.npp])
plot(plotme[biomass.2014<biom.top, dbh.2014], plotme[biomass.2014<biom.top, ann.npp]) ## we don't have dbh

lin.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~plotme[biomass.2014<biom.top, biomass.2014])
summary(lin.mod)
pol.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~poly(plotme[biomass.2014<biom.top, biomass.2014], 2))
summary(pol.mod)

### biomass in 2014 doesn't really predict annual growth much in the preceding years, R2 0.43 for poly2

for(g in 1:length(focus)){
  plot(plotme[genus==focus[g], biomass.2014], plotme[genus==focus[g], ann.npp],
       main=focus[g])
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
plot(biom1m) ## Raciti lists this as MgC/ha, but it's really Mg-biomass/ha
biom1m.dat <- as.data.table(as.data.frame(biom1m))
names(biom1m.dat) <- c("biom.raw")
biom1m.dat[, MgC.ha.corr:=2*biom.raw]
biom1m.dat[, Mgbiom:=biom.raw/(10^4)]
biom1m.dat[, mean(Mgbiom, na.rm=T)]
mean.biom.raw <- biom1m.dat[, mean(biom.raw, na.rm=T)] ## 5.75 Mgbiom/ha
sum.biom.raw <- biom1m.dat[, sum(biom.raw, na.rm=T)]
area <- biom1m.dat[!is.na(biom.raw), length(biom.raw)] ## 124 km2
(340000-315000)*(4695000-4680000)/(10^6) ## size of raw raster approx: 3.75E8 = 375km2
mean.biom.raw*(area/10^4) ## 71412 Mgbiom/ha
sum.biom.raw/(10^4) ## 71.4k Mgbiom
(sum.biom.raw/(10^4))/(area/(10^4)) ## 5.75 Mg biom/ha

## biomass formula, basically just lidar canopy height linear transform h(m)-->bm(kg)
biomass=(2.1015*h)+0.8455
biomass=82 ## max biomass in this 1m raster
h=27 ### actual canopy height in the highest biomass pixels (brief visual search of a single area) is about 27m
(2.1015*h)+0.8455 ### a single pixel of 27m height gives biomass of 57kg (this is consistent with the biom1m raster value)
biomass=57
## at max biomass/pixel
biomass=82
(biomass-0.8455)/2.1015 ### this would imply a 39m canopy height

### test, single open-grown tree
## a tree of 14m diameter, max height ~19m, area = 153m, ~30 kg biomass per pixel average
153*30 ## about 4.5 Mg biomass
### actual extract
library(rgeos)
library(rgdal)
tree.foot <- readOGR(dsn="processed/boston/tree.footprint.shp", layer="tree.footprint")
crs(tree.foot)
crs(biom1m)
tree.foot <- spTransform(tree.foot, crs(biom1m))
### this is a open-grown tree, 14m diameter canopy roughly circular, max height is 19m
a <- extract(biom1m, tree.foot)
a <- unlist(a)
sum(a) ### 3941 kg biomass in this tree
## general equation for biomass (Ian street tree)
biomass=exp(-2.48+(2.4835(log(dbh))))
biomass=sum(a)
exp((log(biomass)+2.48)/2.4835) ## this would imply a 76cm dbh tree with 19m height
## at 19m high, this tree, according to the transform in Raciti
h=19
(2.1015*h)+0.8455 ## about 40 kg/pixel
40*153 ## roughly correct, somewhat too high (max pixel height)

### CONCLUSION: RACITI'S 1m BIOMASS DATA is correct transform of lidar height, yielding kg-biomass/1m.cell
sum.biom.raw <- biom1m.dat[, sum(biom.raw, na.rm=T)]
area <- biom1m.dat[!is.na(biom.raw), length(biom.raw)] ## 124 km2
(sum.biom.raw/10^3)/(area/(10^4)) ### average of 57.5 Mgbiomass/ha
(sum.biom.raw/10^3)/(area/(10^4))/2 ### average of 28.8 MgC/ha --> same as raciti
(sum.biom.raw/10^3)/2/(10^3) ### 357 GgC total storage --> same as raciti
### the scale on the corrected figures from Raciti is what you get when you take
### kgbiomass/m2 * 1E4m2/ha * 1kgC/2kgbiomass * 1MgC/1E3kgC (e.g. 82 kgbiomass/m2 --> 410 MgC/ha)

plot(can1m)
can1m.dat <- as.data.table(as.data.frame(can1m))
names(can1m.dat) <- c("can")
can.totcan <- can1m.dat[!is.na(can), sum(can)]
can.totar <- can1m.dat[!is.na(can), length(can)]
can.totcan/can.totar ## avg cover 25.8%, contrast 25.5% reported 
### summary: canopy data looks about like what was reported in Raciti



#### OK with that out of the way: Let's try to figure out the tree distribution per 30m cell and uptake
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m 
plot(biom)
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
length(unique(street$Species)) ## 77 species!
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
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
focus

### approach 0: get npp biomass increment as function of median tree biomass, FAI growth rates
## take median value of tree dbh, ignore genus (using inline biomass ~ common hardwood allometric)
dbh <- street[standing.dead.2014==0 & genus %in% focus, median(dbh.2014, na.rm=T)]
biomass <- street[standing.dead.2014==0 & genus %in% focus, median(biomass.2014, na.rm=T)]
hist(street[standing.dead.2014==0 & genus %in% focus, biomass.2014]) ## highly skewed towards <1000kg
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,num.trees.med:=bos.biom30m/biomass]
biom.dat[,range(num.trees.med, na.rm=T)] ## max 280 trees per cell i.e. a tree per every 3m2

## the FIA use this equation for volume of wood in a stand by stand age
a=123.67
b=0.04
AGE=40
y = a*(1-exp(-b*AGE))^3
y
### standing live wood C storage in local forest that resemble (?) what we'd have in Boston:
### i.e. N.red oak; red maple/oak; mixed upland hwoods; red maple uplands: Range is  94.7-105.1 MgC/ha
## read in summaries of C stock change in Q.alba/Q.rubra/Carya and Q.rubra forest
## both forest types include values for reforestation and afforestation conditions
ne.summ <- read.csv("docs/FIA_CTMANHRI.csv")
plot(ne.summ$Age.yrs, ne.summ$live.tree.c.inc)

mix.oak.ref <- read.csv("docs/FIA_QUALQURUCATO_Reforest.csv")
plot(mix.oak.ref$Age.yrs, mix.oak.ref$live.tree.c.inc)
mix.oak.aff <- read.csv("docs/FIA_QUALQURUCATO_Afforest.csv")
plot(mix.oak.aff$Age.yrs, mix.oak.aff$live.tree.c.inc)

red.oak.ref <- read.csv("docs/FIA_QURU_Reforest.csv")
plot(red.oak.ref$Age.yrs, red.oak.ref$live.tree.c.inc)
red.oak.aff <- read.csv("docs/FIA_QURU_Afforest.csv")
plot(red.oak.aff$Age.yrs, red.oak.aff$live.tree.c.inc)

fia.summ <- as.data.frame(cbind(ne.summ$Age.yrs, ne.summ$live.tree.c.inc, mix.oak.ref$live.tree.c.inc,
                  mix.oak.aff$live.tree.c.inc, red.oak.ref$live.tree.c.inc, red.oak.aff$live.tree.c.inc))
colnames(fia.summ) <- c("Age", "NE.total", "mix.oak.ref", "mix.oak.aff", "red.oak.ref", "red.oak.aff")
plot(fia.summ$Age, fia.summ$NE.total, pch=15, col="black")
points(fia.summ$Age, fia.summ$mix.oak.aff, pch=16, col="lightblue")
points(fia.summ$Age, fia.summ$mix.oak.ref, pch=16, col="blue")
points(fia.summ$Age, fia.summ$red.oak.aff, pch=17, col="pink")
points(fia.summ$Age, fia.summ$red.oak.ref, pch=17, col="red")
### the only thing that changes between reforestation and afforestation are values for forest floor and soil C

