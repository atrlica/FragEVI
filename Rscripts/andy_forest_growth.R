library(raster)
library(data.table)

######
###### Treat all canopy as Andy-like forest

# ### procecssing 1m edge biomass to 30m cells
# ## identify 1m biomass pixels that are edge canopy (all LULC types)
# ed <- raster("processed/boston/bos.ed10.tif")
# biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# aoi <- raster("processed/boston/bos.aoi.tif")
# # biom <- projectRaster(biom, aoi)
# biom <- crop(biom, aoi)
# 
# ### identify biomass that is edge or not (all LULC included)
# edges.biom <- function(x, y, filename) { # x is edge class, y is biomass, 
#   out <- raster(y)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     e <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## edge class
#     b <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) ## biomass
#     b[e==0] <- 0 ## cancel non-edges
#     out <- writeValues(out, b, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# ed.biom <- edges.biom(ed, biom, filename="processed/boston/bos.biomass.ed10only_allLULC.tif")
# r <- raster("processed/boston/bos.biomass.ed10only_allLULC.tif")
# plot(r)

### convert to 1m edge biomass map to 30m aggregate (arcpy)

### now working in the 30m space
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
ed.can <- raster("processed/boston/bos.ed1030m.tif")
can <- raster("processed/boston/bos.can30m.tif")
ed.biom <- raster("processed/boston/bos.biomass.ed10only_allLULC30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")

biom.dat <- as.data.table(as.data.frame(biom))
aoi.dat <- as.data.table(as.data.frame(aoi))
ed.can.dat <- as.data.table(as.data.frame(ed.can))
can.dat <- as.data.table(as.data.frame(can))
ed.biom.dat <- as.data.table(as.data.frame(ed.biom))
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat <- cbind(biom.dat, aoi.dat, ed.can.dat, can.dat, ed.biom.dat, isa.dat)
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
names(biom.dat) <- c("biom", "aoi", "ed.can", "can", "ed.biom", "isa", "pix.ID")

## figure out the relative area of interior forestcanopy in each cell
biom.dat[,int.can:=can-ed.can]
biom.dat[,range(int.can, na.rm=T)] ## includes negative numbers
biom.dat[int.can<0, length(int.can)] ## 112, some are rounding errors
## kill the artifacts
biom.dat[biom==0, int.can:=0]
biom.dat[biom==0, ed.can:=0]
biom.dat[is.na(biom), ed.can:=NA]
biom.dat[is.na(biom), int.can:=NA]
biom.dat[int.can<0 & int.can>(-0.01), int.can:=0] ## kill rounding errors
biom.dat[,int.can:=can-ed.can]
# biom.dat[int.can<0, length(int.can)] ## 31
# View(biom.dat[int.can<0,]) ## all areas with 0 or partial forest, places with forest biom>0 are 100% edge biomass
# biom.dat[ed.biom==forest.biom & int.can>0,]
# biom.dat[ed.can>can, ] ## almost all partial pixels, usually forest edge is 100% of forest biomass, seem like minor disagreements in canopy area figuring
biom.dat[ed.can>can, int.can:=0]
# biom.dat[ed.can>can,] ## still 31 records where ed.can doesn't make sense, but fuck it
# biom.dat[ed.can==can,] ## 113k records where canopy is all edge canopy (vast bulk of the pixels)
# biom.dat[aoi>500,] ##138k mostly complete cells

## now ID biomass by edge vs int
biom.dat[,int.biom:=biom-ed.biom] # internal forest biomass
# biom.dat[,range(int.biom, na.rm=T)] ## 0-51k 
# biom.dat[,range(ed.biom, na.rm=T)] ## 0-31k (so the peak biomass cells are deep forest somewhere)
hist(biom.dat[aoi>800, int.biom])
hist(biom.dat[aoi>800, ed.biom])


##### ANALYSIS OF TREE CORE RECORDS
## 
andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv") 
andy.bai <- as.data.table(andy.bai)
### the basic allometrics to get biomass
## Jenkins, C.J., D.C. Chojnacky, L.S. Heath and R.A. Birdsey. 2003. Forest Sci 49(1):12-35.
## Chojnacky, D.C., L.S. Heath and J.C. Jenkins. 2014. Forestry 87: 129-151.
# Pinus rigida, Pinus strobus, both spg>0.45
# Acer rubrum, Aceraceae <0.50 spg
# Quercus alba, Quercus coccinea, Quercus rubra, Quercus velutina --> deciduous Fagaceae
b0.l <- c(-3.0506, -2.0470, -2.0705)
b1.l <- c(2.6465, 2.3852, 2.4410)
biom.pred <- function(x, b0, b1){exp(b0+(b1*log(x)))} ## dbh in cm, biom. in kg
## artifact: dbh's have hidden NA's that are marked as repeating numbers, basically when the dbh gets too small -- they are higher than the min dbh though
cleanup <- as.matrix(andy.bai[,7:33])
for(r in 1:nrow(cleanup)){
  bust <- which(diff(cleanup[r,], lag = 1)>0) ## tells you where the NA's are located
  if(length(bust)>0){
    cleanup[r, which(colnames(cleanup)==names(bust)):ncol(cleanup)] <- NA
  }
}
cleanup <- as.data.table(cleanup)
cleanup$num.yrs <- apply(cleanup, 1, function(x){return(sum(!is.na(x)))})
andy.bai <- cbind(andy.bai[,1:6], cleanup)
names(andy.bai)[7:33] <- paste0("dbh", 2016:1990)
andy.bai[,incr.ID:=seq(1:dim(andy.bai)[1])] ## Tree.ID is not unique to each core
names(andy.bai)[1:6] <- c("Plot.ID", "Tree.ID", "Spp", "X", "Y", "Can.class")

### biomass change as a % of previous biomass per year from the increment data
taxa <- list(c("PIRI", "PIST"),
             c("ACRU"),
             c("QUAL", "QUCO", "QURU", "QUVE"))
bop <- data.table()
## figure out biomass growth history applying species-specific allometrics
for(sp in 1:3){
  tmp <- andy.bai[Spp %in% taxa[[sp]],]
  a <- tmp[, lapply(tmp[,7:33], function(x){biom.pred(x, b0.l[sp], b1.l[sp])})]
  a <- cbind(tmp[,.(Tree.ID, incr.ID, Spp, X, Y)], a)
  bop <- rbind(bop, a)
}
bop$Spp <- as.character(bop$Spp) ### bop is the biomass time series for every core
names(bop)[6:32] <- paste0("biom", 2016:1990)

### lagged biomass gains, relative to previous year biomass
biom.rel <- data.frame(c(2015:1990))
for(i in 1:dim(bop)[1]){
  te <- unlist(bop[i,6:32])
  te.di <- (-1*diff(te, lag=1))
  biom.rel[,i+1] <- te.di/te[-1]
}
names(biom.rel)[-1] <- paste0("Tree", bop[,incr.ID])
names(biom.rel)[1] <- "Year" ## this has all the relative forward biomass gains by year (row) for each tree with increment (column)

## exploratory plots of tree growth rates through time
par(mfrow=c(2,2), mar=c(1,2,2,1), oma=c(2,1,1,1))
plot(biom.rel$Year, biom.rel$Tree2)
plot(biom.rel$Year, biom.rel$Tree4)
plot(biom.rel$Year, biom.rel$Tree14)
plot(biom.rel$Year, biom.rel$Tree188)
# for(u in 41:60){  
#   plot(biom.rel$Year, biom.rel[,u], main=paste(colnames(biom.rel)[u]))
# }
## everything looks like a slow decline in relative growth (not clear if that is just a funciton of trees getting bigger or canopy closing)
## most are in the 2-8% range, but some are v. high (20-40%)

### initial approach: take mean of rel. biomass gain and mean of dbh across all the years present in each core sample
### figure the mean relative biomass gain per year across all years in each tree
mean.na <- function(x){mean(x, na.rm=T)}
m <- apply(biom.rel[,-1], FUN=mean.na, 2)
growth.mean <- data.frame(cbind(c(sub("Tree", '', colnames(biom.rel[,-1]))), m))
names(growth.mean) <- c("incr.ID", "growth.mean")
growth.mean$incr.ID <- as.integer(as.character(growth.mean$incr.ID))
andy.bai <- merge(andy.bai, growth.mean, by="incr.ID")
andy.bai$growth.mean <- as.numeric(as.character(andy.bai$growth.mean))
andy.bai[,avg.dbh:=apply(as.matrix(andy.bai[,8:34]), FUN=mean.na, 1)]
## note: final interaction model only applies correction coefficients for the interior segments

### add guide data on tree position viz. edge
andy.bai[Y<10, seg:=10]
andy.bai[Y>=10 & Y<20, seg:=20]
andy.bai[Y>=20, seg:=30]
andy.bai[seg==10, seg.F:="A"]
andy.bai[seg==20, seg.F:="B"]
andy.bai[seg==30, seg.F:="C"]
andy.bai$seg.F <- as.factor(andy.bai$seg.F)
andy.bai[seg.F=="A", seg.Edge:="E"]
andy.bai[seg.F %in% c("B", "C"), seg.Edge:="I"]
andy.bai[,seg.Edge:=as.factor(seg.Edge)]
write.csv(andy.bai, "processed/andy.bai.dbh.avg.results.csv")

####
### upgrade Jul 23: treat five-year time slices as pseudo-replicates
pseudo.A <- matrix(ncol=4, rep(999,4))
pseudo.B <- matrix(ncol=4, rep(999,4))
pseudo.C <- matrix(ncol=4, rep(999,4))
pseudo.D <- matrix(ncol=4, rep(999,4))
pseudo.E <- matrix(ncol=4, rep(999,4))

### for each tree core, pick successive 5-year chunks moving backward until you run out of rings
for(d in unique(andy.bai$incr.ID)){
  pseudo.A <- rbind(pseudo.A, c(((bop[incr.ID==d, biom2016]-bop[incr.ID==d, biom2011])/bop[incr.ID==d, biom2011])/5,
                                andy.bai[incr.ID==d, dbh2011],
                                andy.bai[incr.ID==d, incr.ID],
                                "A"))
  pseudo.B <- rbind(pseudo.B, c(((bop[incr.ID==d, biom2011]-bop[incr.ID==d, biom2006])/bop[incr.ID==d, biom2006])/5,
                                andy.bai[incr.ID==d, dbh2006],
                                andy.bai[incr.ID==d, incr.ID],
                                "B"))
  pseudo.C <- rbind(pseudo.C, c(((bop[incr.ID==d, biom2006]-bop[incr.ID==d, biom2001])/bop[incr.ID==d, biom2001])/5,
                                andy.bai[incr.ID==d, dbh2001],
                                andy.bai[incr.ID==d, incr.ID],
                                "C"))
  pseudo.D <- rbind(pseudo.D, c(((bop[incr.ID==d, biom2001]-bop[incr.ID==d, biom1996])/bop[incr.ID==d, biom1996])/5,
                                andy.bai[incr.ID==d, dbh1996],
                                andy.bai[incr.ID==d, incr.ID],
                                "D"))

  ### pseudo.E can have from 5-6 years in a successsful sample
  tmp <- unlist(bop[incr.ID==d, 26:32, with=F])
  t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
  if(t<6 | !is.finite(t)){ ## if record doesn't end in 1991 or 1990 (i.e. at least 5 year run)
    psuedo.E <- rbind(pseudo.E, c(999,999,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    # print("YOU LIKE IS FUNKY") ## this record not long enough
  }else{
    g <- ((bop[incr.ID==d, biom1996]-bop[incr.ID==d, (26+(t-1)), with=F])/bop[incr.ID==d, (26+(t-1)), with=F])/(sum(is.finite(tmp))-1)
    deeb <- min(unlist(andy.bai[incr.ID==d, 28:34, with=F]), na.rm=T)
    pseudo.E <- rbind(pseudo.E, c(g, deeb,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    # print("I LIKE IS CHUNKY")
  }
}
## collate
ps.contain <- rbind(pseudo.A[-1,], pseudo.B[-1,], pseudo.C[-1,], pseudo.D[-1,], pseudo.E[-1,])
ps.contain <- data.frame(biom.rel.ann=as.numeric(unlist(ps.contain[,1])),
                 dbh.start=as.numeric(unlist(ps.contain[,2])),
                 incr.ID=as.integer(unlist(ps.contain[,3])),
                 interval=as.character(unlist(ps.contain[,4])),
                 stringsAsFactors = F)
ps.contain$interval <- as.factor(ps.contain$interval)
ps.contain <- as.data.table(ps.contain)

ps.contain[is.finite(biom.rel.ann) & dbh.start>5,] ## 900 reasonably big trees 
plot(ps.contain$dbh.start, ps.contain$biom.rel.ann, col=ps.contain$interval, ylim=c(0,1))
ps.contain <- merge(x=ps.contain, y=andy.bai[,.(incr.ID, seg, seg.Edge)], by="incr.ID", all.x=T, all.y=F)
plot(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, dbh.start], 
     ps.contain[is.finite(biom.rel.ann) & dbh.start>5, biom.rel.ann],
     col=ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge])

### dump the pseudo-replicated file to disk
write.csv(ps.contain, "processed/andy.bai.dbh.pseudo.csv")

### anlaysis of growth rate~dbh
par(mfrow=c(1,1))

# ## avg.dbh~avg.growth
# plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean), 
#      col=as.numeric(andy.bai$seg.Edge), main="Andy trees, avg. dbh")
# summary(andy.bai$growth.mean)
# table(andy.bai$seg.Edge) ## 64 edge: 131 interior
# ## vs. pseudoreps
# plot(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, log(dbh.start)], 
#      ps.contain[is.finite(biom.rel.ann) & dbh.start>5, log(biom.rel.ann)],
#      col=ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge],
#      main="Andy trees, dbh pseudoreps")
# table(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge]) ## 268 edge:632 interior
# 
# summary(andy.bai$growth.mean)
# summary(ps.contain$biom.rel.ann) ## comparable in the middle, pseudos have greater range
# plot(ps.contain[incr.ID==1, dbh.start], ps.contain[incr.ID==1, biom.rel.ann])
# plot(ps.contain[incr.ID==4, dbh.start], ps.contain[incr.ID==4, biom.rel.ann])
# plot(ps.contain[incr.ID==11, dbh.start], ps.contain[incr.ID==11, biom.rel.ann])
# plot(ps.contain[incr.ID==100, dbh.start], ps.contain[incr.ID==100, biom.rel.ann]) ## Ok, comforting. Big trees grow slowly, and within tree growth declines with time/dbh

### modeling growth~dbh*edge
### ifull interactive, avg.dbh and pseudo
mod.int.edge <- summary(lm(log(growth.mean)~log(avg.dbh)*seg.Edge, data=andy.bai)) #r2=0.47, no sig. dbh*edge interaction
mod.int.edge.ps <- summary(lm(log(biom.rel.ann)~log(dbh.start)*seg.Edge, data=ps.contain[dbh.start>=5,])) #r2=0.33, no sig. dbh*edge interaction

## removing the slope modifier (i.e. just a different intercept for edge vs. interior)
mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai)) #r2=0.47, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]

mod.int.edge.fin.ps <- summary(lm(log(biom.rel.ann)~log(dbh.start)+seg.Edge, data=ps.contain[dbh.start>=5,])) #r2=0.32, just edge vs. interior
b0.bai.ps <- mod.int.edge.fin.ps$coefficients[1]
b1.bai.ps <- mod.int.edge.fin.ps$coefficients[2]
b2.bai.ps <- mod.int.edge.fin.ps$coefficients[3]
### so really pretty comparable results at individual stem growth level looking at psuedoreps or avg.dbh/growth

### generalized growth coefficients for edge/interior, based on andy.bai
par(mfrow=c(1,2))
col.edge <- c("black", "red")
plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean), col=col.edge[as.numeric(as.factor(andy.bai$seg.Edge))],
     pch=15, cex=0.5, main="Andy Trees, avg. dbh")
abline(b0.bai, b1.bai, col="black")
abline(b0.bai+b2.bai, b1.bai, col="red")

plot(log(ps.contain$dbh.start), log(ps.contain$biom.rel.ann), col=col.edge[as.numeric(as.factor(ps.contain$seg.Edge))],
     pch=15, cex=0.5, main="Andy Trees, Pseudoreps")
abline(b0.bai.ps, b1.bai.ps, col="black")
abline(b0.bai.ps+b2.bai.ps, b1.bai.ps, col="red")
# plot((ps.contain$dbh.start), (ps.contain$biom.rel.ann), col=col.edge[as.numeric(as.factor(ps.contain$seg.Edge))],
#      pch=15, cex=0.5, main="Andy Trees, Pseudoreps")

# ## non-transformed space ### average growth figures
# mod.int.edge.fin.nls <- nls(growth.mean ~ exp(a + b * log(avg.dbh)), data=andy.bai, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
# summary(mod.int.edge.fin.nls)
# col.edge <- c("royalblue", "purple")
# plot(andy.bai$avg.dbh, andy.bai$growth.mean, col=col.edge[as.numeric(as.factor(andy.bai$seg.Edge))],
#      pch=15, cex=0.3, ylab="Growth Rate (kg/kg)", xlab="Stem DBH (cm)", main="Urban Forest stems")
# points(andy.bai$avg.dbh, exp(b0.bai)*exp(b1.bai*log(andy.bai$avg.dbh)), 
#        col=col.edge[1], pch=13, cex=0.4)
# points(andy.bai$avg.dbh, exp(b0.bai+b2.bai)*exp(b1.bai*log(andy.bai$avg.dbh)), 
#        col=col.edge[2], pch=13, cex=0.4)
# legend(x=50, y=0.3, legend=c("Edge (<10m)", "Interior"), fill=col.edge, bty="n")
# 
# ## non-transformed space ### pseudo-replicates
# mod.int.edge.ps.nls <- nls(biom.rel.ann ~ exp(a + b * log(dbh.start)), data=ps.contain[is.finite(biom.rel.ann & dbh.start>5),], start=list(a=0, b=0)) 
# summary(mod.int.edge.ps.nls) ## pretty similar to the original average-dbh assessment
# ee <- summary(mod.int.edge.ps)

### scatter avg.growth~avg.dbh
col.edge <- c("royalblue", "purple")
plot(andy.bai$avg.dbh, andy.bai$growth.mean, 
     col=col.edge[andy.bai$seg.Edge],
     pch=15, cex=0.4, ylab="Growth Rate (kg/kg)", xlab="Stem DBH (cm)", main="Urban Forest stems (avg dbh)")
points(andy.bai$avg.dbh, 
       exp(b0.bai+(b1.bai*log(andy.bai$avg.dbh))), 
       col=col.edge[1], pch=13, cex=0.6)
points(andy.bai$avg.dbh, 
       exp(b2.bai+b0.bai+(b1.bai*log(andy.bai$avg.dbh))), 
       col=col.edge[2], pch=13, cex=0.6)
legend(x=50, y=0.3, legend=c("Edge (<10m)", "Interior"), fill=col.edge, bty="n")

### with pseudoreps
col.edge <- c("royalblue", "purple")
plot(ps.contain[dbh.start>5, dbh.start], ps.contain[is.finite(biom.rel.ann) & dbh.start>5, biom.rel.ann], 
     col=col.edge[ps.contain[dbh.start>5,seg.Edge]],
     pch=15, cex=0.3, ylab="Growth Rate (kg/kg)", xlab="Stem DBH (cm)", main="Urban Forest stems (pseudoreps)")
points(ps.contain[dbh.start>5, dbh.start], 
       exp(b0.bai.ps+(b1.bai.ps*ps.contain[dbh.start>5, log(dbh.start)])), 
       col=col.edge[1], pch=13, cex=0.4)
points(ps.contain[dbh.start>5, dbh.start], 
       exp(b2.bai.ps+b0.bai.ps+(b1.bai.ps*ps.contain[dbh.start>5, log(dbh.start)])), 
       col=col.edge[2], pch=13, cex=0.4)
legend(x=50, y=0.3, legend=c("Edge (<10m)", "Interior"), fill=col.edge, bty="n")
## ok, minimal disturbance to the overall impression of growth rate by stem, should have minor impact on areal calcs

####
### bring in full andy.dbh data set, groom to determine growth rates on area/edge basis
ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv")
andy.dbh <- as.data.table(andy.dbh)
names(andy.dbh) <- c("Tree.ID", "Plot.ID", "Spp", "X", "Y", "dbh", "Can.class")
andy.dbh[Y<10, seg:=10]
andy.dbh[Y>=10 & Y<20, seg:=20]
andy.dbh[Y>=20, seg:=30]
andy.dbh[,ba:=ba.pred(dbh)]
setkey(andy.dbh, Tree.ID)
andy.dbh[Spp=="FRAL", Spp:="FRAM"] # Fraxinus americana

### get biomass (kg) as function of dbh (cm) for all the trees in andy.dbh
### the following coefficients are from Chojnacky et al. 2014
biom.pred.key <- data.frame(unique(andy.dbh$Spp))
# For HAVI (Witch hazel), CRSp ?? and TASp ?? we will use Jenkins generalized hardwood coefficients
b0 <- c(-2.0705, -3.0506, -2.6177, -2.0705, -2.0705, -2.2118, -2.0470, -2.0705, -2.6327, -2.3480, -2.48, -2.2271, -1.8384, -2.48, -2.48, -1.8011, -2.2271)
b1 <- c(2.4410, 2.6465, 2.4638, 2.4410, 2.4410, 2.4133, 2.3852, 2.4410, 2.4757, 2.3876, 2.4835, 2.4513, 2.3524, 2.4835, 2.4835, 2.3852, 2.4513)
biom.pred.key <- cbind(biom.pred.key, b0, b1)
biom.pred.key <- as.data.table(biom.pred.key)
names(biom.pred.key) <- c("Spp", "b0", "b1")
for(b in 1:dim(biom.pred.key)[1]){
  andy.dbh[Spp==biom.pred.key[b, Spp], biom:=exp(biom.pred.key[b, b0]+(biom.pred.key[b,b1]*log(dbh)))]
}
# hist(andy.dbh[, biom]) ## a couple of monsters

### where we get the generalized biomass gain per kg in edge vs. interior forest
andy.dbh[seg==10, growth.kg:=biom*exp(b0.bai+(b1.bai*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth.kg:=biom*exp(b0.bai+(b1.bai*log(dbh))+b2.bai)]
### growth as calculated from the dbh~growth including pseudoreplicates
andy.dbh[seg==10, growth.kg.ps:=biom*exp(b0.bai.ps+(b1.bai.ps*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth.kg.ps:=biom*exp(b0.bai.ps+(b1.bai.ps*log(dbh))+b2.bai.ps)]
## export the processed andy.dbh figures to make life easier elsewhere
write.csv(andy.dbh, "processed/andy.dbh.proc.results.csv")
andy.dbh <- read.csv("processed/andy.dbh.proc.results.csv")
andy.dbh <- as.data.table(andy.dbh)

### what is areal-basis growth, growth~biomass(kg/ha)
par(mar=c(4,4,1,1))
## is there a usable edge/int model for growth on an areal basis?
g <- andy.dbh[, .(sum(growth.kg), sum(growth.kg.ps), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
names(g)[3:5] <- c("growth.kg", "growth.kg.ps", "biom")
g[,rel.gain:=growth.kg/biom] ## this deals in forest growth as a function of biomass, not of forest area
g[,rel.gain.ps:=growth.kg.ps/biom]
plot(g$biom, g$rel.gain, col=as.numeric(g$seg), xlab="total plot biomass", main="avg.dbh")
plot(g$biom, g$rel.gain.ps, col=as.numeric(g$seg), xlab="total plot biomass", main="pseudoreps")
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
andy.growth.log <- lm(log(g$rel.gain.ps)~log(biom), data=g) ## R2 0.16, barely significant
andy.growth.log.edge <- lm(log(g$rel.gain.ps)~log(biom)*seg.F, data=g) ## R2 0.75, all factors signifiant
andy.plot.mod <- summary(andy.growth.log.edge)
test <- seq(2000, 12000, by=100)
points(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))),
     col="red", pch=12, cex=0.3)
points(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))+andy.plot.mod$coefficients[3]+(andy.plot.mod$coefficients[4]*log(test))),
       col="black", pch=12, cex=0.3)

g[,biom.kg.ha:=(biom/(10*20))*1E4] ## get things in a standard biomass density
g$biom.kg.ha/(1000*2) ### OK this is about the range of our forested pixels in terms of MgC/ha
plot(g$biom.kg.ha, g$rel.gain.ps, col=as.numeric(g$seg), xlab="total plot biomass", main="pseudoreps")
andy.growth.log.edge <- lm(log(g$rel.gain.ps)~log(biom.kg.ha)*seg.F, data=g) ## R2 0.75, sig
andy.plot.mod <- summary(andy.growth.log.edge)
test <- seq(1E5, 6E5, by=1E4)
points(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))),
       col="red", pch=12, cex=0.3)
points(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))+andy.plot.mod$coefficients[3]+(andy.plot.mod$coefficients[4]*log(test))),
       col="black", pch=12, cex=0.3)

## well that seems to fit just dandy. FUck. So slope and intercepts vary with edge/interior and biomass density


### ANDY NPP CALC 0: USE MEAN GROWTH RATE FOR EDGE/INTERIOR
edint.fact <- g[,mean(rel.gain), by=seg.F] ## 4.34% edge vs. 2.85% interior, annual gain on biomass per edge category
edint.fact.ps <- g[,mean(rel.gain.ps), by=seg.F] ## 3.91% edge vs. 2.74% interior, annual gain on biomass per edge category
### apply the growth factors to the edge/interior biomass fractions
biom.dat[,npp.edge:=ed.biom*edint.fact[seg.F=="E", V1]] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int:=int.biom*edint.fact[seg.F=="I", V1]] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot:=npp.edge+npp.int]
biom.dat[,range(npp.tot, na.rm=T)] ## 0 1443 kg biomass/cell/yr
biom.dat[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]

### same for pseudorep coefficients
biom.dat[,npp.edge.ps:=ed.biom*edint.fact.ps[seg.F=="E", V1]] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int.ps:=int.biom*edint.fact.ps[seg.F=="I", V1]] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot.ps:=npp.edge.ps+npp.int.ps]
biom.dat[,range(npp.tot.ps, na.rm=T)] ## 0 1404 kg biomass/cell/yr
biom.dat[,npp.tot.ps.MgC.ha:=((npp.tot.ps*1E-3*(1/2))/aoi)*1E4]

#### ANDY NPP CALC 1: APPLY PLOT-LEVEL MODEL OF growth~biomass*Edge
### range of biomass density in Boston data is up to 51k kg/cell, equiv to 567k kg-biomass/ha, 283 MgC/ha
g[,biom.kg.ha:=(biom/(10*20))*1E4] ## get things in a standard biomass density
g$biom.kg.ha/(1000*2) ### OK this is about the range of our forested pixels in terms of MgC/ha
andy.growth.log.edge <- lm(log(g$rel.gain.ps)~log(biom.kg.ha)*seg.F, data=g) ## R2 0.75, sig
andy.plot.mod <- summary(andy.growth.log.edge)

### now apply mdoel to the biomass data
# biom.dat[,range(ed.biom, na.rm=T)]
# hist(biom.dat$ed.biom)
biom.dat[,biom.kg.ha:=(biom/aoi)*1E4]
# hist(biom.dat[aoi>800, biom.kg.ha]) ## OK
biom.dat[,ed.biom.kg.ha:=(ed.biom/(aoi*ed.can))*1E4]
hist(biom.dat[aoi>800, ed.biom.kg.ha])
biom.dat[,int.biom.kg.ha:=(int.biom/(aoi*int.can))*1E4]
biom.dat[,edge.mod.factor:=exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(ed.biom.kg.ha)))] ## apply growth model for edge to edge biomass density
summary(biom.dat$edge.mod.factor) ## 5 to 9% in middle quartiles
hist(biom.dat[aoi>800 & edge.mod.factor<0.4, edge.mod.factor])
biom.dat[,int.mod.factor:=exp(andy.plot.mod$coefficients[1]+andy.plot.mod$coefficients[3]+((andy.plot.mod$coefficients[2]+andy.plot.mod$coefficients[4])*log(int.biom.kg.ha)))] ## apply growth model for interior growth to interior biomass density
summary(biom.dat$int.mod.factor) ## in the range of 2%
hist(biom.dat[aoi>800, int.mod.factor])

biom.dat[,npp.edge.mod:=ed.biom*edge.mod.factor] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int.mod:=int.biom*int.mod.factor] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot.mod:=npp.edge.mod+npp.int.mod]
biom.dat[,range(npp.tot.mod, na.rm=T)] ## 0 1031 kg biomass/cell/yr
biom.dat[,npp.tot.mod.MgC.ha:=((npp.tot.mod*1E-3*(1/2))/aoi)*1E4]


## how do the static (mean productivity edge/interior) vs. the modeled productivity (growth~biomass*edge) compare in aggregate?
biom.dat[,sum(npp.tot.ps, na.rm=T)]*1E-3*(1/2) #12.7k tC, a full 1k tC less than using the average dbh coefficients
biom.dat[,sum(npp.tot.mod, na.rm=T)]*1E-3*(1/2) #8.9k tC, about 5k tC less than using static values (basically we applied too-high growth factor to the densest edge forest)

(biom.dat[,sum(npp.tot.ps, na.rm=T)]*1E-3*(1/2))/(biom.dat[,sum(aoi, na.rm=T)]*1E-4) ## ie 1.02 tC/ha/yr using static
(biom.dat[,sum(npp.tot.mod, na.rm=T)]*1E-3*(1/2))/(biom.dat[,sum(aoi, na.rm=T)]*1E-4) ## ie 0.71 tC/ha/yr using modeled

hist(biom.dat[aoi>800 & biom>10, npp.tot.ps]) ## up to 1500 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.ps.MgC.ha]) ## up to 8 MgC/ha/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod]) ## up to 1000 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.MgC.ha]) ## up to 6 MgC/ha/yr
## so that's interesting. A lot produce not much, but there's a longer tail of high productivity that starts to look more like real forest -- up to ~8 MgC/ha/yr
## an ancillary question: why does FIA have such a pessimistic idea of productivity, maxes out at 2.5 MgC/ha/yr
write.csv(biom.dat, "processed/andy.forest.results.csv")
"processed/andy.dbh.proc.results.csv"
## noise testing to look at spread around the model coefficients here
### select randomly the coefficients to within +/- 1 sterr of the model estimates
mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]

b0.bai.range <- seq(mod.int.edge.fin$coefficients[1,1]-mod.int.edge.fin$coefficients[1,2], 
                    mod.int.edge.fin$coefficients[1,1]+mod.int.edge.fin$coefficients[1,2],
                    length.out=50)
b1.bai.range <- seq(mod.int.edge.fin$coefficients[2,1]-mod.int.edge.fin$coefficients[2,2], 
                    mod.int.edge.fin$coefficients[2,1]+mod.int.edge.fin$coefficients[2,2],
                    length.out=50)
b2.bai.range <- seq(mod.int.edge.fin$coefficients[3,1]-mod.int.edge.fin$coefficients[3,2], 
                    mod.int.edge.fin$coefficients[3,1]+mod.int.edge.fin$coefficients[3,2],
                    length.out=50)
dbh.dump <- andy.dbh
biom.dump <- biom.dat
npp.track <- numeric()
beta.track <- list()
edge.track <- numeric()
int.track <- numeric()
for(i in 1:2000){
  b0.rand <- sample(b0.bai.range, 1)
  b1.rand <- sample(b1.bai.range, 1)
  b2.rand <- sample(b2.bai.range, 1)
  dbh.dump[seg==10, growth:=biom*exp(b0.rand+(b1.rand*log(dbh)))]
  dbh.dump[seg%in%c(20,30), growth:=biom*exp(b0.rand+(b1.rand*log(dbh))+b2.rand)]
  g <- dbh.dump[, .(sum(growth), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
  g[,rel.gain:=V1/V2] ## this deals in forest growth as a function of biomass, not of forest area
  edge.factor <- g[seg==10, mean(rel.gain)]
  int.factor <- g[seg!=10, mean(rel.gain)]
  biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
  biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
  biom.dump[,npp.tot:=npp.edge+npp.int]
  biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
  npp.track <- c(npp.track, biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)])
  beta.track[[i]] <- c(b0.rand, b1.rand, b2.rand)
  edge.track <- c(edge.track, edge.factor)
  int.track <- c(int.track, int.factor)
  print(i)
}

npp.track[which(npp.track==min(npp.track))]/(1E3*2) ### minimum of 8.3k tC
npp.track[which(npp.track==max(npp.track))]/(1E3*2) ## max of 22.7k tC
b0.bai #-0.037
b1.bai #-0.916
b2.bai #-0.490
beta.track[[which(npp.track==min(npp.track))]] # -0.28, -0.99 -0.51
beta.track[[which(npp.track==max(npp.track))]] # 0.21, -0.84, -0.45
a <- edge.track-int.track
max(a) ## can be up to 2.84% apart
min(a) # can be as little as 0.7% apart

## at maximum plausible edge vs. interior productivity
edge.factor <- max(edge.track)
int.factor <- min(int.track)
biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
biom.dump[,npp.tot:=npp.edge+npp.int]
biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)]*1E-3*(1/2) ## 19.8k tC at maximum separation

## at minimum plausible edge vs. interior productivity
edge.factor <- min(edge.track)
int.factor <- max(int.track)
biom.dump[,npp.edge:=ed.biom*edge.factor] ## apply growth of edge trees to edge biomass
biom.dump[,npp.int:=int.biom*int.factor] ## apply growth of interior trees to interior biomass
biom.dump[,npp.tot:=npp.edge+npp.int]
biom.dump[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
biom.dump[aoi>800 & biom>10, sum(npp.tot, na.rm=T)]*1E-3*(1/2) ## 11.7k tC at minimum separation

save(beta.track, file="processed/andy.forest.beta.samples")
write.csv(cbind(npp.track, edge.track, int.track), "processed/andy.forest.npp.edge.int.samples.csv")

