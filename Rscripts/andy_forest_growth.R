library(raster)
library(data.table)


##### ANALYSIS OF TREE CORE RECORDS
## 
andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv") 
andy.bai <- as.data.table(andy.bai)
### the basic allometrics to get biomass
## Jenkins, C.J., D.C. Chojnacky, L.S. Heath and R.A. Birdsey. 2003. Forest Sci 49(1):12-35.
## Chojnacky, D.C., L.S. Heath and J.C. Jenkins. 2014. Forestry 87: 129-151.

#### limited species represented in the tree core samples
# Pinus rigida, Pinus strobus, both spg>0.45, == PIRI, PIST
# Acer rubrum, Aceraceae <0.50 spg, == ACRU
# Quercus alba, Quercus coccinea, Quercus rubra, Quercus velutina --> deciduous Fagaceae == QUAL, QUCO, QURU, QUVE
b0.l <- c(-3.0506, -2.0470, -2.0705) ## allometric coefficients (Pinus, Acer, Quercus)
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
biom.hist <- data.table() ## the running biomass history of each tree, according to the cores
## figure out biomass growth history applying species-specific allometrics
for(sp in 1:3){
  tmp <- andy.bai[Spp %in% taxa[[sp]],]
  a <- tmp[, lapply(tmp[,7:33], function(x){biom.pred(x, b0.l[sp], b1.l[sp])})]
  a <- cbind(tmp[,.(Tree.ID, incr.ID, Spp, X, Y)], a)
  biom.hist <- rbind(biom.hist, a)
}
biom.hist$Spp <- as.character(biom.hist$Spp) 
names(biom.hist)[6:32] <- paste0("biom", 2016:1990)

### lagged biomass gains, relative to previous year biomass
biom.rel <- data.frame(c(2015:1990))
for(i in 1:dim(biom.hist)[1]){
  te <- unlist(biom.hist[i,6:32])
  te.di <- (-1*diff(te, lag=1))
  biom.rel[,i+1] <- te.di/te[-1]
}
names(biom.rel)[-1] <- paste0("Tree", biom.hist[,incr.ID]) 
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

#######
### initial approach to getting a consistent growth~dbh relationship: 
### take mean of rel. biomass gain and mean of dbh across all the years present in each core sample
### figure the mean relative biomass gain per year across all years in each tree
mean.na <- function(x){mean(x, na.rm=T)}
m <- apply(biom.rel[,-1], FUN=mean.na, 2)
growth.mean <- data.frame(cbind(c(sub("Tree", '', colnames(biom.rel[,-1]))), m))
names(growth.mean) <- c("incr.ID", "growth.mean")
growth.mean$incr.ID <- as.integer(as.character(growth.mean$incr.ID))
andy.bai <- merge(andy.bai, growth.mean, by="incr.ID")
andy.bai$growth.mean <- as.numeric(as.character(andy.bai$growth.mean))
andy.bai[,avg.dbh:=apply(as.matrix(andy.bai[,8:34]), FUN=mean.na, 1)]

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
write.csv(andy.bai, "processed/andy.bai.dbh.avg.csv")

### upgrade Jul 23: treat five-year time slices as pseudo-replicates
pseudo.A <- matrix(ncol=4, rep(999,4))
pseudo.B <- matrix(ncol=4, rep(999,4))
pseudo.C <- matrix(ncol=4, rep(999,4))
pseudo.D <- matrix(ncol=4, rep(999,4))
pseudo.E <- matrix(ncol=4, rep(999,4)) ## so we have growth histories going back up to 25 years

### for each tree core, calculate average annualized (forward) relative growth in successive 5-year chunks, moving chunks backward until you run out of rings
for(d in unique(andy.bai$incr.ID)){  ### the intervals below are non-overlapping: each 5 year chunk is a distinct biomass time series moving backward
  pseudo.A <- rbind(pseudo.A, c(((biom.hist[incr.ID==d, biom2016]-biom.hist[incr.ID==d, biom2012])/biom.hist[incr.ID==d, biom2012])/5,
                                andy.bai[incr.ID==d, dbh2012],
                                andy.bai[incr.ID==d, incr.ID],
                                "A"))
  pseudo.B <- rbind(pseudo.B, c(((biom.hist[incr.ID==d, biom2011]-biom.hist[incr.ID==d, biom2007])/biom.hist[incr.ID==d, biom2007])/5,
                                andy.bai[incr.ID==d, dbh2007],
                                andy.bai[incr.ID==d, incr.ID],
                                "B"))
  pseudo.C <- rbind(pseudo.C, c(((biom.hist[incr.ID==d, biom2006]-biom.hist[incr.ID==d, biom2002])/biom.hist[incr.ID==d, biom2002])/5,
                                andy.bai[incr.ID==d, dbh2002],
                                andy.bai[incr.ID==d, incr.ID],
                                "C"))
  pseudo.D <- rbind(pseudo.D, c(((biom.hist[incr.ID==d, biom2001]-biom.hist[incr.ID==d, biom1997])/biom.hist[incr.ID==d, biom1997])/5,
                                andy.bai[incr.ID==d, dbh1997],
                                andy.bai[incr.ID==d, incr.ID],
                                "D"))

  ### pseudo.E can have from 5-7 years in a successsful sample
  tmp <- unlist(biom.hist[incr.ID==d, 26:32, with=F])
  t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
  if(t<5 | !is.finite(t)){ ## if record ends 1992-1996 (i.e. does not give at least 5 year run)
    pseudo.E <- rbind(pseudo.E, c(999,999,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    print("not enough data for 96-90 chunk") ## this record not long enough
  }else{
    g <- ((biom.hist[incr.ID==d, biom1996]-biom.hist[incr.ID==d, (26+(t-1)), with=F])/biom.hist[incr.ID==d, (26+(t-1)), with=F])/sum(is.finite(tmp))
    deeb <- min(unlist(andy.bai[incr.ID==d, 28:34, with=F]), na.rm=T)
    pseudo.E <- rbind(pseudo.E, c(g, deeb,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    print("retreived 96-90 chunk")
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
ps.contain <- ps.contain[order(incr.ID),]
ps.contain <- ps.contain[dbh.start!=999,]
ps.contain <- merge(x=ps.contain, y=andy.bai[,.(incr.ID, seg, seg.Edge)], by="incr.ID", all.x=T, all.y=F)

## how's it look?
par(mfrow=c(1,1))
ps.contain[is.finite(biom.rel.ann) & dbh.start>5,] ## 903 reasonably big trees 
plot(ps.contain$dbh.start, ps.contain$biom.rel.ann, col=ps.contain$interval, ylim=c(0,1))
plot(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, dbh.start], 
     ps.contain[is.finite(biom.rel.ann) & dbh.start>5, biom.rel.ann],
     col=ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge])
boxplot(biom.rel.ann~interval, data=ps.contain, ylim=c(0, 0.5)) ## forest seems to grow faster in its younger years back in the heady early 90's

### dump the pseudo-replicated file to disk
write.csv(ps.contain, "processed/andy.bai.dbh.pseudo.csv")



#########
### STEM GROWTH~DBH ANALYSIS
andy.bai <- read.csv("processed/andy.bai.dbh.avg.csv")
andy.bai <- as.data.table(andy.bai)
ps.contain <- read.csv("processed/andy.bai.dbh.pseudo.csv")
ps.contain <- as.data.table(ps.contain)
ps.contain <- merge(x=ps.contain, y=andy.bai[,.(incr.ID, Plot.ID)], by="incr.ID", all.x=T, all.y=F) ## put the plot IDs in
write.csv(ps.contain, "processed/andy.bai.dbh.pseudo.csv")

## avg.dbh~avg.growth
par(mfrow=c(1,2))
col.edge <- c("black", "red")
plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean),
     col=as.numeric(andy.bai$seg.Edge), main="Andy trees, avg. dbh")
summary(andy.bai$growth.mean) ### around 4%
summary(andy.bai$avg.dbh) ## about 19cm
table(andy.bai$seg.Edge) ## 64 edge: 131 interior

## vs. pseudoreps
plot(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, log(dbh.start)],
     ps.contain[is.finite(biom.rel.ann) & dbh.start>5, log(biom.rel.ann)],
     col=ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge],
     main="Andy trees, dbh pseudoreps")
table(ps.contain[is.finite(biom.rel.ann) & dbh.start>5, seg.Edge]) ## 270 edge:633 interior

summary(ps.contain$biom.rel.ann) ## comparable in the middle, pseudos a bit lower and have greater range
summary(ps.contain$dbh.start) ## about 19cm
plot(ps.contain[incr.ID==1, dbh.start], ps.contain[incr.ID==1, biom.rel.ann])
plot(ps.contain[incr.ID==4, dbh.start], ps.contain[incr.ID==4, biom.rel.ann])
plot(ps.contain[incr.ID==11, dbh.start], ps.contain[incr.ID==11, biom.rel.ann])
plot(ps.contain[incr.ID==100, dbh.start], ps.contain[incr.ID==100, biom.rel.ann]) ## Ok, comforting. Big trees grow slowly, and within tree growth declines with time/dbh

### modeling growth~dbh*edge
### full interactive dbh*edge, avg.dbh and pseudo
mod.int.edge <- summary(lm(log(growth.mean)~log(avg.dbh)*seg.Edge, data=andy.bai)) #r2=0.47, no sig. dbh*edge interaction
mod.int.edge.ps <- summary(lm(log(biom.rel.ann)~log(dbh.start)*seg.Edge, data=ps.contain[dbh.start>=5,])) #r2=0.31, no sig. dbh*edge interaction

## removing the slope modifier (i.e. just a different intercept for edge vs. interior)
mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai)) #r2=0.46, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]

## model using pseudoreplicates
mod.int.edge.fin.ps <- summary(lm(log(biom.rel.ann)~log(dbh.start)+seg.Edge, data=ps.contain[dbh.start>=5,])) #r2=0.30, just edge vs. interior
b0.bai.ps <- mod.int.edge.fin.ps$coefficients[1]
b1.bai.ps <- mod.int.edge.fin.ps$coefficients[2]
b2.bai.ps <- mod.int.edge.fin.ps$coefficients[3]
mod.int.edge.fin.ps$coefficients

#### Try to do some mixed effects modeling to handle the pseudoreplication issue
library(lme4)
par(mfrow=c(1,1))
ps.contain[, incr.ID:=as.factor(incr.ID)]
boxplot(biom.rel.ann~interval, ps.contain, ylim=c(0, 0.4))
boxplot(biom.rel.ann~incr.ID, ps.contain, ylim=c(0,0.4))
boxplot(biom.rel.ann~Plot.ID, ps.contain, ylim=c(0,0.4))


## edge slope, random intercepts on plot increment interval
mmod2 <- lmer(log(biom.rel.ann)~log(dbh.start)*seg.Edge + (1|Plot.ID) + (1|incr.ID) + (1|interval), 
              data=ps.contain[dbh.start>=5,],
              REML=FALSE)
summary(mmod2)
## edge slope, random slopes on plot increment interval
mmod.Rslopes.Finteract <- lmer(log(biom.rel.ann)~log(dbh.start)*seg.Edge + (1+log(dbh.start)|Plot.ID) + (1+log(dbh.start)|incr.ID) + (1+log(dbh.start)|interval), 
              data=ps.contain[dbh.start>=5,],
              REML=FALSE)
summary(mmod.Rslopes.Finteract) 

## edge intercept, random slopes on plot increment interval
mmod.Rslopes.Fadd <- lmer(log(biom.rel.ann)~log(dbh.start)+seg.Edge + (1+log(dbh.start)|Plot.ID) + (1+log(dbh.start)|incr.ID) + (1+log(dbh.start)|interval), 
                               data=ps.contain[dbh.start>=5,],
                               REML=FALSE)
summary(mmod.Rslopes.Fadd)
coef(mod.int.edge.fin.ps) ## not particularly different from the fixed effects model

### significance of dbh*Edge
anova(mmod.Rslopes.Fadd, mmod.Rslopes.Finteract) ## p>0.33, dbh:Edge not significant

mmod.Rslopes.Fnull <- lmer(log(biom.rel.ann)~log(dbh.start) + (1+log(dbh.start)|Plot.ID) + (1+log(dbh.start)|incr.ID) + (1+log(dbh.start)|interval), 
                           data=ps.contain[dbh.start>=5,],
                           REML=FALSE)
summary(mmod.Rslopes.Fnull)
### significance of Edge
anova(mmod.Rslopes.Fnull, mmod.Rslopes.Fadd) ## p < 0.001, edge factor is significant

mmod.Rslopes.null <- lmer(log(biom.rel.ann)~1 + (1+log(dbh.start)|Plot.ID) + (1+log(dbh.start)|incr.ID) + (1+log(dbh.start)|interval), 
                          data=ps.contain[dbh.start>=5,],
                          REML=FALSE)
summary(mmod.Rslopes.null)
### significance of dbh
anova(mmod.Rslopes.null, mmod.Rslopes.Fnull) ## p < 0.01, dbh is significant

boxplot(biom.rel.ann~seg.Edge*Plot.ID, data=ps.contain[dbh.start>=5,])
### Mungo model, random effects on both edge and dbh, full dbh*Edge
mmod.Mungo.Finteract <- lmer(log(biom.rel.ann)~log(dbh.start)*seg.Edge + 
                               (1+log(dbh.start)|Plot.ID) + 
                               (1+log(dbh.start)|incr.ID) + 
                               (1+log(dbh.start)|interval) +
                               (1+seg.Edge|Plot.ID) +
                               # (1+seg.Edge|incr.ID) + ## there are no individual stems that are in both an edge and interior
                               (1+seg.Edge|interval), 
                               data=ps.contain[dbh.start>=5,],
                               REML=FALSE)
summary(mmod.Mungo.Finteract)
summary(mmod.Rslopes.Finteract) ## pretty similar
anova(mmod.Rslopes.Finteract, mmod.Mungo.Finteract) ### Mungo seems a bit better, p<0.02

mmod.Mungo.Fadd <- lmer(log(biom.rel.ann)~log(dbh.start)+seg.Edge + 
                          (1+log(dbh.start)|Plot.ID) + 
                          (1+log(dbh.start)|incr.ID) + 
                          (1+log(dbh.start)|interval) +
                          (1+seg.Edge|Plot.ID) +
                          # (1+seg.Edge|incr.ID) + ## there are no individual stems that are in both an edge and interior
                          (1+seg.Edge|interval), 
                        data=ps.contain[dbh.start>=5,],
                        REML=FALSE)


summary(mmod.Mungo.Fadd)
summary(mmod.Rslopes.Fadd) ### similar coefficients
coef(mod.int.edge.fin.ps)
### does Mungo do better than random effects on dbh?
anova(mmod.Rslopes.Fadd, mmod.Mungo.Fadd) ## bit better, p<0.03
### is dbh*edge significant in Mungo?
anova(mmod.Mungo.Fadd, mmod.Mungo.Finteract) ### not sig, p>0.22
## model diagnostics
plot(ps.contain[dbh.start>=5, log(dbh.start)], residuals(mmod.Mungo.Fadd)) ## homoskedastic
hist(residuals(mmod.Mungo.Fadd)) ### looks normal
qqnorm(residuals(mmod.Mungo.Fadd)); qqline (residuals(mmod.Mungo.Fadd), col=2) ## light tailed



### let's try the non-linear model version of this
## this is fucking harder
library(nlme)
pray <- nlme(biom.rel.ann~exp(a+b*log(dbh.start)), data=ps.contain[dbh.start>=5,],
         random = a ~ 1|incr.ID, fixed=list(a~1, b~1), start=c(0,0))
summary(pray)
t=1:100
plot(t, exp(1.24+(-1.54*log(t))))

plot(log(andy.bai$avg.dbh), log(andy.bai$growth.mean), pch=13)
plot(log(ps.contain$dbh.start), log(ps.contain$biom.rel.ann), col="blue", pch=15, cex=0.6)

library(lme4)
pray <- nlmer(biom.rel.ann ~ exp(a+b*log(dbh.start)) ~ a|incr.ID,
              data=ps.contain[dbh.start>=5,], start=c(0,0))

### generalized growth coefficients for edge/interior, based on andy.bai
par(mfrow=c(1,1), mar=c(4,4,2,1))
col.edge <- c("steelblue3", "orchid")
andy.xlim <- c(4, 100)
andy.ylim <- c(-0.15, 1)

### andy stem growth, avg.dbh analysis
# plot((andy.bai$avg.dbh), (andy.bai$growth.mean), 
#      col=col.edge[as.numeric(as.factor(andy.bai$seg.Edge))],
#      pch=15, cex=0.6, main="Urban Forest, stem avg. dbh", 
#      xlim=andy.xlim, ylim=andy.ylim, xlab="Avg. DBH (cm)", ylab="Avg. rel. growth (kg/kg)")
# test <- seq(andy.xlim[1], andy.xlim[2], length.out=100)
# lines(test, exp(b0.bai+(b1.bai*log(test))), lty=2, col="blue", lwd=2)
# lines(test, exp(b0.bai+(b1.bai*log(test))+b2.bai), lty=2, col="purple", lwd=2)
# abline(v=andy.bai[seg.Edge=="E", median(avg.dbh)], col="steelblue3")
# abline(v=andy.bai[seg.Edge=="I", median(avg.dbh)], col="orchid")
# abline(h=andy.bai[seg.Edge=="E", median(growth.mean)], col="steelblue3")
# abline(h=andy.bai[seg.Edge=="I", median(growth.mean)], col="orchid")
# legend(x=60, y=0.5, legend=c("Edge", "Interior"), fill=c("steelblue3", "orchid"), bty="n")

### andy stem growth, pseudorep analysis
plot((ps.contain$dbh.start), (ps.contain$biom.rel.ann), 
     col=col.edge[as.numeric(as.factor(ps.contain$seg.Edge))],
     pch=15, cex=0.6, main="Urban Forest", 
     xlim=andy.xlim, ylim=andy.ylim, xlab="Stem DBH (cm)", ylab="Relative growth (kg/kg)")
test <- seq(andy.xlim[1], andy.xlim[2], length.out=100)
lines(test, exp(b0.bai.ps+(b1.bai.ps*log(test))), lty=2, col="blue", lwd=2)
lines(test, exp(b0.bai.ps+(b1.bai.ps*log(test))+b2.bai.ps), lty=2, col="purple", lwd=2)
abline(v=ps.contain[seg.Edge=="E", median(dbh.start)], col="steelblue3")
abline(v=ps.contain[seg.Edge=="I", median(dbh.start)], col="orchid")
abline(h=ps.contain[seg.Edge=="E", median(biom.rel.ann)], col="steelblue3")
abline(h=ps.contain[seg.Edge=="I", median(biom.rel.ann)], col="orchid")
legend(x=70, y=0.6, legend=c("Edge (<10 m)", "Interior"), fill=c("steelblue3", "orchid"), bty="n")

# plot(log(ps.contain$dbh.start), log(ps.contain$biom.rel.ann), col=col.edge[as.numeric(as.factor(ps.contain$seg.Edge))],
#      pch=15, cex=0.5, main="Andy Trees, Pseudoreps")
# abline(b0.bai.ps, b1.bai.ps, col="black")
# abline(b0.bai.ps+b2.bai.ps, b1.bai.ps, col="red")
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

#####
#### PLOT-LEVEL GROWTH~BIOMASS ASSESSMENT BASED ON STEM MEASUREMENTS
### bring in full andy.dbh data set and apply the growth factors based on cores to the whole (exhaustive) dbh record for each plot
andy.bai <- read.csv("processed/andy.bai.dbh.avg.csv")
andy.bai <- as.data.table(andy.bai)
ps.contain <- read.csv("processed/andy.bai.dbh.pseudo.csv")
ps.contain <- as.data.table(ps.contain)

mod.int.edge.fin <- summary(lm(log(growth.mean)~log(avg.dbh)+seg.Edge, data=andy.bai[avg.dbh>=5,])) #r2=0.46, just edge vs. interior
b0.bai <- mod.int.edge.fin$coefficients[1]
b1.bai <- mod.int.edge.fin$coefficients[2]
b2.bai <- mod.int.edge.fin$coefficients[3]

mod.int.edge.fin.ps <- summary(lm(log(biom.rel.ann)~log(dbh.start)+seg.Edge, data=ps.contain[dbh.start>=5,])) #r2=0.30, just edge vs. interior
b0.bai.ps <- mod.int.edge.fin.ps$coefficients[1]
b1.bai.ps <- mod.int.edge.fin.ps$coefficients[2]
b2.bai.ps <- mod.int.edge.fin.ps$coefficients[3]

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

### apply growth function to predicted per-stem biomass according to edge/interior position
### using avg.dbh coefficients
andy.dbh[seg==10, growth.kg:=biom*exp(b0.bai+(b1.bai*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth.kg:=biom*exp(b0.bai+(b1.bai*log(dbh))+b2.bai)]

### growth as calculated from the dbh~growth including pseudoreplicates
andy.dbh[seg==10, growth.kg.ps:=biom*exp(b0.bai.ps+(b1.bai.ps*log(dbh)))]
andy.dbh[seg%in%c(20,30), growth.kg.ps:=biom*exp(b0.bai.ps+(b1.bai.ps*log(dbh))+b2.bai.ps)]

## export the processed andy.dbh figures to make life easier elsewhere
write.csv(andy.dbh, "processed/andy.dbh.proc.results.csv")



### Determine areal-basis growth, growth~biomass(kg/ha) relationship
andy.dbh <- read.csv("processed/andy.dbh.proc.results.csv")
andy.dbh <- as.data.table(andy.dbh)
g <- andy.dbh[, .(sum(growth.kg), sum(growth.kg.ps), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
names(g)[3:5] <- c("growth.kg", "growth.kg.ps", "biom")
g[,rel.gain:=growth.kg/biom] ## this deals in forest growth as a function of biomass, not of forest area
g[,rel.gain.ps:=growth.kg.ps/biom]
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
g$seg.F <- as.factor(g$seg.F)
g[,biom.Mg.ha:=(biom/(10*20))/1000*1E4] ## get things in a standard biomass density Mg.ha based on plot size
g$biom.Mg.ha/(2) ### OK this is about the range of our forested pixels in terms of MgC/ha

## what does plot-level relationships look like?
par(mar=c(4,4,1,1), mfrow=c(1,2))
plot(g$biom, g$rel.gain, col=as.numeric(g$seg), xlab="total plot biomass", main="avg.dbh")
plot(g$biom, g$rel.gain.ps, col=as.numeric(g$seg), xlab="total plot biomass", main="pseudoreps")


## plot-level models growth~biomass.Mg.ha
andy.growth.log <- lm(log(g$rel.gain.ps)~log(biom.Mg.ha), data=g) ## R2 0.1, barely significant
andy.growth.log.edge <- lm(log(rel.gain.ps)~log(biom.Mg.ha)*seg.F, data=g) ## R2 0.75, all factors sig
andy.plot.mod <- summary(andy.growth.log.edge)
andy.growth.lin <- lm(rel.gain.ps~biom.Mg.ha, data=g)
summary(andy.growth.lin) # R2 0.13, p>0.08
andy.growth.lin.edge <- lm(rel.gain.ps~biom.Mg.ha+seg.F, data=g) #R2 0.67, only biom*edge not significant
andy.plot.mod.lin <- summary(andy.growth.lin.edge)

par(mar=c(4,4,1,1), mfrow=c(1,1))
# plot(g$biom, g$rel.gain, col=as.numeric(g$seg), xlab="total plot biomass", main="avg.dbh")
col.edge <- c("steelblue3", "orchid")
andy.ylim <- c(-0.015, 0.083) ## trying to match up roughly with the range of the FIA plots
andy.xlim <- c(0, 650)
plot(g$biom.Mg.ha, g$rel.gain.ps, main="Urban Forest, plot growth",
     col=col.edge[as.numeric(g$seg.F)], pch=15, cex=0.9, ylim=andy.ylim, xlim=andy.xlim,
     xlab="Biomass density (Mg/ha)", ylab="Relative growth (Mg/Mg/ha)")
test <- seq(from=g[,min(biom.Mg.ha)], to = g[,max(biom.Mg.ha)], length.out=100)
lines(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))),
       col="blue", lty=2, lwd=2)
lines(test, exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(test))+andy.plot.mod$coefficients[3]+(andy.plot.mod$coefficients[4]*log(test))),
       col="purple",  lty=2, lwd=2)
legend(fill=c("steelblue3", "orchid"), legend=c("Edge (<10m)", "Interior"), x=450, y=0.05, bty="n")
abline(v=250, col="steelblue3") ## 90th percentile upper
abline(v=329, col="purple")## 90th percentile upper
abline(v=524, col="steelblue3", lty=3)
abline(v=597, col="purple", lty=3)

### compare to linear model
lines(test, andy.plot.mod.lin$coefficients[1]+(andy.plot.mod.lin$coefficients[2]*test), 
      col="blue", lty=3, lwd=2)
lines(test, andy.plot.mod.lin$coefficients[1]+andy.plot.mod.lin$coefficients[3]+(andy.plot.mod.lin$coefficients[2]*test), 
      col="purple", lty=3, lwd=2)
### the exponential fits fine but it gets v. nonlinear right in the low range of biomass density in edge pixels
### the linear is more modest about edge growth at low density, but doesn't fit as well and extrapolates to unrealistically low in interior forest at high biomass
### propose a compromise: Use the exponential fit and apply a cap on edge plots below the lowest density recorded = predicted productivity at that level


###### 
##### NPP Calculation treating canopy as Andy-like forest

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

### range of edge and interior biomass, in Mg-biomass/ha
biom.dat[aoi>800, range(ed.biom, na.rm=T)] ## in kg/cell
biom.dat[aoi>800, ed.biom.Mgbiom:=(ed.biom/1000)]
biom.dat[aoi>800, ed.biom.Mgbiom.ha.can:=(ed.biom.Mgbiom/(aoi*ed.can))*1E4]
biom.dat[aoi>800 & ed.can<0.005, ed.biom.Mgbiom.ha.can:=0]
hist(biom.dat[aoi>800 & ed.biom.Mgbiom.ha.can>0, ed.biom.Mgbiom.ha.can])
biom.dat[aoi>800, quantile(ed.biom.Mgbiom.ha.can, probs=c(0.05, 0.95), na.rm=T)]
biom.dat[aoi>800, range(ed.biom.Mgbiom.ha.can, na.rm=T)]
## what fraction below 100 but still there?
biom.dat[aoi>800 & ed.biom.Mgbiom.ha.can<=100 & ed.biom.Mgbiom.ha.can>5 & !is.na(ed.biom.Mgbiom.ha.can), length(ed.biom.Mgbiom.ha.can)]/biom.dat[aoi>800 & !is.na(ed.biom.Mgbiom.ha.can), length(ed.biom.Mgbiom.ha.can)]
## 28%

biom.dat[aoi>800, range(int.biom, na.rm=T)] ## in kg/cell
biom.dat[aoi>800, int.biom.Mgbiom:=(int.biom/1000)]
biom.dat[aoi>800, int.biom.Mgbiom.ha.can:=(int.biom.Mgbiom/(aoi*int.can))*1E4]
biom.dat[aoi>800 & int.can<0.005, int.biom.Mgbiom.ha.can:=0]
hist(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0, int.biom.Mgbiom.ha.can])
biom.dat[aoi>800, quantile(int.biom.Mgbiom.ha.can, probs=c(0.05, 0.95), na.rm=T)]
biom.dat[aoi>800, range(int.biom.Mgbiom.ha.can, na.rm=T)]
biom.dat[aoi>800 & int.biom.Mgbiom.ha.can<=100 & int.biom.Mgbiom.ha.can>5 & !is.na(int.biom.Mgbiom.ha.can), length(int.biom.Mgbiom.ha.can)]/biom.dat[aoi>800 & !is.na(int.biom.Mgbiom.ha.can), length(int.biom.Mgbiom.ha.can)]
## <1%

### ANDY NPP CALC 0: USE MEAN GROWTH RATE FOR EDGE/INTERIOR
andy.dbh <- read.csv("processed/andy.dbh.proc.results.csv")
andy.dbh <- as.data.table(andy.dbh)
g <- andy.dbh[, .(sum(growth.kg), sum(growth.kg.ps), sum(biom)), by=.(seg, Plot.ID)] ## total biomass gain and total biomass for each plot
names(g)[3:5] <- c("growth.kg", "growth.kg.ps", "biom")
g[,rel.gain:=growth.kg/biom] ## kg/kg growth rate judging from avg.dbh~growth stem data
g[,rel.gain.ps:=growth.kg.ps/biom] ## kg/kg growth rate as determined from pseudoreplicated stem data
g[seg==10, seg.F:="E"]
g[seg!=10, seg.F:="I"]
g$seg.F <- as.factor(g$seg.F)
g[,biom.Mg.ha:=(biom/(10*20))/1000*1E4] ## get things in a standard biomass density Mg.ha based on plot size (already at 100% canopy in these plots)

### apply the growth factors from avg.dbh analysis to the edge/interior biomass fractions
edint.fact <- g[,mean(rel.gain), by=seg.F] ## 4.29% edge vs. 2.85% interior, annual gain on biomass per edge category
g[,mean(rel.gain.ps), by=seg.F] ## contrast 3.1% edge vs. 2.2% interior based on pseudorep model for stem growth
biom.dat[,npp.edge:=ed.biom*edint.fact[seg.F=="E", V1]] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int:=int.biom*edint.fact[seg.F=="I", V1]] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot:=npp.edge+npp.int]
biom.dat[,range(npp.tot, na.rm=T)] ## 0 1465 kg biomass/cell/yr
biom.dat[,npp.tot.MgC.ha:=((npp.tot*1E-3*(1/2))/aoi)*1E4]
biom.dat[aoi>800, range(npp.tot.MgC.ha, na.rm=T)] ### up to ~8 MgC/ha
biom.dat[aoi>800, sum(npp.tot, na.rm=T)]/(2*1000) ##13.7k tC
biom.dat[aoi>800, (sum(npp.tot, na.rm=T)/(2*1000))/sum(aoi, na.rm=T)]*1E4 ## mean of 1.12 tC/ha

### same for growth factors derived from pseudorep analysis
edint.fact.ps <- g[,mean(rel.gain.ps), by=seg.F] ## 3.91% edge vs. 2.74% interior, annual gain on biomass per edge category
biom.dat[,npp.edge.ps:=ed.biom*edint.fact.ps[seg.F=="E", V1]] ## apply growth of edge trees to edge biomass
biom.dat[,npp.int.ps:=int.biom*edint.fact.ps[seg.F=="I", V1]] ## apply growth of interior trees to interior biomass
biom.dat[,npp.tot.ps:=npp.edge.ps+npp.int.ps]
biom.dat[,range(npp.tot.ps, na.rm=T)] ## 0 1122 kg biomass/cell/yr
biom.dat[,npp.tot.ps.MgC.ha:=((npp.tot.ps*1E-3*(1/2))/aoi)*1E4]
biom.dat[aoi>800, range(npp.tot.ps.MgC.ha, na.rm=T)] ## up to ~6 MgC/ha
biom.dat[aoi>800, sum(npp.tot.ps, na.rm=T)/(2*1000)] ## 10.1 tC
biom.dat[aoi>800, (sum(npp.tot.ps, na.rm=T)/(2*1000))/sum(aoi, na.rm=T)]*1E4 ### 0.82 tC/ha


#### ANDY NPP CALC 1: APPLY PLOT-LEVEL MODEL OF growth~biomass*Edge
### NOTE HERE: The Andy Forest calcs all make the correction for the canopy area of each cell
### i.e. biomass density is normalized based on the canopy coverage of each pixel (biomass kg/ha-canopy)
### because growth~biomass density was determined in closed-canopy systems, and biomass density as a predictor was determined based on the plot area being totally canopy covered
### To apply the model to partially-forested pixels, we area calculating pixel biomass density on basis of canopy coverage area per cell

### apply EXPONENTIAL model to the biomass data
andy.growth.log.edge <- lm(log(g$rel.gain.ps)~log(biom.Mg.ha)*seg.F, data=g) ## R2 0.75, sig
### not that this model is specified on the biomass density in g, i.e. in full canopy Andy forest plots (does not require canopy area correction)
andy.plot.mod <- summary(andy.growth.log.edge)
### use model to determine growth factor for edge and interior biomass fraction of each pixel
### determine growth factor for edge biomass
biom.dat[,edge.mod.factor:=exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(ed.biom.Mgbiom.ha.can)))] ## apply growth model for edge to edge biomass density CORRECTED FOR CANOPY AREA
biom.dat[ed.biom.Mgbiom.ha.can<0.1, edge.mod.factor:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
summary(biom.dat[aoi>800, edge.mod.factor]) ## 2.7 to 6% in middle quartiles, 0-1.37
### new hotness: cap the exponential prediction in low-density edge pixels
biom.dat[,edge.mod.factor.cap:=edge.mod.factor]
edge.cap <- exp(andy.plot.mod$coefficients[1]+(andy.plot.mod$coefficients[2]*log(g[seg.F=="E",min(biom.Mg.ha)])))
biom.dat[edge.mod.factor.cap>=edge.cap, edge.mod.factor.cap:=edge.cap]
summary(biom.dat[aoi>800, edge.mod.factor.cap]) ## reduces the median somewhat, cuts the top end
## compare
hist(biom.dat[aoi>800 & edge.mod.factor>0, edge.mod.factor])
hist(biom.dat[aoi>800 & edge.mod.factor.cap>0, edge.mod.factor.cap]) ## this will have the effect of treating most edge forest as a static factor
biom.dat[aoi>800 & edge.mod.factor<0.001, length(edge.mod.factor)] ## 33052 register as 0 growth
biom.dat[aoi>800 & edge.mod.factor.cap<0.001, length(edge.mod.factor.cap)] ## 33052 register as 0 growth
### determine growth factor in interior biomass
biom.dat[,int.mod.factor:=exp(andy.plot.mod$coefficients[1]+andy.plot.mod$coefficients[3]+((andy.plot.mod$coefficients[2]+andy.plot.mod$coefficients[4])*log(int.biom.Mgbiom.ha.can)))] ## apply growth model for interior growth to interior biomass density
biom.dat[int.biom.Mgbiom.ha.can<0.1, int.mod.factor:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
summary(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0, int.mod.factor]) ## in the range of 2%
hist(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0, int.mod.factor]) ## for that minority of cells with lots of interior biomass, below about 2.5
### kill a tiny number of dipshit high interior growth factors
int.cap <- exp(andy.plot.mod$coefficients[1]+andy.plot.mod$coefficients[3]+(andy.plot.mod$coefficients[2]*log(g[seg.F=="I",min(biom.Mg.ha)]))+(andy.plot.mod$coefficients[4]*log(g[seg.F=="I",min(biom.Mg.ha)])))
biom.dat[int.mod.factor>=int.cap, int.mod.factor:=int.cap]
summary(biom.dat[aoi>800 & int.mod.factor>0, int.mod.factor]) ## reduces the median somewhat, cuts the top end


### similarly, apply LINEAR model to biomass data
andy.growth.lin.edge <- lm(rel.gain.ps~biom.Mg.ha+seg.F, data=g) #R2 0.67, only biom*edge not significant
andy.plot.mod.lin <- summary(andy.growth.lin.edge)
### use model to determine growth factor for edge and interior biomass fraction of each pixel
### determine growth factor for edge biomass
biom.dat[,edge.mod.factor.lin:=andy.plot.mod.lin$coefficients[1]+(andy.plot.mod.lin$coefficients[2]*ed.biom.Mgbiom.ha.can)] ## apply growth model for edge to edge biomass density CORRECTED FOR CANOPY AREA
biom.dat[ed.biom.Mgbiom.ha.can<0.1, edge.mod.factor.lin:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
summary(biom.dat[aoi>800 & ed.biom.Mgbiom.ha.can, edge.mod.factor.lin]) ## 3.5 to 3.7 in middle quartiles, 2-4.3
### determine growth factor in interior biomass
biom.dat[,int.mod.factor.lin:=andy.plot.mod.lin$coefficients[1]+andy.plot.mod.lin$coefficients[3]+(andy.plot.mod.lin$coefficients[2]*int.biom.Mgbiom.ha.can)] ## apply growth model for interior growth to interior biomass density
biom.dat[int.biom.Mgbiom.ha.can<0.1, int.mod.factor.lin:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
summary(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0, int.mod.factor.lin]) ## 1.6 to 2.2, range 0.5 to 3.0
hist(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0, int.mod.factor]) ## for that minority of cells with lots of interior biomass, below about 2.5

### interrim diagnostics
biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0.1, length(int.biom.Mgbiom.ha.can)] ## only 24k cells have meaningful interior biomass
hist(biom.dat[aoi>800 & int.biom.Mgbiom.ha.can>0.1, int.biom.Mgbiom.ha.can]) ## that interior forest is pretty beefy, peaks 300 Mg/ha
hist(biom.dat[aoi>800 & ed.biom.Mgbiom.ha.can>0.1, ed.biom.Mgbiom.ha.can]) ## in contrast, the edge biomass tends to be less beefy, peaks about 100 Mg/ha

### static factor estimated NPP
biom.dat[,npp.edge.mod:=ed.biom*edge.mod.factor] ## apply growth of edge trees to edge biomass
biom.dat[,npp.edge.mod.cap:=ed.biom*edge.mod.factor.cap] ## apply alternative capped edge factor
biom.dat[,npp.edge.mod.lin:=ed.biom*edge.mod.factor.lin] ## apply linear model edge growth factor
biom.dat[,npp.int.mod:=int.biom*int.mod.factor] ## apply growth of interior trees to interior biomass
biom.dat[,npp.int.mod.lin:=int.biom*int.mod.factor.lin] ## apply linear model of interior growth factor
biom.dat[,npp.tot.mod:=npp.edge.mod+npp.int.mod]
biom.dat[,npp.tot.mod.cap:=npp.edge.mod.cap+npp.int.mod]
biom.dat[,npp.tot.mod.lin:=npp.edge.mod.lin+npp.int.mod.lin]

## range of estimated NPP, different approaches
biom.dat[,range(npp.tot, na.rm=T)] ## 0 1465 kg biomass/cell yr, static avg.dbh factors
biom.dat[,range(npp.tot.ps, na.rm=T)] ## 0 1122 kg biomass/cell/yr, static pseudorep factors
biom.dat[,range(npp.tot.mod, na.rm=T)] ## 0 848 kg biomass/cell/yr, exponential model
biom.dat[,range(npp.tot.mod.cap, na.rm=T)] ## 0 848 kg biomass/cell/yr, exponential model, capped
biom.dat[,range(npp.tot.mod.lin, na.rm=T)] ## 0 833 kg biomass/cell/yr, linear model

summary(biom.dat$npp.tot.mod)
summary(biom.dat$npp.tot.mod.cap) ## lowers the median a fair bit
summary(biom.dat$npp.tot.ps) ## but still a bit higher median than the static factors 
summary(biom.dat$npp.tot.mod.lin) ## really similar to the capped exponential model, but cleaner
biom.dat[,npp.tot.mod.MgC.ha:=((npp.tot.mod*1E-3*(1/2))/aoi)*1E4]
biom.dat[,npp.tot.mod.cap.MgC.ha:=((npp.tot.mod.cap*1E-3*(1/2))/aoi)*1E4]
biom.dat[,npp.tot.mod.lin.MgC.ha:=((npp.tot.mod.lin*1E-3*(1/2))/aoi)*1E4]


## uncapped model
biom.dat[aoi>800, range(npp.tot.mod.MgC.ha, na.rm=T)] ## up to ~4.7 MgC/ha
biom.dat[aoi>800, sum(npp.tot.mod, na.rm=T)/(2*1000)] ## 12.1 tC
biom.dat[aoi>800, (sum(npp.tot.mod, na.rm=T)/(2*1000))/sum(aoi, na.rm=T)]*1E4 ### 0.98 tC/ha

## capped model -- totals are a bit lower than the uncapped model (the extrapolations were having an impact)
biom.dat[aoi>800, range(npp.tot.mod.cap.MgC.ha, na.rm=T)] ## up to ~4.7 MgC/ha
biom.dat[aoi>800, sum(npp.tot.mod.cap, na.rm=T)/(2*1000)] ## 10.8 tC
biom.dat[aoi>800, (sum(npp.tot.mod.cap, na.rm=T)/(2*1000))/sum(aoi, na.rm=T)]*1E4 ### 0.88 tC/ha

## linear model -- damn near same as capped, slighly lower
biom.dat[aoi>800, range(npp.tot.mod.lin.MgC.ha, na.rm=T)] ## up to ~4.6 MgC/ha
biom.dat[aoi>800, sum(npp.tot.mod.lin, na.rm=T)/(2*1000)] ## 10.6 tC
biom.dat[aoi>800, (sum(npp.tot.mod.lin, na.rm=T)/(2*1000))/sum(aoi, na.rm=T)]*1E4 ### 0.86 tC/ha

## how do the static (mean productivity edge/interior) vs. the modeled productivity (growth~biomass*edge) compare in aggregate?
biom.dat[,sum(npp.tot, na.rm=T)]*1E-3*(1/2) ## 13.8 tC, original static avg.dbh growth factors
biom.dat[,sum(npp.tot.ps, na.rm=T)]*1E-3*(1/2) #10.1k tC, static pseudorep factors
biom.dat[,sum(npp.tot.mod, na.rm=T)]*1E-3*(1/2) #12.2k tC, pseudorep factors applied by model exponential growth~biomass*edge
biom.dat[,sum(npp.tot.mod.cap, na.rm=T)]*1E-3*(1/2) #10.8k tC, pseudorep factors applied by *capped* model growth~biomass*edge 
biom.dat[,sum(npp.tot.mod.lin, na.rm=T)]*1E-3*(1/2) #10.6k tC, pseudorep factors applied by *linear* model growth~biomass*edge 


(biom.dat[aoi>800,sum(npp.tot, na.rm=T)]*1E-3*(1/2))/(biom.dat[aoi>800,sum(aoi, na.rm=T)]*1E-4) ## ie 1.12 tC/ha/yr using static/avg.dbh
(biom.dat[aoi>800,sum(npp.tot.ps, na.rm=T)]*1E-3*(1/2))/(biom.dat[aoi>800,sum(aoi, na.rm=T)]*1E-4) ## ie 0.82 tC/ha/yr using static/pseudoreps
(biom.dat[aoi>800,sum(npp.tot.mod, na.rm=T)]*1E-3*(1/2))/(biom.dat[aoi>800,sum(aoi, na.rm=T)]*1E-4) ## ie 0.98 tC/ha/yr using exponential modeled/pseudoreps
(biom.dat[aoi>800,sum(npp.tot.mod.cap, na.rm=T)]*1E-3*(1/2))/(biom.dat[aoi>800,sum(aoi, na.rm=T)]*1E-4) ## ie 0.88 tC/ha/yr using exponential capped modeled/pseudoreps
(biom.dat[aoi>800,sum(npp.tot.mod.lin, na.rm=T)]*1E-3*(1/2))/(biom.dat[aoi>800,sum(aoi, na.rm=T)]*1E-4) ## ie 0.86 tC/ha/yr using linear modeled/pseudoreps


hist(biom.dat[aoi>800 & biom>10, npp.tot.ps]) ## up to 1000 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.ps.MgC.ha]) ## up to 6 MgC/ha/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod]) ## up to 800 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.MgC.ha]) ## up to 4 MgC/ha/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.cap]) ## up to 800 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.cap.MgC.ha]) ## up to 4 MgC/ha/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.lin]) ## up to 800 kg/cell/yr
hist(biom.dat[aoi>800 & biom>10, npp.tot.mod.lin.MgC.ha]) ## up to 4 MgC/ha/yr


## so that's interesting. A lot produce not much, but there's a longer tail of high productivity that starts to look more like real forest -- up to ~8 MgC/ha/yr
write.csv(biom.dat, "processed/andy.forest.results.csv")

## an ancillary question: why does FIA have such a pessimistic idea of productivity, maxes out at 2.5 MgC/ha/yr
## probably because those forests evidently grow much slower








#########
####### supplemental investigation: noise testing
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

