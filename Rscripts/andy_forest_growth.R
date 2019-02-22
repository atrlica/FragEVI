library(raster)
library(data.table)


##### ANALYSIS OF TREE CORE RECORDS
## data processing
##### 
andy.bai <- as.data.table(read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv")) 

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
biom.hist <- data.table() ## the running biomass + dbh history of each tree, according to the cores
## figure out biomass growth history applying species-specific allometrics
for(sp in 1:3){
  tmp <- andy.bai[Spp %in% taxa[[sp]],]
  a <- tmp[, lapply(tmp[,7:33], function(x){biom.pred(x, b0.l[sp], b1.l[sp])})]
  names(a) <- paste0("biom", 2016:1990)
  b <- tmp[, 7:33, with=F]
  c <- cbind(tmp[,.(Tree.ID, incr.ID, Spp, X, Y)], a, b)
  biom.hist <- rbind(biom.hist, c)
}
biom.hist$Spp <- as.character(biom.hist$Spp) 

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

### lagged dbh increments, absolute
dbh.incr <- data.frame(c(2015:1990))
for(i in 1:dim(biom.hist)[1]){
  te <- unlist(biom.hist[i,33:59])
  te.di <- (-1*diff(te, lag=1))
  dbh.incr[,i+1] <- te.di
}
names(dbh.incr)[-1] <- paste0("Tree", biom.hist[,incr.ID]) 
names(dbh.incr)[1] <- "Year" ## this has all the relative forward biomass gains by year (row) for each tree with increment (column)

## exploratory plots of dbh increments through time
par(mfrow=c(2,2), mar=c(1,2,2,1), oma=c(2,1,1,1))
plot(dbh.incr$Year, dbh.incr$Tree2)
plot(dbh.incr$Year, dbh.incr$Tree4)
plot(dbh.incr$Year, dbh.incr$Tree14)
plot(dbh.incr$Year, dbh.incr$Tree188)


#####
# ### initial approach to getting a consistent growth~dbh relationship: 
# ### take mean of rel. biomass gain and mean of dbh across all the years present in each core sample
# ### figure the mean relative biomass gain per year across all years in each tree
# mean.na <- function(x){mean(x, na.rm=T)}
# m <- apply(biom.rel[,-1], FUN=mean.na, 2)
# growth.mean <- data.frame(cbind(c(sub("Tree", '', colnames(biom.rel[,-1]))), m))
# names(growth.mean) <- c("incr.ID", "growth.mean")
# growth.mean$incr.ID <- as.integer(as.character(growth.mean$incr.ID))
# andy.bai <- merge(andy.bai, growth.mean, by="incr.ID")
# andy.bai$growth.mean <- as.numeric(as.character(andy.bai$growth.mean))
# andy.bai[,avg.dbh:=apply(as.matrix(andy.bai[,8:34]), FUN=mean.na, 1)]
# 
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

# ### upgrade Jul 23: treat five-year time slices as pseudo-replicates
# pseudo.A <- matrix(ncol=4, rep(999,4))
# pseudo.B <- matrix(ncol=4, rep(999,4))
# pseudo.C <- matrix(ncol=4, rep(999,4))
# pseudo.D <- matrix(ncol=4, rep(999,4))
# pseudo.E <- matrix(ncol=4, rep(999,4)) ## so we have growth histories going back up to 25 years
# 
# ### for each tree core, calculate average annualized (forward) relative growth in successive 5-year chunks, moving chunks backward until you run out of rings
# for(d in unique(andy.bai$incr.ID)){  ### the intervals below are non-overlapping: each 5 year chunk is a distinct biomass time series moving backward
#   pseudo.A <- rbind(pseudo.A, c(((biom.hist[incr.ID==d, biom2016]-biom.hist[incr.ID==d, biom2012])/biom.hist[incr.ID==d, biom2012])/5,
#                                 andy.bai[incr.ID==d, dbh2012],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "A"))
#   pseudo.B <- rbind(pseudo.B, c(((biom.hist[incr.ID==d, biom2011]-biom.hist[incr.ID==d, biom2007])/biom.hist[incr.ID==d, biom2007])/5,
#                                 andy.bai[incr.ID==d, dbh2007],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "B"))
#   pseudo.C <- rbind(pseudo.C, c(((biom.hist[incr.ID==d, biom2006]-biom.hist[incr.ID==d, biom2002])/biom.hist[incr.ID==d, biom2002])/5,
#                                 andy.bai[incr.ID==d, dbh2002],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "C"))
#   pseudo.D <- rbind(pseudo.D, c(((biom.hist[incr.ID==d, biom2001]-biom.hist[incr.ID==d, biom1997])/biom.hist[incr.ID==d, biom1997])/5,
#                                 andy.bai[incr.ID==d, dbh1997],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "D"))
#   
#   ### pseudo.E can have from 5-6 years in a successsful sample
#   tmp <- unlist(biom.hist[incr.ID==d, 26:32, with=F])
#   t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
#   if(t<5 | !is.finite(t)){ ## if record ends 1992-1996 (i.e. does not give at least 5 year run)
#     pseudo.E <- rbind(pseudo.E, c(999,999,
#                                   andy.bai[incr.ID==d, incr.ID],
#                                   "E"))
#     print("not enough data for 96-90 chunk") ## this record not long enough
#   }else{
#     g <- ((biom.hist[incr.ID==d, biom1996]-biom.hist[incr.ID==d, (26+(t-1)), with=F])/biom.hist[incr.ID==d, (26+(t-1)), with=F])/sum(is.finite(tmp))
#     deeb <- min(unlist(andy.bai[incr.ID==d, 28:34, with=F]), na.rm=T)
#     pseudo.E <- rbind(pseudo.E, c(g, deeb,
#                                   andy.bai[incr.ID==d, incr.ID],
#                                   "E"))
#     print("retreived 96-90 chunk")
#   }
# }
# 
# ## collate
# ps.contain <- rbind(pseudo.A[-1,], pseudo.B[-1,], pseudo.C[-1,], pseudo.D[-1,], pseudo.E[-1,])
# ps.contain <- data.frame(biom.rel.ann=as.numeric(unlist(ps.contain[,1])),
#                          dbh.start=as.numeric(unlist(ps.contain[,2])),
#                          incr.ID=as.integer(unlist(ps.contain[,3])),
#                          interval=as.character(unlist(ps.contain[,4])),
#                          stringsAsFactors = F)
# ps.contain$interval <- as.factor(ps.contain$interval)
# ps.contain <- as.data.table(ps.contain)
# ps.contain <- ps.contain[order(incr.ID),]
# ps.contain <- ps.contain[dbh.start!=999,]
# ps.contain <- merge(x=ps.contain, y=andy.bai[,.(Plot.ID, incr.ID, Spp, seg, seg.Edge)], by="incr.ID", all.x=T, all.y=F)


#### upgrade 16 Oct, also check out annualized dbh increment -- and change the intervals
### upgrade Jul 23: treat five-year time slices as pseudo-replicates
pseudo.A <- matrix(ncol=4, rep(999,4))
pseudo.B <- matrix(ncol=4, rep(999,4))
pseudo.C <- matrix(ncol=4, rep(999,4))
pseudo.D <- matrix(ncol=4, rep(999,4))
pseudo.E <- matrix(ncol=4, rep(999,4)) ## so we have growth histories going back up to 25 years

### for each tree core, calculate average annualized (forward) relative growth in successive 5-year chunks, moving chunks backward until you run out of rings
for(d in unique(andy.bai$incr.ID)){  ### the intervals below are non-overlapping: each 5 year chunk is a distinct biomass time series moving backward
  pseudo.A <- rbind(pseudo.A, c((biom.hist[incr.ID==d, dbh2016]-biom.hist[incr.ID==d, dbh2011])/5,
                                andy.bai[incr.ID==d, dbh2011],
                                andy.bai[incr.ID==d, incr.ID],
                                "A"))
  pseudo.B <- rbind(pseudo.B, c((biom.hist[incr.ID==d, dbh2011]-biom.hist[incr.ID==d, dbh2006])/5,
                                andy.bai[incr.ID==d, dbh2006],
                                andy.bai[incr.ID==d, incr.ID],
                                "B"))
  pseudo.C <- rbind(pseudo.C, c((biom.hist[incr.ID==d, dbh2006]-biom.hist[incr.ID==d, dbh2001])/5,
                                andy.bai[incr.ID==d, dbh2001],
                                andy.bai[incr.ID==d, incr.ID],
                                "C"))
  pseudo.D <- rbind(pseudo.D, c((biom.hist[incr.ID==d, dbh2001]-biom.hist[incr.ID==d, dbh1996])/5,
                                andy.bai[incr.ID==d, dbh1996],
                                andy.bai[incr.ID==d, incr.ID],
                                "D"))
  
  ### pseudo.E can have from 5-6 years in a successsful sample
  tmp <- unlist(biom.hist[incr.ID==d, 53:59, with=F])
  t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
  if(t<5 | !is.finite(t)){ ## if record ends 1992-1996 (i.e. does not give at least 5 year run)
    pseudo.E <- rbind(pseudo.E, c(999,999,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    print("not enough data for 96-90 chunk") ## this record not long enough
  }else{
    g <- (biom.hist[incr.ID==d, dbh1996]-biom.hist[incr.ID==d, (53+(t-1)), with=F])/(t-1)
    deeb <- min(unlist(biom.hist[incr.ID==d, 53:59, with=F]), na.rm=T)
    pseudo.E <- rbind(pseudo.E, c(g, deeb,
                                  andy.bai[incr.ID==d, incr.ID],
                                  "E"))
    print("retreived 96-90 chunk")
  }
}

## collate
ps.dbh.contain <- rbind(pseudo.A[-1,], pseudo.B[-1,], pseudo.C[-1,], pseudo.D[-1,], pseudo.E[-1,])
# names(ps.dbh.contain) <- c("dbh.incr", "dbh.start", "incr.ID", "interval")
ps.dbh.contain <- data.frame(dbh.incr.ann=as.numeric(unlist(ps.dbh.contain[,1])),
                             dbh.start=as.numeric(unlist(ps.dbh.contain[,2])),
                             incr.ID=as.integer(unlist(ps.dbh.contain[,3])),
                             interval=as.character(unlist(ps.dbh.contain[,4])),
                             stringsAsFactors = F)
ps.dbh.contain$interval <- as.factor(ps.dbh.contain$interval)
ps.dbh.contain <- as.data.table(ps.dbh.contain)
ps.dbh.contain <- ps.dbh.contain[order(incr.ID),]
ps.dbh.contain <- merge(x=ps.dbh.contain, y=andy.bai[,.(incr.ID, seg, seg.Edge, Spp, Plot.ID)], by="incr.ID", all.x=T, all.y=F)
ps.dbh.contain <- ps.dbh.contain[dbh.start!=999 & !is.na(dbh.start),]

# ps.final <- merge(ps.contain, ps.dbh.contain[, c("dbh.incr.ann", "dbh.start", "incr.ID", "interval")],
#                   by=c("incr.ID", "interval"), all.x=T, all.y=T)
# names(ps.final)[c(4,10)] <- c("dbh.start.biom", "dbh.start.incr")

### dump the pseudo-replicated file to disk
write.csv(ps.dbh.contain, "processed/andy.bai.ps.dbhincr.csv")
#####

### STEM GROWTH~DBH ANALYSIS
#####
library(raster)
library(data.table)
ps.contain <- as.data.table(read.csv("processed/andy.bai.ps.dbhincr.csv")) ## this is the nicely formatted BAI data from the psueoreplicated tree cores
ps.contain[dbh.start>=5, .(quantile(dbh.start, probs=c(0.05, 0.5, 0.95), na.rm=T),
                           quantile(dbh.incr.ann, probs=c(0.05, 0.5, 0.95), na.rm=T)),
           by=seg.Edge]


library(lme4)

## edge slope, random slopes on plot increment interval
# library(lme4)
# d <- lmer(dbh.incr.ann~dbh.start*seg.Edge + 
#             (1+dbh.start|Plot.ID) + 
#             (1+dbh.start|incr.ID) + 
#             (1+dbh.start|interval), 
#           data=ps.contain[dbh.start>=5,],
#           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# summary(d) 
## specifying random effects as 1+slope|random vs slope|random appears to have no effect
# d.alt <- lmer(dbh.incr.ann~dbh.start*seg.Edge + 
#             (dbh.start|Plot.ID) + 
#             (dbh.start|incr.ID) + 
#             (dbh.start|interval), 
#           data=ps.contain[dbh.start>=5,],
#           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# summary(d.alt) 
# coef(summary(d.alt))
# coef(summary(d))

## same random slopes/intercepts, simple effect for edge
# e <- lmer(dbh.incr.ann~dbh.start+seg.Edge + 
#             (dbh.start|Plot.ID) + 
#             (dbh.start|incr.ID) + 
#             (dbh.start|interval), 
#           data=ps.contain[dbh.start>=5,],
#           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# summary(e) ## definitely no DBH effect
# plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
#      col=ps.contain[dbh.start>=5, seg.Edge])
# points(ps.contain[dbh.start>=5, dbh.start], predict(e),
#        col="blue", pch=16)
# plot(e)
# hist(resid(e)) ## looks pretty fine
# qqnorm(resid(e)) ## good enough
# qqline(resid(e))

### remove dbh.start as a factor for dbh.incr
# e.null <- lmer(dbh.incr.ann~seg.Edge + 
#             (dbh.start|Plot.ID) + 
#             (dbh.start|incr.ID) + 
#             (dbh.start|interval), 
#           data=ps.contain[dbh.start>=5,],
#           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# 
# summary(e.null) 
# anova(e.null, e, test="Chisq") ## no dbh effect
# sigma(e.null); sigma(e) ## pretty much the same
# e.null.null <- lmer(dbh.incr.ann~1+ 
#                                 (1+dbh.start|Plot.ID) + 
#                                 (1+dbh.start|incr.ID) + 
#                                 (1+dbh.start|interval), 
#                               data=ps.contain[dbh.start>=5,],
#                               REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# summary(e.null.null) 
# anova(e.null.null, e.null, test="Chisq") ## definitely an edge effect


# ## what does it mean if there's no dbh effect on increment?
# dbh0 <- seq(5, 90)
# dbh1 <- dbh0+0.5 ## say a tree just grows half an inch per year, like model predicts for edge trees
# b0 <- -2.48
# b1 <- 2.4835 ## these are eastern hardwood defaults, Jenkins et al., 2003 (following approach of Smith et al. 2018)
# biom.pred <- function(x){exp(b0+(b1*log(x)))}
# biom.inv <- function(x){exp((log(x)-b0)/b1)}
# biom0 <- biom.pred(dbh0)
# biom1 <- biom.pred(dbh1)
# plot(dbh0, (biom1-biom0)/biom0) ## Fan fookin tastic, here's our good ol hyperbola and based on a static dbh change

### what if we are taking too much out of the dbh effect by modeling it separately for everyone?
# f <- lmer(dbh.incr.ann~dbh.start+seg.Edge +
#             (1|Plot.ID) +
#             (1|incr.ID) +
#             (1|interval),
#           data=ps.contain[dbh.start>=5,],
#           REML=F)  ## even here, not a lot going on with dbh, probably not significant
# summary(f)
# plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
#      col=ps.contain[dbh.start>=5, seg.Edge])
# points(ps.contain[dbh.start>=5, dbh.start], predict(f),
#        col=ps.contain[dbh.start>=5, seg.Edge], pch=16)
# 
# plot(f)
# hist(resid(f)) ## looks pretty fine
# qqnorm(resid(f)) ## good enough
# qqline(resid(f))
# 
# f.null <- lmer(dbh.incr.ann~seg.Edge +
#                  (1|Plot.ID) +
#                  (1|incr.ID) +
#                  (1|interval),
#                data=ps.contain[dbh.start>=5,],
#                REML=F)
# summary(f.null)
# anova(f, f.null, test="Chisq") ## nope, p>0.26
# f.null.null <- lmer(dbh.incr.ann~1 +
#                  (1|Plot.ID) +
#                  (1|incr.ID) +
#                  (1|interval),
#                data=ps.contain[dbh.start>=5,],
#                REML=F)
# summary(f.null.null)
# anova(f.null, f.null.null, test="Chisq") ## yessir, p<0.001
# coef(summary(e.null)) ## random slopes and intercepts
# coef(summary(f.null)) ## random intercepts only
## both models are v. similar in fixed effects

### consider a model without random effect for year interval (same intervals correspond across all observations)
g <- lmer(dbh.incr.ann~seg.Edge+dbh.start + 
                 (dbh.start|Plot.ID) + 
                 (dbh.start|incr.ID), 
               data=ps.contain[dbh.start>=5,],
               REML=F) 
g.full <- lmer(dbh.incr.ann~seg.Edge*dbh.start + 
                 (dbh.start|Plot.ID) + 
                 (dbh.start|incr.ID), 
               data=ps.contain[dbh.start>=5,],
               REML=F) 
summary(g.full)
summary(g) ## now looks like we have a dbh effect (negative slope, not huge)
anova(g, g.full, test="Chisq") ## the interactive term matters!
plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
     col=ps.contain[dbh.start>=5, seg.Edge])
points(ps.contain[dbh.start>=5, dbh.start], predict(g),
       col=ps.contain[dbh.start>=5, seg.Edge], pch=16)

g.lin <- lmer(dbh.incr.ann~dbh.start + seg.Edge:dbh.start+
                          (dbh.start|Plot.ID) + 
                          (dbh.start|incr.ID), 
                        data=ps.contain[dbh.start>=5,],
                        REML=F) 
anova(g.lin, g.full, test="Chisq")

g.seg <- lmer(dbh.incr.ann~seg.Edge+
                (dbh.start|Plot.ID) + 
                (dbh.start|incr.ID), 
              data=ps.contain[dbh.start>=5,],
              REML=F) 
anova(g.seg, g.full, test="Chisq")


g.null <- lmer(dbh.incr.ann~seg.Edge+
                 (dbh.start|Plot.ID) + 
                 (dbh.start|incr.ID), 
               data=ps.contain[dbh.start>=5,],
               REML=F) 
anova(g.null, g, test="Chisq") ## sig dbh effect

g.null.null <- lmer(dbh.incr.ann~1+
                                (dbh.start|Plot.ID) + 
                                (dbh.start|incr.ID), 
                              data=ps.contain[dbh.start>=5,],
                              REML=F) 
anova(g.null.null, g.null, test="Chisq") ### sig edge effect

### save out the model for stem growth for use in results
save(g.full, file = "processed/mod.andy.final.sav")

## OK g.full is now your operative andy model

###
### DEVELOPING MODEL COEFFICIENTS+ERROR FOR PLOT-LEVEL GROWTH RATE
#####
### Push the mixed effects model of dbh increment into the areal-basis model(s) for growth
# load("processed/mod.andy.final.sav")
# b0.rand <- rnorm(1000, mean=coef(summary(g.full))[1,1], sd=coef(summary(g.full))[1,2]) ## intercept
# b1.rand <- rnorm(1000, mean=coef(summary(g.full))[2,1], sd=coef(summary(g.full))[2,2]) ## Interior
# b2.rand <- rnorm(1000, mean=coef(summary(g.full))[3,1], sd=coef(summary(g.full))[3,2]) ## dbh slope
# b3.rand <- rnorm(1000, mean=coef(summary(g.full))[4,1], sd=coef(summary(g.full))[4,2]) ## dbh:Interior
# plot(ddd, mean(b0.rand)+ddd*mean(b2.rand), col="red")
# points(ddd, mean(b0.rand)+mean(b1.rand)+ddd*(mean(b2.rand)+mean(b3.rand)), col="blue")

### constrain the predictions of dbh increment to the observed increment range
# ps.contain <- as.data.table(read.csv("processed/andy.bai.ps.dbhincr.csv")) ## this is the nicely formatted BAI data from the psueoreplicated tree cores
# dbh.incr.min.edge <- ps.contain[seg.Edge=="E" & dbh.start>5, min(dbh.incr.ann, na.rm=T)]
# dbh.incr.max.edge <- ps.contain[seg.Edge=="E" & dbh.start>5, max(dbh.incr.ann, na.rm=T)]
# dbh.incr.min.int <- ps.contain[seg.Edge=="I" & dbh.start>5, min(dbh.incr.ann, na.rm=T)]
# dbh.incr.max.int <- ps.contain[seg.Edge=="I" & dbh.start>5, max(dbh.incr.ann, na.rm=T)]

## hunting for why the recent runs of the andy forests are so different from the ones I got prior to CO2USA
# b0.rand <- rnorm(100, mean=coef(summary(f.null))[1,1], sd=coef(summary(f.null))[1,2])
# b1.rand <- rnorm(100, mean=coef(summary(f.null))[2,1], sd=coef(summary(f.null))[2,2])
# par(mfrow=c(1,2))
# hist(b0.rand.f)
# hist(b0.rand)
# hist(b1.rand.f)
# hist(b1.rand) ### there really isn't a big difference between the two ME models here (f.null was mid-october, switched to)

## estimate future biomass in each stem for many realizations of the fuzzy dbh-change model
# dump <- copy(andy.dbh)
# dump.plot <- dump[, sum(biom0), by=.(Plot.ID, seg)]
# names(dump.plot)[3] <- "Init.biomass.kg"
# for(i in 1:1000){
#   ### predict next year'd DBH and restrict 
#   dump[, dbh.delt.pred:=9999]
#   dump[seg.F=="E", dbh.delt.pred:=b0.rand[i]+(dbh*b2.rand[i])]
#   dump[seg.F=="E" & dbh.delt.pred>dbh.incr.max.edge, dbh.delt.pred:=dbh.incr.max.edge]
#   dump[seg.F=="E" & dbh.delt.pred<dbh.incr.min.edge, dbh.delt.pred:=dbh.incr.min.edge]
#   
#   dump[seg.F=="I", dbh.delt.pred:=b0.rand[i]+b1.rand[i]+(dbh*(b2.rand[i]+b3.rand[i]))]
#   dump[seg.F=="I" & dbh.delt.pred>dbh.incr.max.int, dbh.delt.pred:=dbh.incr.max.int]
#   dump[seg.F=="I" & dbh.delt.pred<dbh.incr.min.int, dbh.delt.pred:=dbh.incr.min.int]
#   
#   ## apply next year's dbh to get next year's biomass by SPP
#   for(b in 1:dim(biom.pred.key)[1]){
#     dump[Spp==biom.pred.key[b, Spp], biom.t:=exp(biom.pred.key[b, b0]+(biom.pred.key[b, b1]*log(dbh.delt.pred+dbh)))]
#   }
#   dump[,growth:=biom.t-biom0] ## record the annual biomass change per stem
#   andy.dbh <- cbind(andy.dbh, dump$growth)
#   names(andy.dbh)[11+i] <- paste0("growth.iter", i, ".kg") ## rename col
#   print(paste("plot growth iteration", i))
#   nerd <- dump[,sum(growth), by=.(Plot.ID, seg)]
#   dump.plot <- merge(dump.plot, nerd, by=c("Plot.ID", "seg"))
#   names(dump.plot)[i+3] <- paste("plot.growth.iter", i, ".kg", sep="")
# }
# 
# write.csv(andy.dbh, file="processed/andy.dbh.growth.iter.csv")
# write.csv(dump.plot, file="processed/andy.plot.growth.iter.csv")
### up to here, totally same October vs. November


###
### now develop your plot model coefficients from random draws/projection of your stem growth model
andy.dbh <- read.csv("processed/andy.dbh.growth.iter.csv")
andy.dbh <- andy.dbh[,-1]
andy.dbh <- as.data.table(andy.dbh)
### now estimate coefficients for the area-basis model using the invariant present-day biomass and varying estimates of future biomass
## --> ends up being a simple linear model of MgC-growth/MgC-biomass vs. MgC/ha-can with a single correction factor for edge/interior
plot.mod.b0 <- numeric()
plot.mod.b1 <- numeric()
plot.mod.b2 <- numeric()
edge.max <- numeric()
edge.min <- numeric()
int.max <- numeric()
int.min <- numeric()

## this loop uses the prepared estimates of per-stem growth to iteratively estimate model of plot-basis biomass growth vs. biomass denisty
for(i in 1:1000){ 
  # ### original formulation
  # yip <- andy.dbh[,.(Tree.ID, Plot.ID, Spp, X, Y, dbh, Can.class, seg, seg.F, biom0)]
  # target <- names(andy.dbh)[11+i]
  # yip <- cbind(yip, as.numeric(andy.dbh[,get(target)]))
  # k <- yip[, .(sum(biom0), sum(V2)), by=.(seg, Plot.ID)]
  # k[,V3:=V2/V1]
  # k[,V1:=((V1/2000)/(10*30))*1E4] ## MgC/ha
  # k[,V2:=V2/2000] ## MgC
  # plot(k$V1, k$V2, col=k$seg, main=paste("iter", i, "absolute growth"))
  # plot(k$V1, k$V3, col=k$seg, main=paste("iter", i, "relative growth"))
  # k[,seg.F:="I"]
  # k[seg==10, seg.F:="E"]
  # # aa <- lm(V2~V1+seg.F, data=k)
  # aa <- lm(V3~V1+seg.F, data=k) ## MgC-growth/MgC-biomass(MgC)~density(MgC/ha)+interior
  # points(k$V1, predict(aa), col="black", cex=0.4, pch=16)
  # plot.mod.b0 <- c(plot.mod.b0, coef(aa)[1])
  # plot.mod.b1 <- c(plot.mod.b1, coef(aa)[2])
  # plot.mod.b2 <- c(plot.mod.b2, coef(aa)[3])
  # ### first finding: this appears to be just like the below, only shittier code

  ### November formulation
  g <- andy.dbh[, .(sum(biom0), sum(get(names(andy.dbh)[11+i]))), by=.(seg, Plot.ID)]
  g[,growth.rel:=V2/V1]
  g[,MgC.ha.can:=((V1/2000)/(10*30))*1E4] ## MgC/ha
  g[,MgC.growth:=V2/2000] ## MgC growth
  g[,seg.F:="I"]
  g[seg==10, seg.F:="E"]
  aa <- lm(growth.rel~MgC.ha.can+seg.F, data=g) ## MgC-growth/MgC-biomass(MgC)~density(MgC/ha)+interior
  plot.mod.b0 <- c(plot.mod.b0, rnorm(n=1, mean=summary(aa)$coefficients[1,1], sd=summary(aa)$coefficients[1,2]))
  plot.mod.b1 <- c(plot.mod.b1, rnorm(n=1, mean=summary(aa)$coefficients[2,1], sd=summary(aa)$coefficients[2,2]))
  plot.mod.b2 <- c(plot.mod.b2, rnorm(n=1, mean=summary(aa)$coefficients[3,1], sd=summary(aa)$coefficients[3,2]))
  # plot(g[,MgC.ha.can], g[,growth.rel], col=as.numeric(g[,seg]))
 
   ### need some basic stats to constrain estimates in the npp calculation step
  edge.max <- c(edge.max, g[seg.F=="E", max(growth.rel)])
  edge.min <- c(edge.min, g[seg.F=="E", min(growth.rel)])
  int.max <- c(int.max, g[seg.F=="I", max(growth.rel)])
  int.min <- c(int.min, g[seg.F=="I", min(growth.rel)])
  print(paste("plot model iteration", i))
}

# 
# ## intercept 5.7%, interior -1.6%, -2.6% per 100 MgC+
# t.test(plot.mod.b1, mu=0)
# t.test(plot.mod.b2, mu=0)
## the whole purpose of the above is to produce a vector of model coefficients for use below in an interative NPP estimate
#####

### NPP Calculation treating canopy as Andy-like forest, with model error
library(raster)
library(data.table)
### load up 30m data
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
can <- raster("processed/boston/bos.can.redux30m.tif")
ed.can <- raster("processed/boston/bos.ed10m.redux30m.tif")
ed.biom <- raster("processed/boston/bos.biomass.ed10only30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
ed.can <- crop(ed.can, aoi)
ed.biom <- crop(ed.biom, aoi)

biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=getValues(aoi)]
biom.dat[,can:=getValues(can)]
biom.dat[,ed.can:=getValues(ed.can)]
biom.dat[,ed.biom:=getValues(ed.biom)]
biom.dat[,isa:=getValues(isa)]
biom.dat[,lulc:=getValues(lulc)]
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]

names(biom.dat)[1] <- c("biom")

dim(biom.dat[biom>10 & !is.na(biom) & aoi>800,]) ## 106659 valid biomass pixels to do

## Figure out working parameters in the map data and prep for NPP calc
biom.dat[,int.can:=can-ed.can]
# biom.dat[,range(int.can, na.rm=T)] ## includes negative numbers
# biom.dat[int.can<0 & aoi>800, length(int.can)] ## 11 complete pixels with bad int can coverage
## kill the artifacts
biom.dat[biom==0, int.can:=0]
biom.dat[biom==0, ed.can:=0]
biom.dat[is.na(biom), ed.can:=NA]
biom.dat[is.na(biom), int.can:=NA]
biom.dat[int.can<0 & int.can>(-0.01), int.can:=0] ## kill rounding errors
# biom.dat[int.can<0, length(int.can)] ## 70, only 2 in complete pixels

## now ID biomass by edge vs int
biom.dat[,int.biom:=biom-ed.biom] # internal forest biomass
# biom.dat[,range(int.biom, na.rm=T)] ## 0-51k
# biom.dat[,range(ed.biom, na.rm=T)] ## 0-31k (so the peak biomass cells are deep forest somewhere)
# hist(biom.dat[aoi>800, int.biom])
# hist(biom.dat[aoi>800, ed.biom])

### range of edge and interior biomass, in Mg-biomass/ha
biom.dat[aoi>800, ed.biom.MgC:=(ed.biom/2000)]
biom.dat[aoi>800, ed.biom.MgC.ha.can:=(ed.biom.MgC/(aoi*ed.can))*1E4] ### the model injests density MgC/ha-can and spits out growth factor MgC-growth/MgC-biomass
biom.dat[aoi>800 & ed.can<0.005, ed.biom.MgC.ha.can:=0]
# hist(biom.dat[aoi>800 & ed.biom.MgC.ha.can>0, ed.biom.MgC.ha.can])
# biom.dat[aoi>800, quantile(ed.biom.MgC.ha.can, probs=c(0.05, 0.95), na.rm=T)] ## 0 to 147 MgC/ha
# biom.dat[aoi>800, range(ed.biom.MgC.ha.can, na.rm=T)]
# biom.dat[aoi>800 & ed.biom.MgC.ha.can<=100 & ed.biom.MgC.ha.can>5 & !is.na(ed.biom.MgC.ha.can), length(ed.biom.MgC.ha.can)]/biom.dat[aoi>800 & !is.na(ed.biom.MgC.ha.can), length(ed.biom.MgC.ha.can)]
## 49% below 100 MgC/ha-can

# biom.dat[aoi>800, range(int.biom, na.rm=T)] ## in kg/cell
biom.dat[aoi>800, int.biom.MgC:=(int.biom/2000)]
biom.dat[aoi>800, int.biom.MgC.ha.can:=(int.biom.MgC/(aoi*int.can))*1E4]
biom.dat[aoi>800 & int.can<0.005, int.biom.MgC.ha.can:=0]
hist(biom.dat[aoi>800 & int.biom.MgC.ha.can>0, int.biom.MgC.ha.can])
biom.dat[aoi>800, quantile(int.biom.MgC.ha.can, probs=c(0.05, 0.95), na.rm=T)] ## 0 to 160 MgC/ha, higher peak
biom.dat[aoi>800, range(int.biom.MgC.ha.can, na.rm=T)]
biom.dat[aoi>800 & int.biom.MgC.ha.can<=100 & int.biom.MgC.ha.can>1 & !is.na(int.biom.MgC.ha.can), length(int.biom.MgC.ha.can)]/biom.dat[aoi>800 & !is.na(int.biom.MgC.ha.can) & int.biom.MgC.ha.can>0, length(int.biom.MgC.ha.can)]
# ## of pixels with interior biomass, 11% is below 100 MgC/ha


###
### VERSION 2: ITERATE ESTIMATION USING ERROR DISTRIBUTION OF COEFFICIENTS IN (LINEAR) PLOT-LEVEL MODEL
### VERSION 3: same as version 2 but using random slopes+intercepts model for dbh increment (edge segment but no dbh.start effect) (e.null)
### recall: this model predicts growth factor (MgC-growth/MgC-biomass) a f(MgC/ha) in closed canopy, with correction factor for interior canopy
### VERSION 4: Version 3, but uses a somewhat simpler g.full model for stem growth rate
### VERSION 5: Version 4, but iterates 1000 times AND uses the biomass>0 canopy map and edge can/biomass map

## establish limits for estimated growth factors based on the distribution of estimated productivities seen in the plots
edge.hi <- mean(edge.max)+sd(edge.max) ## these were the records we got from the iterative plot biomasses obtained above
edge.lo <- mean(edge.min)-sd(edge.min)
int.hi <- mean(int.max)+sd(int.max)
int.lo <- mean(int.min)-sd(int.min)

### these are the distributions of the coefficients for the plot growth model (intercept different for interior)
b0.sel <- rnorm(n = 1000, mean=mean(plot.mod.b0), sd=sd(plot.mod.b0))
b1.sel <- rnorm(n = 1000, mean=mean(plot.mod.b1), sd=sd(plot.mod.b1))
b2.sel <- rnorm(n = 1000, mean=mean(plot.mod.b2), sd=sd(plot.mod.b2))

med.gfact.edge <- numeric()
med.gfact.int <- numeric() ## record the median growth factors for each realization
pixnum.gfact.edge.max <- numeric()
pixnum.gfact.edge.min <- numeric()
pixnum.gfact.int.max <- numeric()
pixnum.gfact.int.min <- numeric()
for(i in 1:1000){
  dump <- copy(biom.dat)
  ## October code
  ## assign growth factors based on position and biomass density
  # dump[,edge.fact:=b0.sel+(b1.sel*ed.biom.MgC.ha.can)] ## get a vector of edge factors
  # dump[ed.biom.MgC.ha.can<0.1, edge.fact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  # dump[,int.fact:=b0.sel+b2.sel+(b1.sel*int.biom.MgC.ha.can)] ## apply growth model for interior growth to interior biomass density
  # dump[int.biom.MgC.ha.can<0.1, int.fact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  # dump[int.fact<int.lo, int.fact:=int.lo] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  # dump[edge.fact<edge.lo, edge.fact:=edge.lo] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  # dump[int.fact>int.hi, int.fact:=int.hi] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  # dump[edge.fact>edge.hi, edge.fact:=edge.hi] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  # summary(dump[aoi>800 & edge.fact>0, edge.fact]); hist(dump[aoi>800 & edge.fact>0, edge.fact])
  # summary(dump[aoi>800 & int.fact>0, int.fact]); hist(dump[aoi>800 & int.fact>0, int.fact])
  
  
  ### November code
  ## assign growth factors based on position and biomass density
  dump[,edge.fact:=b0.sel[i]+(b1.sel[i]*ed.biom.MgC.ha.can)] ## get a vector of edge factors
  dump[edge.fact<edge.lo, edge.fact:=edge.lo] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  dump[edge.fact>edge.hi, edge.fact:=edge.hi] ## eliminate unrealistically high growth factors
  dump[ed.biom.MgC.ha.can<0.1, edge.fact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  dump[,int.fact:=b0.sel[i]+b2.sel[i]+(b1.sel[i]*int.biom.MgC.ha.can)] ## apply growth model for interior growth to interior biomass density
  dump[int.fact<int.lo, int.fact:=int.lo] 
  dump[int.fact>int.hi, int.fact:=int.hi] 
  dump[int.biom.MgC.ha.can<0.1, int.fact:=0] 
  # summary(dump[aoi>800 & edge.fact>0, edge.fact]); hist(dump[aoi>800 & edge.fact>0, edge.fact])
  # summary(dump[aoi>800 & int.fact>0, int.fact]); hist(dump[aoi>800 & int.fact>0, int.fact])
  ### compared to the October code, the edge factors are just slightly higher when you calculate this way
  ## the interior factors come up across the board
  
  ## in summary, I can detect no substantive difference in the code from Oct-Nov other than the slight differences
  ## in how the final growth factors are calculated
  
  ## calculate NPP and store
  dump[,npp.tot:=(ed.biom*edge.fact)+(int.biom*int.fact)] ## apply growth factor to biomass present
  names(dump)[17] <- paste0("andy.npp.iter.", i, ".kg")
  
  write.csv(dump[,17], paste0("processed/results/andy/mod/andy.forest.results.iter", i, ".csv"))
  # sav.ground <- cbind(sav.ground, dump[,npp.kg.hw.ground])
  # names(sav.ground)[9+i] <- paste("npp.fia.ground.iter.", i, ".kg", sep="")
  # sav.forest <- cbind(sav.forest, dump[,npp.kg.hw.forest])
  # names(sav.forest)[9+i] <- paste("npp.fia.forest.iter.", i, ".kg", sep="")
  # sav.perv <- cbind(sav.perv, dump[,npp.kg.hw.perv])
  # names(sav.perv)[9+i] <- paste("npp.fia.perv.iter.", i, ".kg", sep="")
  print(paste("iteration",i))
  med.gfact.edge <- c(med.gfact.edge, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, median(edge.fact, na.rm=T)])
  med.gfact.int <- c(med.gfact.int, dump[aoi>800 & int.biom.MgC.ha.can>=0.1, median(int.fact, na.rm=T)])
  pixnum.gfact.edge.max <- c(pixnum.gfact.edge.max, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & edge.fact==edge.hi, length(edge.fact)])
  pixnum.gfact.edge.min <- c(pixnum.gfact.edge.min, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & edge.fact==edge.lo, length(edge.fact)])
  pixnum.gfact.int.max <- c(pixnum.gfact.int.max, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & int.fact==int.hi, length(int.fact)])
  pixnum.gfact.int.min <- c(pixnum.gfact.int.min, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & int.fact==int.lo, length(int.fact)])
}
frac.gfact.edge.max <- pixnum.gfact.edge.max/dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, length(edge.fact)]
frac.gfact.edge.min <- pixnum.gfact.edge.min/dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, length(edge.fact)]
frac.gfact.int.max <- pixnum.gfact.int.max/dump[aoi>800 & int.biom.MgC.ha.can>=0.1, length(int.fact)]
frac.gfact.int.min <- pixnum.gfact.int.min/dump[aoi>800 & int.biom.MgC.ha.can>=0.1, length(int.fact)]
summary(med.gfact.edge)
summary(med.gfact.int) ## median growth factor for interior is minimum set through the median of all 1000 runs
summary(frac.gfact.edge.max); hist(frac.gfact.edge.max) ## almost all 0%; at most 2% get the max factor
summary(frac.gfact.edge.min); hist(frac.gfact.edge.min) ### very few get minimum, but up to 93%
summary(frac.gfact.int.max); hist(frac.gfact.int.max) ## nearly none get to max
summary(frac.gfact.int.min); hist(frac.gfact.int.min) ## like most get the minimum
summary(dump[aoi>800 & ed.biom.MgC.ha.can>0.1, edge.fact]); hist(dump[aoi>800 & ed.biom.MgC.ha.can>0.1, edge.fact]) ## peaking in this realization about 0.05, over by 0.07
summary(dump[aoi>800 & int.biom.MgC.ha.can>0.1, int.fact]); hist(dump[aoi>800 & int.biom.MgC.ha.can>0.1, int.fact]) ## most are stuck at the minimum

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1014)
tmp[,1:14] <- unlist(biom.dat)
for(g in 1:1000){
  aa <- fread(paste0("processed/results/andy/mod/andy.forest.results.iter", g, ".csv"))
  tmp[,g+14] <- unlist(aa[,2])
  print(paste("built iteration", g))
}
fwrite(tmp, "processed/andy.forest.results.V5.csv") ## with biomass>0 canopy, g.full model, no dbh.incr|interval random effect, dbh*seg.E fixed effects, stem growth predictions empirically constrained


### look at tmp -- make sure it's cool
tmp.andy <- fread("processed/andy.forest.results.V5.csv")
hist(tmp.andy[biom>10 & aoi>800, ed.biom.MgC.ha.can]) ## nice bell curve centers about 80 MgC/ha
hist(tmp.andy[biom>10 & aoi>800 & int.biom.MgC.ha.can>30, int.biom.MgC.ha.can]) ### vast majority is like low-density incidental shit, but a minor symmetric peak at 150 also
b0.mn <- 0.056
b1.mn <- -2.52E-04
b2.mn <- -0.017
## predicted uptake rates then, using this shit
ed.gfact.mn <- tmp.andy[biom>10 & aoi>800, (b0.mn+(ed.biom.MgC.ha.can*b1.mn))]
hist(ed.gfact.mn) ## about 0.04 peak
int.gfact.mn <- tmp.andy[biom>10 & aoi>800 & int.biom.MgC.ha.can>30, (b0.mn+(int.biom.MgC.ha.can*b1.mn)+b2.mn)] ## a lot of negative values
hist(int.gfact.mn)

### What is median productivity of andy pixels, edge vs. int?
dim(tmp.andy) ## 354068
dim(tmp.andy[!is.na(biom) & aoi>800 & biom>10,]) ## 106659 with valid biomass and complete
dim(tmp.andy[!is.na(andy.npp.iter.3.kg),]) ## 135705 valid growth retreivals

median.na <- function(x){median(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp.andy[,15:1014]), MARGIN=1, FUN=median.na)
rrr <- raster(biom)
rrr <- setValues(rrr, npp.med)
writeRaster(rrr, filename = "processed/results/andy/bos.andy.forest.V5.npp.tif", format="GTiff", overwrite=T)
sum(npp.med, na.rm=T)/2000 ## 8.05 ktC



##
### in summary, I can't replicate the results that were produced in October re Andy Forest
### code review from the Git commits doesn't even produce good looking results
### I suspect it's something to do with either the data that was input (though a quick inspection
### of the associated data files appear identical), or has something to do with
### how growth factors were distributed over the landscape based on the modeled growth rates
### I am inclined to believe the V3 results as I have gone over their workings very carefully
### and they do not appear to have any errors and get the basic idea right
# par(mfrow=c(1,2))
# hist(dump[edge.fact>0, edge.fact])## as an example, higher factors overall
# hist(dump[int.fact>0, int.fact])## as an example, higher factors overall
# sum.na <- function(x){sum(x, na.rm=T)}
# hur <- apply(biom.dat[,15:114], MARGIN=2, FUN=sum.na)
# hist(hur/2000) ## ca 10k tC in V4 model, 12.5k tC in the V3 model
# median(hur/2000)
## an ancillary question: why does FIA have such a pessimistic idea of productivity, maxes out at 2.5 MgC/ha/yr
## probably because those forests evidently grow much slower

#####