### respiration 

### Annual CO2 efflux from Decina 2016
rdat <- data.frame(mo=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"),
                   umolCO2.m2.s=c(0,0,0,0,2.25,2.75,3.7,3.12,2.6,1.75,0,0))

rdat$MgC.m2.s <- rdat$umolCO2.m2.s*1E-6*44.01*1E-6*12/44 ## umol/mol, g.CO2/mol.CO2, MgCO2/gCO2, MgC/MgCO2
rdat$MgC.ha.mo <- rdat$MgC.m2.s*1E4*60*60*24*30 ## m2/ha, s/min, min/hr, hr/d, d/mo
plot(rdat$mo, rdat$MgC.ha.mo)
sum(rdat$MgC.ha.mo) ## 5 MgC/ha/yr, not counting anything in spring/winter


## savage and davidson NE study, at HF: 50-400 mgC m-2 hr-1 --> 0.37 to 0.99 kgC m-2 yr-2 depending on year and drainage class

2.25*(1E-6)*44.01*(12/44)*1000*60*60  ### so 2.25 umol/m2/s is 97 mgC m-2 hr-1
3.7*(1E-6)*44.01*(12/44)*1000*60*60 ## gets up to 160 mgC m-2 hr-1

## Davidson 1998 was mean 7.2 (5.3-8.5 depending on drainage) Mg CO2/ha annual, Q10 about 3.9 overall (apparent)
7.2*(12/44) ### 1.96 contrast decina reported growing season total was ~5 MgC/ha

## templer sinusoid curve
## decina linear model based on soil conditions


