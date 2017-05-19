rm(list=ls())
setwd("~/research/miller_et_al/")

#libraries
library(ggplot2)
library(dplyr)

#functions

boxes <- function(Values)
{
  L <- length(Values)
  par(mfrow=c(2,ceiling(L/2)), mai=c(0.3,0.3,0.2,0.2))
  
  for(i in 1:L) {boxplot(Values[,i], main=names(Values)[i])}
}

hists <- function(Values)
{
  L <- length(Values)
  par(mfrow=c(2,ceiling(L/2)), mai=c(0.3,0.3,0.2,0.2))
  
  for(i in 1:L) {hist(Values[,i], breaks = 10, main=names(Values)[i])}
}

# extra functions for pairs multicollinearity exploration
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

# lets go!

jm_data<-read.csv("input_data/dens_abs_1km_extract_correct.csv", h=T)


nrow(jm_data)
length(which(is.na(jm_data$SST)==TRUE))
length(which(is.na(jm_data$CHLA)==TRUE))
jm_data<-na.omit(jm_data) ## remove NAs, wow ok good I did this! 
# na.omit is integral to information theoretic approach1
jm_data[jm_data$BATHY>0,]$BATHY<-0

jm_data$YEAR<-substr(jm_data$St_date, 7,10)
jm_data$YEAR<-factor(jm_data$YEAR)
jm_data$MONTH<-substr(jm_data$St_date, 4,5)
jm_data$MONTH<-factor(jm_data$MONTH)

jm_data$DAYNIGHT<-"DAY"
jm_data[jm_data$Survey=="Birojima5" | jm_data$Survey=="Birojima6",]$DAYNIGHT<-"NIGHT"
jm_data$DAYNIGHT<-factor(jm_data$DAYNIGHT)

jm_data$Location<-factor(jm_data$Location)
jm_data$Survey<-factor(jm_data$Survey)

### Add correct effort offset as per reviewers comments

# just tidy up, there are some points in Survey: Izu3 that don't belong there, we will remove
# DONE ###########################################
#j1<-jm_data[jm_data$Survey=="Izu3",]
#pres_ras<-rasterize(data.frame(j1$Longitude,
#                               j1$Latitude),
#                    rasterTemplate, background=NA, field=1)
#surv_pol<-rasterToPolygons(pres_ras, na.rm=T)
#row.names(j1)<-1:nrow(j1)
#surv_spdf<-SpatialPolygonsDataFrame(surv_pol, data=j1)

#surv<-surveys[surveys$Survey=="Izu3",]
#plot(surv_spdf)
#plot(surv, add=T, col=2)
#surv_gd<-surv_spdf[surv,]
#plot(surv_gd, add=T, border=3)

#j_temp<-jm_data[jm_data$Survey!="Izu3",]
#out<-rbind(j_temp, surv_gd@data)
#write.csv(out, "input_data/dens_abs_1km_extract_correct.csv", quote=F, row.names=F)

##############################################

library(rgdal)
library(raster)
library(rgeos)
surveys<-readOGR(dsn="GIS", layer="all_cols_survey_tracks")
rasterTemplate<-raster("extract_rasters/bathy_1km.tif")


jm_data$perc1km_surv<-0
for(i in 1:nrow(jm_data))
{
d1<-jm_data[i,]  
pres_ras<-rasterize(data.frame(d1$Longitude, d1$Latitude),
                    rasterTemplate, background=NA, field=1)
sp1<-SpatialPoints(data.frame(d1$Longitude, d1$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

pres_ras<-crop(pres_ras, extent(sp1)+0.02)
surv<-surveys[surveys$Survey==d1$Survey,]

DgProj <- CRS(paste("+proj=laea +lon_0=135 +lat_0=34", sep=""))
TshapeProj <- spTransform(surv, CRS=DgProj)
TBuffProj <- gBuffer(TshapeProj, width=100, quadsegs=50)
#plot(TBuffProj)
#plot(TshapeProj, add=T, col="grey")
surv_buff <- spTransform(TBuffProj, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#from ras

onekmsq<-rasterToPolygons(pres_ras, na.rm=T)
#plot(onekmsq)
#plot(surv_buff, add=T)
area_surv<-gIntersection(surv_buff, onekmsq, byid = T)
#plot(area_surv, add=T, border=2)
jm_data[i,]$perc1km_surv<-(area(area_surv)/area(onekmsq))*100
print(i)
}

#jm_data[jm_data$perc1km_surv==999,]$perc1km_surv<-0.5
#jm_data[jm_data$perc1km_surv>100,]$perc1km_surv<-100
#summary(jm_data$perc1km_surv)

#write.csv(jm_data, "input_data/dens_abs_1km_extract_correct.csv", quote=F, row.names=F)



############~~~~~~~~~~ MULTICOLLINEARITY ANALYSES ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########

pairs(data.frame(jm_data$BATHY, jm_data$SLOPE, jm_data$D_COL, jm_data$COL_INF,
                 jm_data$D_LAND, jm_data$D_KURO, jm_data$SST, jm_data$G_SST, jm_data$CHLA,
                 jm_data$G_CHLA, jm_data$YEAR, jm_data$MONTH, jm_data$Density), upper.panel = panel.smooth,lower.panel=panel.cor)

# CHLA and G_CHLA correlated at 0.8 (pearsons)

cMod<-glm(Density~CHLA, data=jm_data, family=poisson(link=log), na.action=na.omit)

gMod<-glm(Density~G_CHLA, data=jm_data, family=poisson(link=log), na.action=na.omit)

summary(cMod)
exp(-0.234174)
exp(-0.003086)

summary(gMod)
exp(-0.245068)
exp(-0.002432)

## will remove G_CHLA from further analyses as it is a product of the original CHLA data

############~~~~~~~~~~  VARIABLE TRANSFORMATIONS  ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########

rawdf<-jm_data

trfdf<-rawdf

# remove outliers
boxes(tr_dens[,7:16])

tr_dens<-tr_dens[tr_dens$G_SST<4,]
tr_dens<-tr_dens[tr_dens$SLOPE<22,]
tr_dens<-tr_dens[tr_dens$CHLA<14,]

# standardize for lme4
library(vegan)
enviro_std<-decostand(tr_dens[,7:16], method="standardize")
tr_dens_std<-data.frame(tr_dens[,1:6], enviro_std, tr_dens[,17:21])


#trfdf$CHLA<-log(trfdf$CHLA)
#trfdf$G_CHLA<-log(trfdf$G_CHLA)
#trfdf$BATHY<-log((trfdf$BATHY+trfdf$BATHY^2)+1)
#trfdf$SLOPE<-sqrt(trfdf$SLOPE)
#trfdf$D_KURO<-sqrt(trfdf$D_KURO)
#trfdf$D_COL<-sqrt(trfdf$D_COL)
#trfdf$D_LAND<-sqrt(trfdf$D_LAND)

hists(trfdf[,6:16])

##

surv_meta<-cbind(aggregate(effort~Survey, jm_data, mean), aggregate(Density~Survey, jm_data, mean), 
                 aggregate(St_date~Survey, jm_data, function(x) (x[1])))

cor(surv_meta$effort, surv_meta$tot_mur) #0.75

qplot(x=surv_meta$effort, y=surv_meta$tot_mur,colour=surv_meta$Survey, geom="point")

############~~~~~~~~~~  Modelling  ~~~~~~~~~~~#########
############~~~~~~~~~~ *********** ~~~~~~~~~~~#########
library(lme4)
library(car)
#library(gamm4)
source("scripts/glmm_helper_funcs.R")

## Night survey data sig different? ##

biro_temp<-jm_data[jm_data$Location=="Birojima",]

levels(biro_temp$Survey)

biro_temp$Survey<-relevel(biro_temp$Survey, ref="Birojima5")

night_mod<-glmer(Density~Survey + (1|YEAR) + (1|MONTH), data=biro_temp, family=poisson)

summary(night_mod)

night_mod<-glm(Density~Survey, data=biro_temp, family=poisson)

summary(night_mod)
sum(resid(night_mod, type='pearson')^2)/df.residual(night_mod)

library(MASS)
night_mod2<-glm.nb(Density~Survey, data=biro_temp)
sum(resid(night_mod2, type='pearson')^2)/df.residual(night_mod2) #1
summary(night_mod2) 
# the night surveys are pres-abs only so wouldnt affect density!


##########******** DENSITY MODEL ********##########
############~~~~~~~************~~~~~~~~~~~#########

tr_dens<-trfdf[-which(trfdf$Survey=="Oki1" | trfdf$Survey=="Birojima5" | trfdf$Survey=="Birojima6"),]

## drop Hachijojima colony as per reviewers comments, no actually keep - its basically part of izu

tr_dens$Location<-factor(tr_dens$Location)
tr_dens$Survey<-factor(tr_dens$Survey) #  drop omitted factor levels
#tr_dens$YEAR<-factor(tr_dens$YEAR) dont need year as Survey RF splits perfectly

# Investigation of outliers in Density data

ggplot(jm_data, aes(x=Density, y=Density))+ geom_point() + facet_wrap(~Survey)

g1<-glmer(Density~1 +
            (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log))
g2<-glmer(Density~1 +
            (1+Density|Location/Survey), data=tr_dens[tr_dens$Density<50,],family=poisson(link=log))
g3<-glmer(Density~1 +
            (1+Density|Location/Survey), data=tr_dens[tr_dens$Density<25,],family=poisson(link=log))

sum(resid(g1, type='pearson')^2)/df.residual(g1) #1.73
sum(resid(g2, type='pearson')^2)/df.residual(g2) #1.36
sum(resid(g3, type='pearson')^2)/df.residual(g3) #1.02

# density cutoff trial

co_seq<-sort(seq(14, 30, 2), decreasing =T)

for ( i in co_seq)
{
  tmod<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                (1|Location/Survey), data=tr_dens[tr_dens$Density<=i,], family=poisson(link=log), na.action=na.omit)
  
  
  qqPlot(resid(tmod))
  print(i)
  print(sum(resid(tmod, type='pearson')^2)/df.residual(tmod))
  print(AIC(tmod))
  
}


### Density set to < 25 !! ###
tr_dens<-tr_dens[tr_dens$Density<25,]

##############################

#par(fig=c(0,1,0,0.35))
boxplot(tr_dens$Density, horizontal=T, bty="n", xlab="JM counts", ylim=c(0,20), cex.lab=1.2)
#par(fig=c(0,1,0.15,1), new=T)
x <- sort(tr_dens$Density)
hist(tr_dens$Density, freq=F, main="", col="darkgray", ylim=c(0,0.9),breaks=seq(0,25,1), xlab="");box()
legend("topright", c("log-normal distribution", "exponential distribution","poisson distribution","nbinom distribution", "density counts"),
       lty=c(1,1,3,1,2), lwd=c(1.5,1,3,2,1), col=c(2, "grey30","darkgreen","black", 4), bty="n")
#distributions
lines(density(tr_dens$Density), lty=2, col=4)
#negative binomial
curve(dnbinom(round(x), size=0.5, mu=1), lty=1,lwd=2, col="black", add=T)
#lognormal
curve(dlnorm(x, meanlog=0, sdlog=1), lwd=1.5, add=T, col=2)
#exponential
curve(dexp(x), add=T, lty=1, col="grey30", lwd=1.5)
#poisson
lines(dpois(x,mean(tr_dens$Density)), lty=3, lwd=3, col="darkgreen")
rug(tr_dens$Density, ticksize = 0.02)


## Go for Density<25 cutoff to make modelling better

table(tr_dens$Survey)
tr_dens$Survey<-relevel(tr_dens$Survey, ref="Hchijojima1")
# we drop 'Location' variable as arbitrary.
#Also drop D_LAND

glmden1<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + 
                 D_KURO + SST + G_SST + CHLA + MONTH + offset(log(perc1km_surv))+
                 (1|Survey), data=tr_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(glmden1, type='pearson')^2)/df.residual(glmden1)) #4.73
cor(fitted(glmden1), tr_dens$Density) #0.382
plot(glmden1) #average


glmden_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                      G_SST + CHLA + MONTH + offset(log(perc1km_surv))+
                      (1|Survey), data=tr_dens) 
# does not converge

# try without RE

glmden_nb<-glm.nb(Density~BATHY + SLOPE + D_COL +  D_KURO + poly(SST,2) +
                    G_SST + CHLA + MONTH + offset(log(perc1km_surv))+
                    Survey, data=tr_dens) 

sum(resid(glmden_nb, type='pearson')^2)/df.residual(glmden_nb) #1.27
cor(fitted(glmden_nb), tr_dens$Density)#0.3

#try glmmadmb
library(glmmADMB)
# D_LAND and COL_INF dropped cos theyre crap and collinear
admb_mod<-glmmadmb(Density~BATHY + SLOPE + D_COL 
                    D_KURO + SST + G_SST + CHLA + MONTH + offset(log(perc1km_surv))+
                     (1|Survey), zeroInflation=T, data=tr_dens, family="poisson")

sum(resid(admb_mod, type='pearson')^2)/df.residual(admb_mod) #2.43
cor(fitted(admb_mod), tr_dens$Density) #0.31
# still a little overdispursed but can live with it


# see if any varibs need a poly
library(reshape2)
d1<-melt(tr_dens, id.vars=c("Density", "St_date","Location", "Survey", "MONTH","DAYNIGHT", "YEAR"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=Density, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")



m_lin<-glm(Density~SST,
                data=tr_dens, family="poisson")
m_pol<-glm(Density~poly(SST,2),
             data=tr_dens, family="poisson")
d1<-data.frame(dens=tr_dens$Density, SST=tr_dens$SST, m_lin=predict(m_lin, new_data=tr_dens, type="response"), m_pol=predict(m_pol, new_data=tr_dens, type="response"))
g1<-ggplot(data=d1, aes(y=dens, x=SST))
g1+geom_point()+geom_line(aes(y=m_lin, x=SST), colour=2)+
  geom_line(aes(y=m_pol, x=SST), colour=3)
anova(m_lin, m_pol)
#poly model better for sst

#update model
admb_final<-glmmadmb(Density~BATHY + SLOPE + D_COL +
                     D_KURO + poly(SST,2) + G_SST + CHLA + MONTH +
                      offset(log(perc1km_surv))+ (1|Survey),
                   zeroInflation=T, data=tr_dens, family="poisson")

sum(resid(admb_mod, type='pearson')^2)/df.residual(admb_mod) #2.62
cor(fitted(admb_mod), tr_dens$Density) #0.23
summary(admb_final)

# Negative binomial

admb_nb<-glmmadmb(Density~BATHY + SLOPE + D_COL +
                       D_KURO + poly(SST,2) + G_SST + CHLA + MONTH +
                       offset(log(perc1km_surv))+ (1|Survey),
                     zeroInflation=F, data=tr_dens, family="nbinom")

sum(resid(admb_nb, type='pearson')^2)/df.residual(admb_nb) #1.23
cor(fitted(admb_nb), tr_dens$Density) #0.313
summary(admb_nb)

# negative binomial the way to go :)


### check spatial autocorrelation
library(ncf)
spac_final <- spline.correlog(tr_dens$Longitude, tr_dens$Latitude,
                                 residuals(admb_nb, type="pearson")[,1],
                             xmax=200, resamp=5,latlon=TRUE)
plot(spac_final)
#crank up resamp
# doesnt seem to be much SPAC
write.csv(data.frame(
CI_5perc=spac_final$boot$boot.summary$predicted$y[3,],
mean=spac_final$boot$boot.summary$predicted$y[6,],
CI_95perc=spac_final$boot$boot.summary$predicted$y[9,]),
"~/research/miller_et_al/remodelling_2017/global_nb_SPAC.csv",
quote=F, row.names=F)
# seems to be fine

## global model - used to make sure model is good enough for IT inference
## Model evaluation of global model


## Calculate Pseudo R squared use the function Jarred Byrners wrote.

r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  
  summary(lmfit)$r.squared  }

r2.corr.mer(gamden_global$mer)
# 0.9955187

r2.corr.mer(gamden_global_RAC$mer)
#0.9937483


# fixed effects only

lmfit <-  lm(model.response(model.frame(gamden_global_RAC$gam)) ~ predict(gamden_global_RAC$gam, newdata=tr_dens))

summary(lmfit)$r.squared #0.08086817 

#confirm
library(MuMIn)
r.squaredGLMM(admb_nb)

cor.test(fitted(gamden_global$mer), tr_dens$Density)
#95 percent confidence interval:
#0.8000917 0.8218037
#sample estimates:
#cor 
#0.9977568

cor.test(fitted(gamden_global_RAC$mer), tr_dens$Density)
#95 percent confidence interval:
#0.7994455 0.8212207
#sample estimates:
#cor 
#0.9968692

# fixed effects only
cor.test(predict(gamden_global_RAC$gam, newdata=tr_dens), tr_dens$Density)
#95 percent confidence interval:
#  0.2549438 0.3132765
#sample estimates:
#  cor 
#0.2843733 

plot(admb_nb, resid(., type='pearson')~fitted(.), type=c('smooth', 'p'), abline=0)


### NOW for individual IT models ###

## time model - Q.Is JM density driven by timing through incubation and between years?
gamden_time<-gamm4(Density~ + YEAR + MONTH + offset(log(effort)),
                   random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_time,ras=rasTempl, dat=tr_dens)

gamden_time<-gamm4(Density~ + YEAR + MONTH + gam_RAC_term + offset(log(effort)),
                   random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## runoff model - Q. Is JM density determined by runoff from land and coastal productivity? 
gamden_runoff<-gamm4(Density~ s(CHLA, k=5) + s(D_LAND, k=5) + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_runoff,ras=rasTempl, dat=tr_dens)

gamden_runoff<-gamm4(Density~ s(CHLA, k=5) + s(D_LAND, k=5) + gam_RAC_term + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## seabed model - Q. Is JM density driven by seabed habitat and associated prey accesibility (diving)? 
gamden_seabed<-gamm4(Density~ s(BATHY, k=5) + s(SLOPE, k=5) + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_seabed,ras=rasTempl, dat=tr_dens)

gamden_seabed<-gamm4(Density~ s(BATHY, k=5) + s(SLOPE, k=5) + gam_RAC_term + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## static model - Q. Is JM density driven by long term static features that control access to food etc during incubation? 
gamden_static<-gamm4(Density~ s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + s(D_LAND, k=5) + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_static,ras=rasTempl, dat=tr_dens)

gamden_static<-gamm4(Density~ s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + s(D_LAND, k=5) + gam_RAC_term + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## bird model - Q. Is JM density driven by foraging range and colony interractions with other JM? 
gamden_bird<-gamm4(Density~ s(D_COL, k=5) + COL_INF + offset(log(effort)),
                   random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_bird,ras=rasTempl, dat=tr_dens)

gamden_bird<-gamm4(Density~  gam_RAC_term + offset(log(effort)),
                   random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## competition model - Q. Is JM density driven by inter-colony competition for long term resources JM? 
gamden_compet<-gamm4(Density~ s(CHLA, k=5) + s(SST, k=5) + COL_INF + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_compet,ras=rasTempl, dat=tr_dens)

gamden_compet<-gamm4(Density~ s(CHLA, k=5) + s(SST, k=5) + COL_INF + gam_RAC_term + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## dynamic sst - Q. Is JM density driven by SST preference that they track throughout incubation and between years? 
gamden_dyn_sst<-gamm4(Density~ s(SST, k=5) + s(G_SST, k=5) +MONTH + YEAR + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_dyn_sst,ras=rasTempl, dat=tr_dens)

gamden_dyn_sst<-gamm4(Density~ s(SST, k=5) + s(G_SST, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## dynamic chla - Q. Is JM density driven by primary productivity or lower order prey that they track throughout incubation and between years? 
gamden_dyn_chla<-gamm4(Density~ s(CHLA, k=5) + MONTH + YEAR + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_dyn_chla,ras=rasTempl, dat=tr_dens)

gamden_dyn_chla<-gamm4(Density~ s(CHLA, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## dynamic oceanographic - Q. Is JM density driven by dynamic oceanographic conditions throughout incubation and between years? 
gamden_dyn_oceo<-gamm4(Density~ s(CHLA, k=5) + s(SST, k=5) + s(G_SST, k=5) + s(D_KURO, k=5) + MONTH + YEAR + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_dyn_oceo,ras=rasTempl, dat=tr_dens)

gamden_dyn_oceo<-gamm4(Density~ s(CHLA, k=5) + s(SST, k=5) + s(G_SST, k=5) + s(D_KURO, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## dynamic front col limited - Q. Is JM density driven by dynamic current/fronts, limited by foraging range throughout incubation and between years? 
gamden_col_front<-gamm4(Density~ s(D_COL, k=5) + s(D_KURO) + MONTH + YEAR + offset(log(effort)),
                        random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_col_front,ras=rasTempl, dat=tr_dens)

gamden_col_front<-gamm4(Density~ s(D_COL, k=5) + s(D_KURO, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                        random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## bathy col limited - Q. Is JM density driven by primary productivity and water depth (prey?), limited by foraging range throughout incubation and between years? 
gamden_col_bat<-gamm4(Density~ s(D_COL, k=5) + s(BATHY, k=5) + MONTH + YEAR + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_col_bat,ras=rasTempl, dat=tr_dens)

gamden_col_bat<-gamm4(Density~ s(D_COL, k=5) + s(BATHY, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

## dynamic land col limited - Q. Is JM density determined by coastal proximity, limited by foraging range throughout incubation and between years? 
gamden_col_land<-gamm4(Density~ s(D_COL, k=5) + s(D_LAND, k=5) + MONTH + YEAR + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_col_land,ras=rasTempl, dat=tr_dens)

gamden_col_land<-gamm4(Density~ s(D_COL, k=5) + s(D_LAND, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## sst col - Q. Is JM density driven by preference for sst and gradient throughout incubation and between years? 
gamden_sst_col<-gamm4(Density~ s(SST, k=5) + s(G_SST, k=5) + s(D_COL, k=5) + MONTH + YEAR + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_sst_col,ras=rasTempl, dat=tr_dens)

gamden_sst_col<-gamm4(Density~ s(SST, k=5) + s(G_SST, k=5) + s(D_COL, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                      random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## chla col limited - Q. Is JM density driven by preference for sst, limited by foraging range throughout incubation and between years? 
gamden_chla_col<-gamm4(Density~ s(CHLA, k=5) + s(D_COL, k=5) + MONTH + YEAR + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_chla_col,ras=rasTempl, dat=tr_dens)

gamden_chla_col<-gamm4(Density~ s(CHLA, k=5) + s(D_COL, k=5) + MONTH + YEAR + gam_RAC_term + offset(log(effort)),
                       random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


temp<-AIC(      gamden_runoff$mer, 
                gamden_seabed$mer,    gamden_bird$mer, 
                gamden_compet$mer,    gamden_dyn_sst$mer,
                gamden_dyn_chla$mer,  gamden_dyn_oceo$mer,
                gamden_col_front$mer, gamden_col_bat$mer,
                gamden_col_land$mer,  gamden_sst_col$mer,
                gamden_chla_col$mer, gamden_time$mer)


mod.list<-list(      gamden_runoff$mer, 
                     gamden_seabed$mer,    gamden_bird$mer, 
                     gamden_compet$mer,    gamden_dyn_sst$mer,
                     gamden_dyn_chla$mer,  gamden_dyn_oceo$mer,
                     gamden_col_front$mer, gamden_col_bat$mer,
                     gamden_col_land$mer,  gamden_sst_col$mer,
                     gamden_chla_col$mer, gamden_time$mer)


it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]
it.out

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)

write.csv(it.out, "results/IT_habitatuse_model_results_obsRE.csv", quote=F, row.names=F)


#gamden_sst_col model summary
summary(gamden_col_front$mer)

#Parametric coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)   -6.0356     0.4031 -14.974  < 2e-16 ***
#  MONTH04        1.6408     0.4829   3.398 0.000679 ***
#  MONTH05        3.1278     0.5441   5.748 9.02e-09 ***
#  YEAR2009      -1.7521     0.9442  -1.856 0.063488 .  
#  YEAR2010      -4.2124     0.4473  -9.417  < 2e-16 ***
#  YEAR2011      -4.9549     0.4347 -11.399  < 2e-16 ***
#  YEAR2012      -4.4437     0.4787  -9.283  < 2e-16 ***
#  gam_RAC_term  24.0321     1.5415  15.590  < 2e-16 ***

gamden_sst_col$mer
#Random effects:
#  Groups          Name        Variance  Std.Dev.
#obs             (Intercept) 7.274e+00  2.69710
#Survey:Location (Intercept) 1.327e-04  0.01152
#Location        (Intercept) 5.791e-05  0.00761
#Xr.0            s(D_KURO)   6.207e+02 24.91366
#Xr              s(D_COL)    1.059e-01  0.32535

# variation between random effects:

exp((1.327e-04^2)) # variation of JM density between surveys within location= 1?
exp((5.791e-05^2)) # variation of JM density between locations = 1?

## response curve plots ##

#make bar and error bar plots from model predictions for YEAR and MONTH

factor_varb_df<-data.frame(Variable=)

#               Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)  -3.99710    0.23197 -17.231  < 2e-16 ***
#  MONTH04       0.32457    0.25892   1.254   0.2100    
#  MONTH05       1.53097    0.30577   5.007 5.53e-07 ***
#  YEAR2009     -0.97171    0.40441  -2.403   0.0163 *  
#  YEAR2010     -3.14351    0.31018 -10.135  < 2e-16 ***
#  YEAR2011     -3.50652    0.30843 -11.369  < 2e-16 ***
#  YEAR2012     -1.76935    0.27901  -6.342 2.27e-10 ***
#  gam_RAC_term  1.05556    0.04141  25.490  < 2e-16 ***

## year

pred_year<-NULL
for(i in unique(tr_dens$YEAR))
{
  pred_frame<-data.frame(D_KURO=median(tr_dens$D_KURO), D_COL=median(tr_dens$D_COL), 
                         YEAR="2009", MONTH="05", effort=median(tr_dens$effort),
                         gam_RAC_term=0)  
  
  pred_frame$YEAR=i
  
  print(i)
  pred_frame$mod1<-predict(gamden_col_front$gam, type='response', newdata=pred_frame, re.form=~0)
  pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
  pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
  
  pred_frame$YEAR_ID=i
  pred_year<-rbind(pred_year, pred_frame)
}

pred_year$mod1_low[5]<-0.15 ## bodge!!

limits <- aes(ymax = mod1_upp, ymin=mod1_low)
pred_yr <- ggplot(pred_year, aes( y=mod1, x=YEAR_ID)) + geom_bar( stat="identity")+geom_errorbar(limits, width=0.25) + + scale_y_sqrt(breaks=c(0,1,3,6,10,20,30,40)) + xlab("Year") + ylab(expression("JM density"~(birds~'/'~km^{2}))) + theme_bw()

## month

pred_month<-NULL
for(i in unique(tr_dens$MONTH))
{
  pred_frame<-data.frame(D_KURO=median(tr_dens$D_KURO), D_COL=median(tr_dens$D_COL), 
                         YEAR="2008", MONTH="04", effort=median(tr_dens$effort),
                         gam_RAC_term=0)  
  
  pred_frame$MONTH=i
  
  print(i)
  pred_frame$mod1<-predict(gamden_col_front$gam, type='response', newdata=pred_frame, re.form=~0)
  pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
  pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
  
  pred_frame$MONTH_ID=i
  pred_month<-rbind(pred_month, pred_frame)
}

limits <- aes(ymax = mod1_upp, ymin=mod1_low)
pred_mn <- ggplot(pred_month, aes( y=mod1, x=MONTH_ID)) + geom_bar( stat="identity", fill="darkgrey", colour="black")+geom_errorbar(limits, width=0.25) + xlab("Month") + ylab(expression("JM density"~(birds~'/'~km^{2}))) + theme_bw()

# pred_dkuro
pred_frame<-data.frame(id=seq(1, 500), D_KURO=median(tr_dens$D_KURO), D_COL=median(tr_dens$D_COL), 
                       YEAR="2008", MONTH="04", effort=median(tr_dens$effort),
                       gam_RAC_term=0)  


pred_frame$D_KURO=seq(0,13, length.out=500)

pred_frame$mod1<-predict(gamden_col_front$gam, type='response', newdata=pred_frame, re.form=~0)
pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit

pred_frame$D_KURO<-(pred_frame$D_KURO)^2 # reverse log transformation from original data munipulation

pred_dkuro<-ggplot(data=pred_frame, aes(x=D_KURO, y=mod1)) + 
  geom_point(data=tr_dens, aes(y=Density, x=D_KURO^2), size=1.2) +
  geom_smooth(stat='identity',aes(ymin=mod1_low, ymax=mod1_upp), fill="darkgrey",colour="black", size=1) +  scale_x_continuous(limits=c(0,120)) + scale_y_continuous(limits=c(0,20)) + xlab("Distance to Kuroshio current (km)")  + ylab(expression("JM density"~(birds~'/'~km^{2}))) +theme_bw()


# d_col
pred_frame<-data.frame(id=seq(1, 500), D_KURO=median(tr_dens$D_KURO), D_COL=median(tr_dens$D_COL), 
                       YEAR="2008", MONTH="04", effort=median(tr_dens$effort),
                       gam_RAC_term=0)  


pred_frame$D_COL=seq(0.4,6.9, length.out=500)

pred_frame$mod1<-predict(gamden_col_front$gam, type='response', newdata=pred_frame, re.form=~0)
pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_col_front$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit

pred_frame$D_COL<-(pred_frame$D_COL)^2 # reverse log transformation from original data munipulation

pred_dcol<-ggplot(data=pred_frame, aes(x=D_COL, y=mod1)) + 
  geom_point(data=tr_dens, aes(y=Density, x=D_COL^2), size=1.2) +
  geom_smooth(stat='identity',aes(ymin=mod1_low, ymax=mod1_upp), fill="darkgrey", colour="black", size=1)+ scale_x_continuous(limits=c(0,40))  + xlab("Distance to colony (km)")  + ylab(expression("JM density"~(birds~'/'~km^{2}))) +theme_bw()


library(gridExtra)

grid.arrange( pred_dcol, pred_dkuro, pred_mn, widths = c(1,1), ncol=3,nrow=1)

png("D:/BIRDLIFE/miller_et_al/results/habitat_pref_response_plots_obsRE_3varib.png", width = 9, height =6 , units ="in", res =600)

grid.arrange( pred_dcol, pred_dkuro, pred_mn, widths = c(1,1), ncol=3,nrow=1)

dev.off()

## extra looking at impact of d_kuro

pred_frame<-data.frame(v1=rep(1, 500), SST=median(tr_dens$SST), G_SST=median(tr_dens$G_SST), D_KURO=median(tr_dens$D_KURO), CHLA=median(tr_dens$CHLA), 
                       YEAR="2009", MONTH="04", effort=median(tr_dens$effort),
                       gam_RAC_term=0)  

pred_frame$D_KURO=seq(0,13, length.out=500)

pred_frame$mod1<-predict(gamden_dyn_oceo$gam, type='response', newdata=pred_frame, re.form=~0)
pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_dyn_oceo$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_dyn_oceo$gam, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit

pred_frame$D_KURO<-(pred_frame$D_KURO)^2 # reverse log transformation from original data munipulation

pred_dcol<-ggplot(data=pred_frame, aes(x=D_KURO, y=mod1)) + 
  geom_point(data=tr_dens, aes(y=Density, x=D_KURO^2), size=1.2) +
  geom_smooth(stat='identity',aes(ymin=mod1_low, ymax=mod1_upp), colour="black", size=1) + scale_x_continuous(limits=c(0,45)) + scale_y_continuous(limits=c(0,20)) + xlab("Distance to colony (km)")  + ylab(expression("JM density"~(birds~'/'~km^{2}))) +theme_bw()

## bit of gam playing

# google: nabble mgcv:gamm: predict to reflect random s() effects

gamden_gam_1<-gam(Density~ s(SST, k=5) + s(G_SST, k=5) + s(D_COL, k=5) + MONTH + YEAR + offset(log(effort)) +s(Survey, bs="re", by=dum_1), data=tr_dens, family=poisson(link=log), method="REML") 
gamden_gam_2<-gam(Density~ s(SST, k=5) + s(G_SST, k=5) + s(D_COL, k=5) + MONTH + YEAR + offset(log(effort)) +s(Location,  Survey, bs="re", by=dum_1), data=tr_dens, family=poisson(link=log), method="REML") 
gamden_gam_3<-gam(Density~ s(SST, k=5) + s(G_SST, k=5) + s(D_COL, k=5) + MONTH + YEAR + offset(log(effort)) +s(Density, Location,  Survey, bs="re", by=dum_1), data=tr_dens, family=poisson(link=log), method="REML") 

#including density makes really good - just observation level random effect?? equivilent to slope? in any case cannot predict RE with it


# change Location/survey to see random effects

pred_frame<-data.frame(v1=rep(1, 500), SST=median(tr_dens$SST), G_SST=median(tr_dens$G_SST), D_COL=median(tr_dens$D_COL), 
                       YEAR="2008", MONTH="04", effort=median(tr_dens$effort),
                       gam_RAC_term=0, Location="Birojima", Survey="Birojima1", dum_1=1 )  

pred_frame$SST=seq(14.8,20, length.out=500)

pred_frame$mod1<-predict(gamden_gam_1, type='response', newdata=pred_frame, re.form=~0)
pred_frame$mod1_upp<-pred_frame$mod1 + predict(gamden_gam_1, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit
pred_frame$mod1_low<-pred_frame$mod1 - predict(gamden_gam_1, se.fit=T, type='response', newdata=pred_frame, re.form=~0)$se.fit

pred_sst<-ggplot(data=pred_frame, aes(x=SST, y=mod1)) + 
  geom_point(data=tr_dens, aes(y=Density, x=SST), size=1.2) +
  geom_smooth(stat='identity',aes(ymin=mod1_low, ymax=mod1_upp), colour="black", size=1) + scale_x_continuous(breaks=c(15,16,17,18,19,20)) + xlab(expression(paste("Sea surface temperature (", degree ~ C, " )")))  + ylab(expression("JM density"~(birds~'/'~km^{2}))) +theme_bw()
pred_sst

### ahhh cooool nooo waaaay

