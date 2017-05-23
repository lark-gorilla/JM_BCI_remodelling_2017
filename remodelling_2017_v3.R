rm(list=ls())
setwd("~/research/miller_et_al/")

#libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(MASS)
library(vegan)
library(piecewiseSEM)
library(ncf)

boxes <- function(Values)
{
  L <- length(Values)
  par(mfrow=c(2,ceiling(L/2)), mai=c(0.3,0.3,0.2,0.2))
  
  for(i in 1:L) {boxplot(Values[,i], main=names(Values)[i])}
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


#lets go!

jm_data<-read.csv("input_data/dens_abs_1km_extract_correct.csv", h=T)


nrow(jm_data)
length(which(is.na(jm_data$SST)==TRUE))
length(which(is.na(jm_data$CHLA)==TRUE))
jm_data<-na.omit(jm_data) ## remove NAs, wow ok good I did this! 
# na.omit is integral to information theoretic approach1
#jm_data[jm_data$BATHY>0,]$BATHY<-0

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


########## Invidual dataset MODELLING ###########
########################################

#IZU
izu_dens<-jm_data[jm_data$Location=="Izu",]

# check for outliers

boxes(izu_dens[,7:16])

izu_dens<-izu_dens[izu_dens$SLOPE<5,]
izu_dens<-izu_dens[izu_dens$G_SST<2,]
#izu_dens<-izu_dens[izu_dens$CHLA<2.5,]
izu_raw<-izu_dens
boxes(izu_dens[,7:16])

# standardize for lme4 - will need to untransform for prediction
enviro_std<-decostand(izu_dens[,7:16], method="standardize")
izu_dens<-data.frame(izu_dens[,1:6], enviro_std, izu_dens[,17:21])

#reverse standardization
#dcol<-decostand(izu_raw$D_COL, method="standardize")
#head(dcol*attr(dcol, "scaled:scale")+attr(dcol, "scaled:center"))

############~~~~~~~~~~ MULTICOLLINEARITY ANALYSES ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########

pairs(data.frame(izu_dens$BATHY, izu_dens$SLOPE, izu_dens$D_COL, izu_dens$COL_INF,
                 izu_dens$D_LAND, izu_dens$D_KURO, izu_dens$SST, izu_dens$G_SST, izu_dens$CHLA,
                 izu_dens$G_CHLA, izu_dens$YEAR, izu_dens$MONTH, izu_dens$Density), upper.panel = panel.smooth,lower.panel=panel.cor)

#remove G_CHLA and CHLA :(, D_LAND and COL_INF
# remember that D_KURO aned SST are 0.6 corr..

#Check my if sst poly is needed?

library(reshape2)
d1<-melt(izu_dens, id.vars=c("Density", "St_date","Location", "Survey", "MONTH","DAYNIGHT", "YEAR"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=Density, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")

m_lin<-glm(Density~SST,
           data=izu_dens, family="poisson")
m_pol<-glm(Density~poly(SST,2),
           data=izu_dens, family="poisson")
d1<-data.frame(dens=izu_dens$Density, SST=izu_dens$SST, 
               m_lin=predict(m_lin, new_data=izu_dens, type="response"),
               m_pol=predict(m_pol, new_data=izu_dens, type="response"))
g1<-ggplot(data=d1, aes(y=dens, x=SST))
g1+geom_point()+geom_line(aes(y=m_lin, x=SST), colour=2)+  geom_line(aes(y=m_pol, x=SST), colour=3)
anova(m_lin, m_pol)

# actually poly might not be that good..
hist(izu_dens[izu_dens$Density==0,]$SST)
hist(izu_dens[izu_dens$Density>0,]$SST)

############~~~~~~~~~~  Modelling  ~~~~~~~~~~~#########
############~~~~~~~~~~ *********** ~~~~~~~~~~~#########

hist(izu_dens$Density)

izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                 G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #9.98
#overdispursed
cor(fitted(izu_poi), izu_dens$Density) #0.418
plot(izu_poi) #bad


# Use negative binomial

izu_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                   G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                   (1|Survey), data=izu_dens) 
# does not converge
print(sum(resid(izu_nb, type='pearson')^2)/df.residual(izu_nb)) #0.93
cor(fitted(izu_nb), izu_dens$Density) #0.4
plot(izu_nb) #ok

#try glmmadmb
library(glmmADMB)
izu_nb2<-glmmadmb(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                    G_SST + MONTH + YEAR + offset(log(perc1km_surv))+
                    (1|Survey), zeroInflation=F,
                  data=izu_dens, family="nbinom")

sum(resid(izu_nb2, type='pearson')^2)/df.residual(izu_nb2) #0.93
cor(fitted(izu_nb2), izu_dens$Density) #0.31
plot(izu_nb2)
# fine

# Model selection forwards stepwise using D_col first, then adding year then month,
# then testing each oceanographic varib and another combo

drop1(izu_nb, test="Chisq") #get idea of who to drop

izu_nb1<-glmer.nb(Density~ D_COL +offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens) 
izu_nb2<-glmer.nb(Density~ D_COL + YEAR+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens) 
izu_nb3<-glmer.nb(Density~ D_COL + MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens) 
izu_nb4<-glmer.nb(Density~ D_COL +YEAR+ MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens) 

izu_nb3a<-glmer.nb(Density~ D_COL * MONTH+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens) 
izu_nb3b<-glmer.nb(Density~ D_COL : MONTH+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens) 
izu_nb2a<-glmer.nb(Density~ D_COL * YEAR+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens) 

anova(izu_nb1, izu_nb2, izu_nb2a, izu_nb3, izu_nb3a, izu_nb3b, izu_nb4) # ok month best (izu_nb3)

# evidence for month and month:D_col interaction, not year
Anova(izu_nb2)
Anova(izu_nb2a)

# check if interaction is doing what I think it is

# reverse standardization
dcol<-decostand(izu_raw$D_COL, method="standardize")
# (x - mean(x)) / sd(x)

D_COL_std<-(((1:50)-attr(dcol, "scaled:center"))/attr(dcol, "scaled:scale"))

new_dat<-rbind(data.frame(D_COL=D_COL_std, MONTH='04', perc1km_surv=50),
               data.frame(D_COL=D_COL_std, MONTH='05', perc1km_surv=50))

p1<-predict(izu_nb3, newdata=new_dat, type="response", re.form=~0)
p2<-predict(izu_nb3a, newdata=new_dat, type="response", re.form=~0)
p3<-predict(izu_nb3b, newdata=new_dat, type="response", re.form=~0)


out<-rbind(data.frame(new_dat, mod="non_intr", pred=p1),
           data.frame(new_dat, mod="intr", pred=p2),
           data.frame(new_dat, mod="intr_only", pred=p3))

qplot(data=out, x=D_COL, y=pred, colour=MONTH, linetype=mod, geom="line")

# also
ggplot(data=izu_raw[izu_raw$Density<20,], aes(y=Density, x=D_COL))+
  +   geom_jitter( aes(y=Density, colour=MONTH), width=0, height=1)


# see which other variables could be added to the model

add1(izu_nb3a,scope= ~D_COL*MONTH +BATHY + SLOPE + D_KURO + SST +G_SST,
     test="Chisq")
# evidence for G_SST and Slope
# lets investigate

qplot(data=izu_raw[izu_raw$G_SST<1,], y=Density, x=G_SST,
      geom="point")+geom_smooth(method="glm")+
  geom_jitter(width=0, height=1) # outliers could be contributing to the pattern..

ggplot(data=izu_raw, aes(y=Density, x=G_SST))+
  geom_density(data=izu_raw[izu_raw$Density<1,],
               aes(x=G_SST, (..scaled..)))+
  geom_density(data=izu_raw[izu_raw$Density>0,],
               aes(x=G_SST, (..scaled..)/2, colour=2))
# not overly pschyed with G_SST

ggplot(data=izu_raw, aes(y=Density, x=SLOPE))+
  geom_density(data=izu_raw[izu_raw$Density<1,],
               aes(x=SLOPE, (..scaled..)))+
  geom_density(data=izu_raw[izu_raw$Density>0,],
               aes(x=SLOPE, (..scaled..)/2, colour=2))+
  geom_point(data=izu_raw[izu_raw$Density>0,], aes(y=Density/90))
# not overly pschyed with G_SST
# not overly pschyed with Slope either

# model with SLOPE

izu_nb5<-glmer.nb(Density~ D_COL * MONTH+ SLOPE+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens) 

anova(izu_nb3a, izu_nb5) # slope adds a bit

# any others

add1(izu_nb5,scope= ~D_COL*MONTH +SLOPE+BATHY + D_KURO + SST +G_SST,
     test="Chisq")

izu_nb6<-glmer.nb(Density~ D_COL * MONTH+ SLOPE+G_SST+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens) 

anova(izu_nb3a, izu_nb5, izu_nb6) 

summary(izu_nb6)

Anova(izu_nb6) # final mod

print(sum(resid(izu_nb6, type='pearson')^2)/
        df.residual(izu_nb6))

summary(izu_nb3a)

Anova(izu_nb3a) # final mod

print(sum(resid(izu_nb3a, type='pearson')^2)/
        df.residual(izu_nb3a)) 


sem.model.fits(izu_nb3a)
sem.model.fits(izu_nb6)


# test spatial autocorrelation

spac_izu <- spline.correlog(izu_dens$Longitude, izu_dens$Latitude,
                            residuals(izu_nb3a, type="pearson"),
                            xmax=250, resamp=5,latlon=TRUE)
plot(spac_izu)
# doesnt seem to be much SPAC
write.csv(data.frame(
  CI_5perc=spac_izu$boot$boot.summary$predicted$y[3,],
  mean=spac_izu$boot$boot.summary$predicted$y[6,],
  CI_95perc=spac_izu$boot$boot.summary$predicted$y[9,]),
  "~/research/miller_et_al/remodelling_2017/_col_month_izu_SPAC.csv",
  quote=F, row.names=F)
#######################################################

#BIRO

#IZU
bir_dens<-jm_data[jm_data$Location=="Birojima",]

# check for outliers

boxes(bir_dens[,7:16])

bir_dens<-bir_dens[bir_dens$SLOPE<4,]
bir_dens<-bir_dens[bir_dens$G_SST<2,]
#bir_dens<-bir_dens[bir_dens$CHLA<2.5,]
bir_raw<-bir_dens
boxes(bir_dens[,7:16])

# standardize for lme4 - will need to untransform for prediction
enviro_std<-decostand(bir_dens[,7:16], method="standardize")
bir_dens<-data.frame(bir_dens[,1:6], enviro_std, bir_dens[,17:21])

#reverse standardization
#dcol<-decostand(izu_raw$D_COL, method="standardize")
#head(dcol*attr(dcol, "scaled:scale")+attr(dcol, "scaled:center"))

############~~~~~~~~~~ MULTICOLLINEARITY ANALYSES ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########

pairs(data.frame(bir_dens$BATHY, bir_dens$SLOPE, bir_dens$D_COL, bir_dens$COL_INF,
                 bir_dens$D_LAND, bir_dens$D_KURO, bir_dens$SST, bir_dens$G_SST, bir_dens$CHLA,
                 bir_dens$G_CHLA, bir_dens$YEAR, bir_dens$MONTH, bir_dens$Density), upper.panel = panel.smooth,lower.panel=panel.cor)

#remove G_CHLA and CHLA :(, D_LAND and COL_INF
# remember that YEAR aned SST are 0.7 corr..

############~~~~~~~~~~  Modelling  ~~~~~~~~~~~#########
############~~~~~~~~~~ *********** ~~~~~~~~~~~#########

hist(bir_dens$Density)

bir_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                 G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=bir_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(bir_poi, type='pearson')^2)/df.residual(bir_poi)) #16.14
#overdispursed
cor(fitted(bir_poi), bir_dens$Density) #0.49
plot(bir_poi) #bad big outlier get rid

library(sp)
resglm<-residuals(bir_poi, type="pearson")
sp2<-SpatialPointsDataFrame(bir_dens[,4:5], data=data.frame(resglm), 
                            proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')

# Use negative binomial

bir_nb<-glmer.nb(Density~D_COL + MONTH + YEAR+offset(log(perc1km_surv))+
                   (1|Survey), data=bir_dens) 

print(sum(resid(bir_nb, type='pearson')^2)/df.residual(bir_nb)) #2.52
cor(fitted(bir_nb), bir_dens$Density) #0.4
plot(bir_nb) #ok, but still big outliers

summary(bir_nb)
# no need for random effect with month and year in model

# Use negative binomial

bir_nb2<-glm.nb(Density~D_COL + MONTH + YEAR+offset(log(perc1km_surv)),
                   data=bir_dens) 

print(sum(resid(bir_nb2, type='pearson')^2)/df.residual(bir_nb2)) #2.51 no diff from glmm
cor(fitted(bir_nb2), bir_dens$Density) #0.4
plot(bir_nb2) # not great still

resglm<-residuals(bir_nb2, type="pearson")
sp2<-SpatialPointsDataFrame(bir_dens[,4:5], data=data.frame(resglm), 
                            proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')


bir_dens[which(resid(bir_nb2, type="pearson")>10),]

bir_dens<-bir_dens[-which(resid(bir_poi, type="pearson")>15),]



 #hmm lets check other models and see the 


# Model selection forwards stepwise using D_col first, then adding year then month,
# then testing each oceanographic varib and another combo

drop1(bir_nb, test="Chisq") #get idea of who to drop

bir_nb1<-glmer.nb(Density~ D_COL +offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens) 
bir_nb2<-glmer.nb(Density~ D_COL + YEAR+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens) 
bir_nb3<-glmer.nb(Density~ D_COL + MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens) 
bir_nb4<-glmer.nb(Density~ D_COL +YEAR+ MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens) 

bir_nb3a<-glmer.nb(Density~ D_COL * MONTH+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens) 
bir_nb3b<-glmer.nb(Density~ D_COL : MONTH+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens) 
bir_nb2a<-glmer.nb(Density~ D_COL * YEAR+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens) 

summary(bir_nb1);summary(bir_nb2);summary(bir_nb3);summary(bir_nb4);

# most of the RE variance explained by year, month supports a little, try this
bir_nb4a<-glmer.nb(Density~ D_COL * YEAR+MONTH+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens) 


anova(bir_nb1, bir_nb2, bir_nb2a, bir_nb3, bir_nb3a, bir_nb3b, bir_nb4, bir_nb4a) # ok month best (bir_nb3)

# most of the RE variance explained by year, month supports a little, try this

#lets switch to non RE models for the best ones

bir_glm0<-glm.nb(Density~ D_COL +
                   offset(log(perc1km_surv)),data=bir_dens) 

bir_glm1<-glm.nb(Density~ D_COL +YEAR+ MONTH+
                     offset(log(perc1km_surv)),data=bir_dens) 

bir_glm2<-glm.nb(Density~ D_COL*YEAR+ MONTH+
                     offset(log(perc1km_surv)),data=bir_dens) 

bir_glm3<-glm.nb(Density~ D_COL*YEAR+
                     offset(log(perc1km_surv)),data=bir_dens) 

anova(bir_glm0, bir_glm1, bir_glm2, bir_glm3)

anova(bir_glm1, bir_glm2) # simplist appears better (no intereactions)

print(sum(resid(bir_glm1, type='pearson')^2)/df.residual(bir_glm1)) 
# still overdispursed 2.5

plot( fitted(bir_glm1),resid(bir_glm1, type='pearson'))
bir_dens[which(resid(bir_glm1, type="pearson")>8),] # seems a good cutoff

bir_dens<-bir_dens[-which(resid(bir_glm1, type="pearson")>8),]

# refit and test again
bir_glm0<-glm.nb(Density~ D_COL +
                   offset(log(perc1km_surv)),data=bir_dens) 
bir_glm1<-glm.nb(Density~ D_COL +YEAR+ MONTH+
                   offset(log(perc1km_surv)),data=bir_dens) 
bir_glm2<-glm.nb(Density~ D_COL*YEAR+ MONTH+
                   offset(log(perc1km_surv)),data=bir_dens) 
bir_glm3<-glm.nb(Density~ D_COL*YEAR+
                   offset(log(perc1km_surv)),data=bir_dens) 
bir_glm4<-glm.nb(Density~ D_COL*YEAR+D_COL*MONTH+
                   offset(log(perc1km_surv)),data=bir_dens) 
bir_glm5<-glm.nb(Density~ YEAR+D_COL*MONTH+
                   offset(log(perc1km_surv)),data=bir_dens) 

anova(bir_glm0, bir_glm1, bir_glm2, bir_glm3, bir_glm4)
AIC(bir_glm0, bir_glm1, bir_glm2, bir_glm3, bir_glm4)

print(sum(resid(bir_glm1, type='pearson')^2)/df.residual(bir_glm1)) 
print(sum(resid(bir_glm2, type='pearson')^2)/df.residual(bir_glm2)) 
print(sum(resid(bir_glm3, type='pearson')^2)/df.residual(bir_glm3)) 
print(sum(resid(bir_glm4, type='pearson')^2)/df.residual(bir_glm4)) 
# all look good

# just check no RE needed
bir_nb4<-glmer.nb(Density~ D_COL +YEAR+ MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
summary(bir_nb4)

anova(bir_glm1, bir_glm5) # no evidence that D_COL * MONTH intr is needed
anova(bir_glm1, bir_glm2) # no evidence that D_COL * YEAR intr is needed

anova(bir_glm1, bir_glm3) # Month is needed.

# check if interaction is doing what I think it is

# reverse standardization
dcol<-decostand(bir_raw$D_COL, method="standardize")
# (x - mean(x)) / sd(x)

D_COL_std<-(((1:50)-attr(dcol, "scaled:center"))/attr(dcol, "scaled:scale"))

new_dat<-rbind(data.frame(D_COL=D_COL_std, YEAR='2008', MONTH= "04", perc1km_surv=50),
               data.frame(D_COL=D_COL_std, YEAR='2009', MONTH= "04",perc1km_surv=50),
               data.frame(D_COL=D_COL_std, YEAR='2012', MONTH= "04",perc1km_surv=50))

p1<-predict(bir_glm1, newdata=new_dat, type="response")
p2<-predict(bir_glm2, newdata=new_dat, type="response")
p3<-predict(bir_glm3, newdata=new_dat, type="response")


out<-rbind(data.frame(new_dat, mod="non_intr_and_mnt", pred=p1),
           data.frame(new_dat, mod="intr_year_and_mnt", pred=p2),
           data.frame(new_dat, mod="intr_year", pred=p3))

qplot(data=out, x=D_COL, y=pred, colour=YEAR, linetype=mod, geom="line")

# also
ggplot(data=izu_raw[izu_raw$Density<20,], aes(y=Density, x=D_COL))+
  +   geom_jitter( aes(y=Density, colour=MONTH), width=0, height=1)


# see which other variables could be added to the model

add1(bir_glm1,scope= ~D_COL+YEAR+MONTH +BATHY + SLOPE + D_KURO +SST +G_SST,
     test="Chisq")
# evidence for nothing really, SST is corr with Year so not included

# final mod

summary(bir_glm1)

Anova(bir_glm1) # final mod

print(sum(resid(bir_glm1, type='pearson')^2)/
        df.residual(bir_glm1)) 

sem.model.fits(bir_glm1)
cor(fitted(bir_glm1), bir_dens$Density)



# test spatial autocorrelation

spac_bir <- spline.correlog(bir_dens$Longitude, bir_dens$Latitude,
                            residuals(bir_glm1, type="pearson"),
                            xmax=100, resamp=10,latlon=TRUE)
plot(spac_bir)
# doesnt seem to be much SPAC at all 
write.csv(data.frame(
  CI_5perc=spac_bir$boot$boot.summary$predicted$y[3,],
  mean=spac_bir$boot$boot.summary$predicted$y[6,],
  CI_95perc=spac_bir$boot$boot.summary$predicted$y[9,]),
  "~/research/miller_et_al/remodelling_2017/_col__year_month_bir_SPAC.csv",
  quote=F, row.names=F)

########## Invidual dataset MODELLING ###########
########################################

#IZU foraging model
izu_dens<-jm_data[jm_data$Location=="Izu",]

qplot(data=izu_dens, x=D_COL, y=Density, geom="point")
#> 20 km cutoff for foraging  model?

izu_dens$for_cut<-"OUT"
izu_dens[izu_dens$D_COL>20,]$for_cut<-"IN"
qplot(data=izu_dens, x=Longitude, y=Latitude, colour=for_cut, geom="point")
# cool

izu_dens<-izu_dens[izu_dens$D_COL>20,]


# check for outliers

boxes(izu_dens[,7:16])

izu_dens<-izu_dens[izu_dens$SLOPE<10,]
izu_dens<-izu_dens[izu_dens$G_SST<1.5,]
izu_dens<-izu_dens[izu_dens$CHLA<15,]
izu_raw<-izu_dens
boxes(izu_dens[,7:16])


# standardize for lme4 - will need to untransform for prediction
enviro_std<-decostand(izu_dens[,7:16], method="standardize")
# unimportant COl_INF warning
izu_dens<-data.frame(izu_dens[,1:6], enviro_std, izu_dens[,17:21])

#reverse standardization
#dcol<-decostand(izu_raw$D_COL, method="standardize")
#head(dcol*attr(dcol, "scaled:scale")+attr(dcol, "scaled:center"))

############~~~~~~~~~~ MULTICOLLINEARITY ANALYSES ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########

pairs(data.frame(izu_dens$BATHY, izu_dens$SLOPE, izu_dens$D_COL, 
                 izu_dens$D_LAND, izu_dens$D_KURO, izu_dens$SST, izu_dens$G_SST, izu_dens$CHLA,
                 izu_dens$G_CHLA, izu_dens$YEAR, izu_dens$MONTH, izu_dens$Density), upper.panel = panel.smooth,lower.panel=panel.cor)

#remove G_CHLA :(, D_LAND and COL_INF
# remember that D_COL, BATHY and SST are 0.6 corr..

#Check my if sst poly is needed?

library(reshape2)
d1<-melt(izu_dens, id.vars=c("Density", "St_date","Location", "Survey", "MONTH","DAYNIGHT", "YEAR"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=Density, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")


############~~~~~~~~~~  Modelling  ~~~~~~~~~~~#########
############~~~~~~~~~~ *********** ~~~~~~~~~~~#########

hist(izu_dens$Density)

izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                 G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #2.4 not bad
#overdispursed
cor(fitted(izu_poi), izu_dens$Density) #0.19
plot(izu_poi) #not bad


# Use negative binomial

izu_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + D_KURO +
                   G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                   (1|Survey), data=izu_dens) 
# failed to converge
print(sum(resid(izu_nb, type='pearson')^2)/df.residual(izu_nb)) #0.76
cor(fitted(izu_nb), izu_dens$Density) #0.4
plot(izu_nb) #ok

#ok use poisson instead but strip out the big residuals

izu_dens[which(resid(izu_poi, type="pearson")>7),]

izu_dens<-izu_dens[-which(resid(izu_poi, type="pearson")>7),]

#refit
izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                 G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #1.6 fine
cor(fitted(izu_poi), izu_dens$Density) #0.24
plot(izu_poi) #not bad


# Model selection forwards and backwards stepwise using D_col first, then adding year then month,
# then testing each oceanographic varib and another combo


izu_nb1<-glmer(Density~ D_COL +offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens, family=poisson) 
izu_nb2<-glmer(Density~ D_COL + YEAR+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens, family=poisson) 
izu_nb3<-glmer(Density~ D_COL + MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens, family=poisson) 
izu_nb4<-glmer(Density~ D_COL +YEAR+ MONTH+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens, family=poisson) 

anova(izu_nb1, izu_nb2, izu_nb3, izu_nb4) # not much going on there

drop1(izu_poi, test="Chisq") # not very useful

izu_poi1<-glmer(Density~ 1 +offset(log(perc1km_surv))+
                           (1|Survey), data=izu_dens, family=poisson) 

add1(izu_poi1, Density~1+BATHY + SLOPE + D_COL + D_KURO + SST +
       G_SST + MONTH + YEAR+offset(log(perc1km_surv)) ) #hmm

izu_p0<-glmer(Density~ 1+offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson) 
izu_p1<-glmer(Density~ SST +offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson) 
izu_p2<-glmer(Density~ G_SST +offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson) 
izu_p3<-glmer(Density~ CHLA +offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson) 
izu_p4<-glmer(Density~ BATHY +offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson) 
izu_p5<-glmer(Density~ SLOPE +offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson) 
izu_p6<-glmer(Density~ D_KURO +offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson) 
izu_p7<-glmer(Density~ D_LAND +offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson) 

anova(izu_p0, izu_p1, izu_p2, izu_p3, izu_p4, izu_p5, izu_p6, izu_p7)

# Distance to land looking pretty important, but its coor with D_COL.
# check if D_COl is important using CHLA
izu_p3a<-glmer(Density~ CHLA +D_COL+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson)
anova(izu_p3, izu_p3a)

add1(izu_p7, scope= ~D_LAND+YEAR+MONTH + SLOPE + D_KURO +CHLA +G_SST)

izu_p8<-glmer(Density~ D_LAND+CHLA+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson)
anova(izu_p7, izu_p8) # mild evidence for chla inclusion

izu_p8<-glmer(Density~ D_LAND+CHLA+offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson)

table(izu_dens$MONTH)
table(izu_dens$YEAR)

izu_p8a<-glmer(Density~ D_LAND+CHLA+MONTH+offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson)
izu_p8b<-glmer(Density~ D_LAND+CHLA+YEAR+offset(log(perc1km_surv))+
                (1|Survey), data=izu_dens, family=poisson)

anova(izu_p8, izu_p8a)
anova(izu_p8, izu_p8b)






