rm(list=ls())
setwd("~/research/miller_et_al/")

#libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(MASS)
library(vegan)

boxes <- function(Values)
{
  L <- length(Values)
  par(mfrow=c(2,ceiling(L/2)), mai=c(0.3,0.3,0.2,0.2))
  
  for(i in 1:L) {boxplot(Values[,i], main=names(Values)[i])}
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

trfdf<-jm_data
trfdf$G_SST<-log(trfdf$G_SST+0.001)
trfdf$CHLA<-log(trfdf$CHLA)
trfdf$G_CHLA<-log(trfdf$G_CHLA)
trfdf$BATHY<-log((trfdf$BATHY+trfdf$BATHY^2)+1)
trfdf$SLOPE<-sqrt(trfdf$SLOPE)
trfdf$D_KURO<-sqrt(trfdf$D_KURO)
trfdf$D_COL<-sqrt(trfdf$D_COL)
trfdf$D_LAND<-sqrt(trfdf$D_LAND)



############~~~~~~~~~~  Modelling  ~~~~~~~~~~~#########
############~~~~~~~~~~ *********** ~~~~~~~~~~~#########

#IZU

izu_dens<-trfdf[trfdf$Location=="Izu",]

# check for outliers

boxes(izu_dens[,7:16])

izu_dens<-izu_dens[izu_dens$SLOPE<5,]
izu_dens<-izu_dens[izu_dens$CHLA<2.5,]

boxes(izu_dens[,7:16])

# standardize for lme4 - will need to untransform for prediction

#enviro_std<-decostand(izu_dens[,7:16], method="standardize")
#izu_dens<-data.frame(izu_dens[,1:6], enviro_std, izu_dens[,17:21])


############~~~~~~~~~~ MULTICOLLINEARITY ANALYSES ~~~~~~~~~~~#########
############~~~~~~~~~~ ************************** ~~~~~~~~~~~#########


## will remove G_CHLA from further analyses as it is a product of the original CHLA data
# also D_LAND and COL_INF as >0.5

hist(izu_dens$Density)

izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                 G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #15
cor(fitted(izu_poi), izu_dens$Density) #0.382
plot(izu_poi) #bad

# kill big residuals (v high counts)
izu_dens[which(resid(izu_poi, type="pearson")>15),]

izu_dens<-izu_dens[-which(resid(izu_poi, type="pearson")>15),]

#refit
izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                 G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #4.58
cor(fitted(izu_poi), izu_dens$Density) #0.382
plot(izu_poi) #still overdisp

# Use negative binomial

izu_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                   G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                   (1|Survey), data=izu_dens) 
# does not converge
print(sum(resid(izu_nb, type='pearson')^2)/df.residual(izu_nb)) #0.93
cor(fitted(izu_nb), izu_dens$Density) #31
plot(izu_nb) #ok

#try glmmadmb
library(glmmADMB)
# D_LAND and COL_INF dropped cos theyre crap and collinear
izu_nb2<-glmmadmb(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                    G_SST + CHLA + MONTH + YEAR + offset(log(perc1km_surv))+
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

anova(izu_nb1, izu_nb2, izu_nb3, izu_nb3a, izu_nb4) # ok month best (izu_nb3)


izu_nb5<-glmer.nb(Density~ D_COL + MONTH+ BATHY+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb6<-glmer.nb(Density~ D_COL + MONTH+ SLOPE+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb7<-glmer.nb(Density~ D_COL + MONTH+ D_KURO+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb8<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb9<-glmer.nb(Density~ D_COL + MONTH+ G_SST+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb10<-glmer.nb(Density~ D_COL + MONTH+ CHLA+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)


anova(izu_nb3, izu_nb5, izu_nb6, izu_nb7, izu_nb8, izu_nb9, izu_nb10) # poly(SST) izu_nb8 best

izu_nb11<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ BATHY+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb12<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ SLOPE+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
izu_nb13<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ G_SST+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
izu_nb14<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ D_KURO+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
izu_nb15<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ CHLA+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)

anova(izu_nb8, izu_nb11, izu_nb12, izu_nb13, izu_nb14, izu_nb15) # slope izu_nb12 best

izu_nb16<-glmer.nb(Density~ D_COL + MONTH+ poly(SST,2)+ SLOPE+G_SST+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
# only combo worth testing 

anova(izu_nb12, izu_nb16)

drop1(izu_nb16, test="Chisq") #make sure we cant get rid of one

print(sum(resid(izu_nb16, type='pearson')^2)/df.residual(izu_nb16)) # 0.91
library(piecewiseSEM)
sem.model.fits(izu_nb16)
#     Class            Family Link    n  Marginal Conditional
#1 glmerMod Negative Binomial  log 2981 0.6363992    0.671492
summary(izu_nb16)

#######################################################

#BIRO

bir_dens<-trfdf_std[trfdf_std$Location=="Birojima",]

hist(bir_dens$Density)

bir_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                 G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=bir_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(bir_poi, type='pearson')^2)/df.residual(bir_poi)) #21
cor(fitted(bir_poi), bir_dens$Density) #0.382
plot(bir_poi) #bad

# kill big residuals (v high counts)
bir_dens[which(resid(bir_poi, type="pearson")>15),]

bir_dens<-bir_dens[-which(resid(bir_poi, type="pearson")>15),]

#refit
bir_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                 G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=bir_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(bir_poi, type='pearson')^2)/df.residual(bir_poi)) #3.67
cor(fitted(bir_poi), bir_dens$Density) #0.58
plot(bir_poi) #still overdisp

# Use negative binomial

bir_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                   G_SST + CHLA + MONTH + YEAR+offset(log(perc1km_surv))+
                   (1|Survey), data=bir_dens) 

print(sum(resid(bir_nb, type='pearson')^2)/df.residual(bir_nb)) #1.56
cor(fitted(bir_nb), bir_dens$Density) #0.45
plot(bir_nb) #ok

#try glmmadmb
library(glmmADMB)
# D_LAND and COL_INF dropped cos theyre crap and collinear
bir_nb2<-glmmadmb(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                    G_SST + CHLA + MONTH + YEAR + offset(log(perc1km_surv))+
                    (1|Survey), zeroInflation=F,
                  data=bir_dens, family="nbinom")

sum(resid(bir_nb2, type='pearson')^2)/df.residual(bir_nb2) 
cor(fitted(bir_nb2), bir_dens$Density)
plot(bir_nb2)
# fine

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

anova(bir_nb1, bir_nb2, bir_nb3, bir_nb4) # ok month best (bir_nb4)


bir_nb5<-glmer.nb(Density~ D_COL+YEAR + MONTH+ BATHY+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
bir_nb6<-glmer.nb(Density~ D_COL+YEAR + MONTH+ SLOPE+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
bir_nb7<-glmer.nb(Density~ D_COL+YEAR + MONTH+ D_KURO+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
bir_nb8<-glmer.nb(Density~ D_COL+YEAR + MONTH+ poly(SST,2)+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
bir_nb9<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ G_SST+offset(log(perc1km_surv))+
                    (1|Survey), data=bir_dens)
bir_nb10<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ CHLA+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)

anova(bir_nb3, bir_nb5, bir_nb6, bir_nb7, bir_nb8, bir_nb9, bir_nb10) # poly(SST) bir_nb8 best
anova(bir_nb8,bir_nb5) # taking nb8 over 5


bir_nb11<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ poly(SST,2)+ BATHY+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)
bir_nb12<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ poly(SST,2)+ SLOPE+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)
bir_nb13<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ poly(SST,2)+ G_SST+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)
bir_nb14<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ poly(SST,2)+ D_KURO+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)
bir_nb15<-glmer.nb(Density~ D_COL +YEAR+ MONTH+ poly(SST,2)+ CHLA+offset(log(perc1km_surv))+
                     (1|Survey), data=bir_dens)

anova(bir_nb8, bir_nb11, bir_nb12, bir_nb13, bir_nb14, bir_nb15) # bathy nb11 is best

drop1(bir_nb11, test="Chisq") #make sure we cant get rid of one

print(sum(resid(bir_nb11, type='pearson')^2)/df.residual(bir_nb11)) # 1.6
library(piecewiseSEM)
sem.model.fits(bir_nb11)
#     Class            Family Link    n  Marginal Conditional
#1 glmerMod Negative Binomial  log 593 0.8595552   0.8595552
summary(bir_nb11)


###################################

#Hjijojima

haj_dens<-trfdf_std[trfdf_std$Location=="Hchijojima",]

hist(haj_dens$Density)

haj_poi<-glm(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                 G_SST + CHLA +offset(log(perc1km_surv)),
                data=haj_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(haj_poi, type='pearson')^2)/df.residual(haj_poi)) #15
cor(fitted(haj_poi), haj_dens$Density) #0.382
plot(haj_poi) #bad

# kill big residuals (v high counts) nope, barely any
#haj_dens[which(resid(haj_poi, type="pearson")>15),]

#haj_dens<-haj_dens[-which(resid(haj_poi, type="pearson")>15),]

#refit
haj_poi<-glm(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
               G_SST + CHLA +offset(log(perc1km_surv)),
             data=haj_dens, family=poisson(link=log))

print(sum(resid(haj_poi, type='pearson')^2)/df.residual(haj_poi)) #4.58
cor(fitted(haj_poi), haj_dens$Density) #0.382
plot(haj_poi) #still overdisp

# Use negative binomial

haj_nb<-glm.nb(Density~BATHY + SLOPE + D_COL + D_KURO + poly(SST,2) +
                   G_SST + CHLA +offset(log(perc1km_surv)),
                  data=haj_dens) 

print(sum(resid(haj_nb, type='pearson')^2)/df.residual(haj_nb)) #0.93
cor(fitted(haj_nb), haj_dens$Density) #31
plot(haj_nb) #ok

# NOPE to few data points!!
table(haj_dens$Density)

########## Oki

oki<-trfdf_std[trfdf_std$Location=="Oki",]

table(oki$Density)# v few points, the count of 3 jm is spurious 
# just an artefact of 1km grid aggregation, data is actually PA
oki$PA<-0
oki[oki$Density>0,]$PA<-1
table(oki$PA)

oki_bn<-glm(PA~BATHY + SLOPE + D_COL + SST +
               G_SST + CHLA,
             data=oki, family=binomial)
print(sum(resid(oki_bn, type='pearson')^2)/df.residual(oki_bn))

plot(oki_bn) #o

drop1(oki_bn, test="Chisq") #get idea of who to drop

bir_bn1<-glm(PA~ D_COL, data=oki, family=binomial)
bir_bn2<-glm(PA~ D_COL +BATHY, data=oki, family=binomial)
bir_bn3<-glm(PA~ D_COL +SLOPE, data=oki, family=binomial)
bir_bn4<-glm(PA~ D_COL +SST, data=oki, family=binomial)
bir_bn5<-glm(PA~ D_COL +G_SST, data=oki, family=binomial)
bir_bn6<-glm(PA~ D_COL +CHLA, data=oki, family=binomial)

#D_col is the only one worth keeping.

print(sum(resid(bir_bn1, type='pearson')^2)/df.residual(bir_bn1))
summary(bir_bn1)



############## end for today


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

