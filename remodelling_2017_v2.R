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


########## DATASET MODELLING ###########
########################################

#IZU
izu_dens<-jm_data[jm_data$Location=="Izu",]

# check for outliers

boxes(izu_dens[,7:16])

izu_dens<-izu_dens[izu_dens$SLOPE<5,]
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
cor(fitted(izu_poi), izu_dens$Density) #0.418
plot(izu_poi) #bad

# kill big residuals (v high counts)
izu_dens[which(resid(izu_poi, type="pearson")>15),]

izu_dens<-izu_dens[-which(resid(izu_poi, type="pearson")>15),]

#refit
izu_poi<-glmer(Density~BATHY + SLOPE + D_COL + D_KURO + SST +
                 G_SST + MONTH + YEAR+offset(log(perc1km_surv))+
                 (1|Survey), data=izu_dens, family=poisson(link=log))
# sweet convereges!

print(sum(resid(izu_poi, type='pearson')^2)/df.residual(izu_poi)) #4.95
cor(fitted(izu_poi), izu_dens$Density) #0.29
plot(izu_poi) #still overdisp

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

anova(izu_nb1, izu_nb2, izu_nb2a, izu_nb3, izu_nb3a, izu_nb4) # ok month best (izu_nb3)

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


add1(izu_nb3a,scope= ~D_COL*MONTH +BATHY + SLOPE + D_KURO + SST +G_SST,
     test="Chisq")

izu_nb5<-glmer.nb(Density~ D_COL *MONTH+ BATHY+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb6<-glmer.nb(Density~ D_COL * MONTH+ SLOPE+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb7<-glmer.nb(Density~ D_COL * MONTH+ D_KURO+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb8<-glmer.nb(Density~ D_COL * MONTH+ SST+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb9<-glmer.nb(Density~ D_COL * MONTH+ G_SST+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)


anova(izu_nb3a, izu_nb5, izu_nb6, izu_nb7, izu_nb8, izu_nb9) # poly(SST) izu_nb8 best

izu_nb11<-glmer.nb(Density~ D_COL * MONTH+ poly(SST,2)+ BATHY+offset(log(perc1km_surv))+
                    (1|Survey), data=izu_dens)
izu_nb12<-glmer.nb(Density~ D_COL * MONTH+ poly(SST,2)+ SLOPE+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
izu_nb13<-glmer.nb(Density~ D_COL * MONTH+ poly(SST,2)+ G_SST+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)
izu_nb14<-glmer.nb(Density~ D_COL * MONTH+ poly(SST,2)+ D_KURO+offset(log(perc1km_surv))+
                     (1|Survey), data=izu_dens)

anova(izu_nb8, izu_nb11, izu_nb12, izu_nb13, izu_nb14) # slope izu_nb12 best

# not trying adding D_kuro as its corr with SST

anova(izu_nb12) # final mod

drop1(izu_nb12, test="Chisq") #make sure we cant get rid of one

print(sum(resid(izu_nb12, type='pearson')^2)/
        df.residual(izu_nb12)) # 1.01

sem.model.fits(izu_nb12)
#     Class            Family Link    n  Marginal Conditional
#1 glmerMod Negative Binomial  log 2981 0.6363992    0.671492
summary(izu_nb12)

# test spatial autocorrelation

spac_izu <- spline.correlog(izu_dens$Longitude, izu_dens$Latitude,
                              residuals(izu_nb12, type="pearson"),
                              xmax=250, resamp=5,latlon=TRUE)
plot(spac_izu)

# doesnt seem to be much SPAC
write.csv(data.frame(
  CI_5perc=spac_izu$boot$boot.summary$predicted$y[3,],
  mean=spac_izu$boot$boot.summary$predicted$y[6,],
  CI_95perc=spac_izu$boot$boot.summary$predicted$y[9,]),
  "~/research/miller_et_al/remodelling_2017/izu_SPAC.csv",
  quote=F, row.names=F)
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

