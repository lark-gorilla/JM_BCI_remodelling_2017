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

jm_data<-read.csv("input_data/dens_abs_1km_extract.csv", h=T)

nrow(jm_data)
length(which(is.na(jm_data$SST)==TRUE))
length(which(is.na(jm_data$CHLA)==TRUE))

jm_data<-na.omit(jm_data) ## remove NAs
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
jm_data$effort<-0
for(i in unique(jm_data$Survey)){
  jm_data[jm_data$Survey==i,]$effort<-nrow(jm_data[jm_data$Survey==i,])
}


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
trfdf$G_SST<-log(trfdf$G_SST+0.001)
trfdf$CHLA<-log(trfdf$CHLA)
trfdf$G_CHLA<-log(trfdf$G_CHLA)
trfdf$BATHY<-log((trfdf$BATHY+trfdf$BATHY^2)+1)
trfdf$SLOPE<-sqrt(trfdf$SLOPE)
trfdf$D_KURO<-sqrt(trfdf$D_KURO)
trfdf$D_COL<-sqrt(trfdf$D_COL)
trfdf$D_LAND<-sqrt(trfdf$D_LAND)

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


##########******** DENSITY MODEL ********##########
############~~~~~~~************~~~~~~~~~~~#########

tr_dens<-filter(trfdf, Survey!="Oki1" & Survey!="Birojima5" & Survey!="Birojima6")

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

## Go for Density<25 cutoff to make modelling better


glmden1<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                 D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                 (1|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)
g_temp<-glmden1
#failed to converge
glmden1 <- refit(glmden1)
glmerConverge(glmden1) #slightly different from g_temp - does work dod on final model for precise coeficients

print(sum(resid(glmden1, type='pearson')^2)/df.residual(glmden1)) #8.93
cor(fitted(glmden1), tr_dens$Density) #0.317
plot(glmden1) #average

tr_dens$obs<-seq(1:nrow(tr_dens)) # add observation level random effect

glmden2<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                 D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                 (1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log), na.action=na.omit)

print(sum(resid(glmden2, type='pearson')^2)/df.residual(glmden2)) #0.0383
cor(fitted(glmden2), tr_dens$Density) #0.999 ?!!!! amazing!
plot(glmden2) #good but zero infl

glmden2.5<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                   D_KURO + SST + G_SST + CHLA + offset(log(effort))+
                   (1|Location/Survey) + (1|YEAR) +(1|MONTH), data=tr_dens, family=poisson(link=log), na.action=na.omit)

print(sum(resid(glmden2.5, type='pearson')^2)/df.residual(glmden2.5)) #8.919
cor(fitted(glmden2.5), tr_dens$Density) #0.31
plot(glmden2.5) #


glmden3<-glmmadmb(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND + D_KURO + SST + G_SST + CHLA + offset(log(effort))+(1|Location/Survey) + (1|YEAR) + (1|MONTH), zeroInflation=T, data=tr_dens, family="poisson")


##testing cor on independent data

train<-tr_dens[sample(1:nrow(tr_dens), (nrow(tr_dens)/100*75)),]
test<-tr_dens[sample(1:nrow(tr_dens), (nrow(tr_dens)/100*25)),]


glm_tr<-update(glmden2, data=train)

p1<-predict(glm_tr, train, type="response", re.form=~0)
cor(p1, train$Density) #0.255
p2<-predict(glm_tr, train, type="response", re.form=~(1|Location/Survey))
cor(p2, train$Density) #0.256
p3<-predict(glm_tr, train, type="response", re.form=~(1|Location/Survey)+(1|obs))
cor(p3, train$Density) #0.999
p4<-predict(glm_tr, train, type="response", re.form=~(1|obs))
cor(p4, train$Density) #0.999

p5<-predict(glmden2, tr_dens, type="response", re.form=~(1|Location/Survey)+(1|obs))
cor(p5, tr_dens$Density) #0.999

p6<-predict(glmden2, tr_dens, type="response", re.form=~0)
cor(p6, tr_dens$Density) #0.19

gamm1<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
               s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5), random=~
               (1|Location/Survey) +(1|YEAR) +(1|MONTH), data=tr_dens, family=poisson(link=log))



glmden2<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                 D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                 (1|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)
#glmerConverge(glmden2)

glmden3<-glmer(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                 D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                 (1|Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)
#glmerConverge(glmden3)

anova(glmden3, glmden2, glmden1)

sum(resid(glmden1, type='pearson')^2)/df.residual(glmden1) #1.02
sum(resid(glmden2, type='pearson')^2)/df.residual(glmden2) #4.86
sum(resid(glmden3, type='pearson')^2)/df.residual(glmden3) #4.84

## mixed modelling alternatives

glm1<-glm(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
            D_KURO + SST + G_SST + CHLA + YEAR + MONTH + Location + Survey + offset(log(effort))
          , data=tr_dens, family=poisson(link=log), na.action=na.omit)

print(sum(resid(glm1, type='pearson')^2)/df.residual(glm1)) #4.89
cor(fitted(glm1), tr_dens$Density) # 0.3128945

library(MASS)
glm2<-glm.nb(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
               D_KURO + SST + G_SST + CHLA + YEAR + MONTH + Location + Survey + offset(log(effort))
             , data=tr_dens, na.action=na.omit)

print(sum(resid(glm2, type='pearson')^2)/df.residual(glm2)) #1.114596
cor(fitted(glm2), tr_dens$Density) # 0.2905831

library(pscl)
f1<-formula(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
              D_KURO + SST + G_SST + CHLA + YEAR + MONTH  +
              offset(log(effort)) |
              BATHY + SLOPE + D_COL + COL_INF + D_LAND +
              D_KURO + SST + G_SST + CHLA + YEAR +MONTH) # Survey, Location, YEAR and MONTH are collinear. choose which. Survey on own or month + year are best


Zip1<-zeroinfl(f1, dist="poisson", link="logit", data=tr_dens)

Zinb1<-zeroinfl(f1, dist="negbin", link="logit", data=tr_dens)

print(sum(resid(Zip1, type='pearson')^2)/df.residual(Zip1))
cor(fitted(Zip1), tr_dens$Density)

print(sum(resid(Zinb1, type='pearson')^2)/df.residual(Zinb1))
cor(fitted(Zinb1), tr_dens$Density)

Hup<-hurdle(f1, dist="poisson", link="logit", data=tr_dens)

Hub<-hurdle(f1, dist="negbin", link="logit", data=tr_dens)

print(sum(resid(Hup, type='pearson')^2)/df.residual(Hup))
cor(fitted(Hup), tr_dens$Density)

print(sum(resid(Hub, type='pearson')^2)/df.residual(Hub))
cor(fitted(Hub), tr_dens$Density)

AIC(Zip1, Zinb1, Hup, Hub)


glmden_nb<-glmer.nb(Density~BATHY + SLOPE + D_COL + COL_INF + D_LAND +
                      D_KURO + SST + G_SST + CHLA + YEAR + MONTH + offset(log(effort))+
                      (1|Location/Survey), data=tr_dens, na.action=na.omit) # have to remove random slope

sum(resid(glmden_nb, type='pearson')^2)/df.residual(glmden_nb) #1.07 
cor(fitted(glmden_nb), tr_dens$Density)

#gam test
library(mgcv)
library(MASS)

gam1<-gam(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
            s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR +MONTH
          +s(Location,  Survey, bs="re"), data=tr_dens, family=poisson(link=log))

sum(resid(gam1, type='pearson')^2)/df.residual(gam1) #4.6 
cor(fitted(gam1), tr_dens$Density) #0.35
qplot(x=fitted(gam1), y=resid(gam1, type="pearson"))

gam1_nb<-gam(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
               s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR +MONTH
             +s(Location,  Survey, bs="re"), data=tr_dens, family=nb())

sum(resid(gam1_nb, type='pearson')^2)/df.residual(gam1_nb) #1.0403 
cor(fitted(gam1_nb), tr_dens$Density) #0.29912
qplot(x=fitted(gam1_nb), y=resid(gam1_nb, type="pearson"))

gam1_zip<-gam(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
                s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR +MONTH
              +s(Location,  Survey, bs="re"), data=tr_dens, family=ziplss)

sum(resid(gam1_zip)^2)/df.residual(gam1_zip) #1.0403 
cor(fitted(gam1_zip), tr_dens$Density) #0.29912
qplot(x=fitted(gam1_zip), y=resid(gam1_zip))

library(nlme)
gamm1<-gamm(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
              s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5)+ YEAR +MONTH
            +s(Location,  Survey, bs="re"), data=tr_dens, family=nb())




#test for non-linearity in variables 

tr_dens$glmer1_resid<-resid(glmer1)
p<-ggplot(data=tr_dens, aes(x=BATHY, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(BATHY), data=tr_dens)) # 0.1 sig

p<-ggplot(data=tr_dens, aes(x=SLOPE, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(SLOPE), data=tr_dens)) # 0.01 sig

p<-ggplot(data=tr_dens, aes(x=D_COL, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(D_COL), data=tr_dens)) #0.01 sig

p<-ggplot(data=tr_dens, aes(x=COL_INF, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(COL_INF), data=tr_dens)) #factor?

p<-ggplot(data=tr_dens, aes(x=D_LAND, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(D_LAND), data=tr_dens)) #nope

p<-ggplot(data=tr_dens, aes(x=D_KURO, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(D_KURO), data=tr_dens)) #0.0001 sig

p<-ggplot(data=tr_dens, aes(x=SST, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(SST), data=tr_dens)) #0.0001 sig

p<-ggplot(data=tr_dens, aes(x=G_SST, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(G_SST), data=tr_dens)) #nope

p<-ggplot(data=tr_dens, aes(x=CHLA, y=glmer1_resid));p+geom_point()+stat_smooth(formulas=y~s(x),method="gam")
summary(gam(glmer1_resid~s(CHLA), data=tr_dens)) #0.0001 sig

# non-linear gamm version of glmer1

gamm1<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
               s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR + MONTH, random=~
               (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log))

# non-linear alternatives
library(splines)

glmden1_spline<-glmer(Density~ns(BATHY, df=2) + ns(SLOPE, df=2) + ns(D_COL, df=2) + COL_INF + D_LAND +
                        ns(D_KURO, df=2) + ns(SST, df=2) + G_SST + ns(CHLA, df=2) + YEAR + MONTH + offset(log(effort))+
                        (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)

glmden1_spline1.5<-glmer(Density~ns(BATHY, df=3) + ns(SLOPE, df=3) + ns(D_COL, df=3) + COL_INF + D_LAND +
                           ns(D_KURO, df=3) + ns(SST, df=3) + G_SST + ns(CHLA, df=3) + YEAR + MONTH + offset(log(effort))+
                           (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)


glmden1_spline2<-glmer(Density~ns(BATHY, df=4) + ns(SLOPE, df=4) + ns(D_COL, df=4) + COL_INF + D_LAND +
                         ns(D_KURO, df=4) + ns(SST, df=4) + G_SST + ns(CHLA, df=4) + YEAR + MONTH + offset(log(effort))+
                         (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)

glmden1_spline_bs<-glmer(Density~bs(BATHY, df=2) + bs(SLOPE, df=2) + bs(D_COL, df=2) + COL_INF + D_LAND +
                           bs(D_KURO, df=2) + bs(SST, df=2) + G_SST + bs(CHLA, df=2) + YEAR + MONTH + offset(log(effort))+
                           (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)

glmden1_spline2_bs<-glmer(Density~bs(BATHY, df=4) + bs(SLOPE, df=4) + bs(D_COL, df=4) + COL_INF + D_LAND +
                            bs(D_KURO, df=4) + bs(SST, df=4) + G_SST + bs(CHLA, df=4) + YEAR + MONTH + offset(log(effort))+
                            (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)


glmden1_poly<-glmer(Density~poly(BATHY, 2) + poly(SLOPE, 2) + poly(D_COL,2) + COL_INF + D_LAND +
                      poly(D_KURO, 2) + poly(SST, 2) + G_SST + poly(CHLA, 2) + YEAR + MONTH + offset(log(effort))+
                      (1+Density|Location/Survey), data=tr_dens, family=poisson(link=log), na.action=na.omit)

gamden1<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
                 s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR + MONTH + offset(log(effort)),
               random=~(1+Density|Location/Survey), data=tr_dens, family=poisson(link=log)) #worked with random slope??

gamden_nb<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
                   s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR + MONTH + offset(log(effort)),
                 random=~(1|Location/Survey), data=tr_dens, family=nb), link = "log")) # try nb also
# note no random slope in negbin. doesnt work.

AIC(glmden3, glmden2, glmden1, glmden_nb, glmden1_spline, glmden1_spline2, glmden1_spline_bs,  glmden1_spline2_bs, glmden1_poly, gamden1$mer)

sum(resid(gamden1$mer, type='pearson')^2)/df.residual(gamden1$mer) #1.03 # almost as gd as gamden1 

plot(gamden1$mer, resid(., type='pearson')~fitted(.), type=c('smooth', 'p'), abline=0)
plot(glmden1, resid(., type='pearson')~fitted(.), type=c('smooth', 'p'), abline=0)

library(car)
qqPlot(resid(gamden1$mer))
qqPlot(resid(glmer1))

library(lattice)
shapiro.test(ranef(gamden1$mer)$'Survey:Location'[,1])
qqmath(ranef(gamden1$mer)$'Survey:Location'[,1], type=c('p', 'r'))

shapiro.test(ranef(gamden1$mer)$'Survey:Location'[,2])
qqmath(ranef(gamden1$mer)$'Survey:Location'[,2], type=c('p', 'r'))

#trial best 2: glmden1 and; gamden1 with RAC term


glm_RAC_term<-makeRAC(mod=glmden1, ras=rasTempl, dat=tr_dens)

tr_dens<-cbind(tr_dens, glm_RAC_term)

glmden1_RAC<-update(glmden1, ~. + glm_RAC_term)

gam_RAC_term<-makeRAC(mod=gamden1$mer, ras=rasTempl, dat=tr_dens)

tr_dens<-cbind(tr_dens, gam_RAC_term)

gamden1_RAC<-update(gamden1, ~. + gam_RAC_term)

#validation.. kinda

AIC(glmden3, glmden2, glmden1, glmden_nb, gamden1$mer, glmden1_RAC, gamden1_RAC$mer)

glmden1_pred<- predict(glmden1, type="response", newdata=tr_dens, re.form=~0)
glmden_nb_pred<- predict(glmden_nb, type="response", newdata=tr_dens, re.form=~0)
gamden1_pred<- predict(gamden1$gam, type="response", newdata=tr_dens, re.form=~0)

glmden1_re1_pred<- predict(glmden1, type="response", newdata=tr_dens, re.form=~(1+Density|Location)) # no nesing
glmden1_re_pred<- predict(glmden1, type="response", newdata=tr_dens, re.form=~(1+Density|Location/Survey))
glmden2_re_pred<- predict(glmden1, type="response", newdata=tr_dens, re.form=~(1|Location/Survey)) # no density
glmden_nb_re_pred<- predict(glmden_nb, type="response", newdata=tr_dens, re.form=~(1|Location/Survey))
gamden1_re_pred<- predict(gamden1$gam, type="response", newdata=tr_dens, re.form=~(1+Density|Location/Survey)) # only predicts fixed effects

qplot(x=glmden1_pred, y=tr_dens$Density)
cor(glmden1_pred, tr_dens$Density) #0.20

qplot(x=glmden1_re_pred, y=tr_dens$Density)
cor(glmden1_re_pred, tr_dens$Density) #0.76! Density v. important!
cor(glmden2_re_pred, tr_dens$Density) #0.24

qplot(x=glmden_nb_pred, y=tr_dens$Density)
cor(glmden_nb_pred, tr_dens$Density) #0.25
qplot(x=glmden_nb_re_pred, y=tr_dens$Density)
cor(glmden_nb_re_pred, tr_dens$Density) #0.27

qplot(x=gamden1_pred, y=tr_dens$Density)
cor(gamden1_pred, tr_dens$Density) #0.27

qplot(x=fitted(gamden1$mer), y=tr_dens$Density)
cor(fitted(gamden1$mer), tr_dens$Density) #0.81! use fitted. fitted returns prediction with RF effects

library(lmodel2)
caltest<-lmodel2(log(tr_dens$Density+0.000001)~log(fitted(glmden1)), nperm=100)
caltest<-lmodel2(log(tr_dens$Density+0.000001)~log(fitted(gamden1$mer)), nperm=100)


# with full model

#SST=rep(seq(14.5, 20.1, 0.01),3)
varb=rep(seq(min(tr_dens$D_LAND), max(tr_dens$D_LAND), length.out=500),3) #G_SST


pred_frame<-data.frame(SST=median(tr_dens$SST), G_SST=median(tr_dens$G_SST), CHLA=median(tr_dens$CHLA),
                       D_KURO=median(tr_dens$D_KURO), BATHY=median(tr_dens$BATHY), SLOPE=median(tr_dens$SLOPE), D_COL=median(tr_dens$D_COL), 
                       COL_INF=median(tr_dens$COL_INF),D_LAND=varb, YEAR="2008",
                       MONTH="04", effort=median(tr_dens$effort), Location=c(rep("Birojima", 500), rep("Izu", 500), rep("Hchijojima", 500)))


pred_frame$glmden1<-predict(glmden1, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_spline<-predict(glmden1_spline, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_spline2<-predict(glmden1_spline2, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_spline1.5<-predict(glmden1_spline1.5, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_poly<-predict(glmden1_poly, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_spline_bs<-predict(glmden1_spline_bs, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$glmden1_spline2_bs<-predict(glmden1_spline2_bs, type='response', newdata=pred_frame, re.form=~(1|Location))
pred_frame$gamden1<-predict(gamden1$gam, type='response', newdata=pred_frame, re.form=~(1|Location))

ggplot(data=pred_frame, aes(x=D_LAND, y=glmden1, colour=Location)) + 
  scale_x_discrete(tr_dens$D_LAND^2) + 
  geom_point(data=tr_dens, aes(y=jitter(tr_dens$Density, 2), x=tr_dens$D_LAND)) +
  geom_smooth(stat='identity') +theme_bw() 



for (i in c("G_SST",  "D_KURO",  "SLOPE", "D_COL"))
  
{
  varb=rep(seq(min(tr_dens[,which(names(tr_dens)==i)]), max(tr_dens[,which(names(tr_dens)==i)]), length.out=500))
  
  
  pred_frame<-data.frame(SST=median(tr_dens$SST), G_SST=median(tr_dens$G_SST), CHLA=median(tr_dens$CHLA),
                         D_KURO=median(tr_dens$D_KURO), BATHY=median(tr_dens$BATHY), SLOPE=median(tr_dens$SLOPE), D_COL=median(tr_dens$D_COL), 
                         COL_INF=median(tr_dens$COL_INF),D_LAND=median(tr_dens$D_LAND), YEAR="2009",
                         MONTH="05", effort=median(tr_dens$effort), Location=rep("Izu", 500),
                         gam_RAC_term=0)  # Location RE doesnt seem to impact but year and month certainly do!
  
  pred_frame[,which(names(pred_frame)==i)]<-varb
  
  print(i)
  pred_frame$glmden1_spline_full<-predict(gamden_col_front$gam, type='response', newdata=pred_frame, re.form=~0)
  
  names(pred_frame)[which(names(pred_frame)==i)]<-"myvarib"
  names(tr_dens)[which(names(tr_dens)==i)]<-"myvarib"
  
  print(ggplot(data=pred_frame, aes(x=myvarib, y=glmden1_spline_full, colour=Location)) + 
          geom_point(data=tr_dens, aes(y=jitter(tr_dens$Density, 2), x=myvarib)) +
          geom_smooth(stat='identity') +theme_bw())
  
  names(tr_dens)[which(names(tr_dens)=="myvarib")]<-i
  
  readline()
}

#global_model

tr_dens$obs<-seq(1:nrow(tr_dens))

gamden1<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + D_LAND +
                 s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR + MONTH + offset(log(effort)),
               random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) #worked with random slope??

#checking k values

gamtemp<-gam(Density~s(BATHY, k=3) + s(SLOPE, k=3) + s(D_COL, k=3) + COL_INF + s(D_LAND, k=3) +
               s(D_KURO, k=3) + s(SST, k=3) + G_SST + s(CHLA, k=3) + YEAR + MONTH + offset(log(effort)),
             data=tr_dens, family=poisson(link=log))

gam.check(gamtemp) # for most varib, k=3 too low


gamtemp2<-gam(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + s(D_LAND, k=5) +
                s(D_KURO, k=5) + s(SST, k=5) + G_SST + s(CHLA, k=5) + YEAR + MONTH + offset(log(effort)),
              data=tr_dens, family=poisson(link=log))

gam.check(gamtemp2) 

gamtemp3<-gam(Density~s(BATHY, k=10) + s(SLOPE, k=10) + s(D_COL, k=10) + COL_INF + s(D_LAND, k=10) +
                s(D_KURO, k=10) + s(SST, k=10) + s(G_SST, k=10) + s(CHLA, k=10) + YEAR + MONTH + offset(log(effort)),
              data=tr_dens, family=poisson(link=log))

gam.check(gamtemp3) # probs go with 10, nah it takes too long and we want 'ecologically' interpretable outputs: go with 5

## Information theoretic approach to model selection ##

tr_dens$gam_RAC_term<-0
tr_dens$obs<-seq(1, nrow(tr_dens))

library(raster)
rasTempl<-raster("extract_rasters/bathy_1km.tif")

makeRAC2<-function(mod=mod, ras=rasTempl, dat=tr_dens){ 
  values(rasTempl)<-NA
  xy_residuals <-cbind(tr_dens$Longitude, tr_dens$Latitude, resid(mod$mer))
  rasTempl[cellFromXY(rasTempl,xy_residuals[,1:2])]<-xy_residuals[,3]
  focal_rac_rast<-focal(rasTempl, w=matrix(1,3,3), fun = mean, na.rm = TRUE)
  focal_rac_vect<-extract(focal_rac_rast,xy_residuals[,1:2])
  return(focal_rac_vect)}


## global model - used to make sure model is good enough for IT inference
gamden_global<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + s(D_LAND, k=5) +
                       s(D_KURO, k=5) + s(SST, k=5) + s(G_SST, k=5) + s(CHLA, k=5) + YEAR + MONTH + offset(log(effort)),
                     random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 

tr_dens$gam_RAC_term<-makeRAC2(mod=gamden_global,ras=rasTempl, dat=tr_dens)

gamden_global_RAC<-gamm4(Density~s(BATHY, k=5) + s(SLOPE, k=5) + s(D_COL, k=5) + COL_INF + s(D_LAND, k=5) +
                           s(D_KURO, k=5) + s(SST, k=5) + s(G_SST, k=5) + s(CHLA, k=5) + YEAR + MONTH + gam_RAC_term + offset(log(effort)),
                         random=~(1|Location/Survey) + (1|obs), data=tr_dens, family=poisson(link=log)) 


## Model evaluation of global model


#Now variance explained (r2) and  Pearson coefficient 
library(MuMIn)

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
r.squaredGLMM(gamden_global_RAC$mer)

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

sum(resid(gamden_global$mer, type='pearson')^2)/
  df.residual(gamden_global$mer) # gamden_global = 1.034
sum(resid(gamden_global_RAC$mer, type='pearson')^2)/
  df.residual(gamden_global_RAC$mer) # gamden_global_RAC = 0.891

plot(gamden_global$mer, resid(., type='pearson')~fitted(.), type=c('smooth', 'p'), abline=0)
plot(gamden_global_RAC$mer, resid(., type='pearson')~fitted(.), type=c('smooth', 'p'), abline=0)

##spatial autocorrelation
library(ncf)

CorRes.gamden <- spline.correlog(tr_dens$Longitude, tr_dens$Latitude, residuals(gamden_global$mer), resamp=10,latlon=TRUE, quiet=TRUE, filter=TRUE)
#crank up resamp
plot(CorRes.gamden)

CorRes.gamden_RAC <- spline.correlog(tr_dens$Longitude, tr_dens$Latitude, residuals(gamden_global_RAC$mer), resamp=10,latlon=TRUE, quiet=TRUE, filter=TRUE)
#crank up resamp
plot(CorRes.gamden_RAC)


Cor_trial <- correlog(tr_dens$Longitude, tr_dens$Latitude, residuals(gamden_global$mer), increment=2,resamp=10,latlon=TRUE)
Cor_trial_RAC <- correlog(tr_dens$Longitude, tr_dens$Latitude, residuals(gamden_global_RAC$mer), increment=2,resamp=10,latlon=TRUE)

summary(Cor_trial$correlation)
#     Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.1190000 -0.0400700 -0.0002804  0.0068970  0.0329100  0.3395000 

summary(Cor_trial_RAC$correlation)

#     Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0752300 -0.0224400  0.0007324  0.0027070  0.0276500  0.1427000 

ct_df<-data.frame(Distance=c(Cor_trial$mean.of.class, Cor_trial_RAC$mean.of.class),
                  Correlation=c(Cor_trial$correlation, Cor_trial_RAC$correlation),
                  Model=c(rep("GAMM", 177), rep("GAMM_RAC", 177)))

p<-ggplot(ct_df[ct_df$Distance<300,], aes(x=Distance, y=Correlation, group=Model, colour=Model))
p+geom_line()+geom_hline(yintercept=0)+xlab("Distance Km")+ylab("Moran's I value")+ggtitle("Spatial autocorrelation of habitat-use model")+theme_bw()

##save plot for appendix

png("D:/BIRDLIFE/miller_et_al/results/FIG_appendix1_habitat_use_mod_SPAC_obsRE.png", 
    width =8, height = 4, units ="in", res =600)

p+geom_line()+geom_hline(yintercept=0)+xlab("Distance Km")+ylab("Moran's I value")+ggtitle("Spatial autocorrelation of habitat-use model")+theme_bw()

dev.off() # close device to save image


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

