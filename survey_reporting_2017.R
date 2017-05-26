rm(list=ls())
setwd("~/research/miller_et_al/")

jm_data<-read.csv("input_data/dens_abs_RAW_extract.csv", h=T)

jm_data$count_data<-as.numeric(as.character(jm_data$count_data))

table(jm_data$count_data) #1-90

sd(jm_data$count_data, na.rm=T) #4.23
mean(jm_data$count_data, na.rm=T)#2.89
sum(jm_data$count_data, na.rm=T) #3161

nrow(jm_data[jm_data$pres_abs==1,])

mean(jm_data[jm_data$Location=="Izu",]$count_data, na.rm=T) # 3.1 +- 4.7
mean(jm_data[jm_data$Location=="Birojima",]$count_data, na.rm=T) # 2.2 +-1.6
mean(jm_data[jm_data$Location=="Hchijojima",]$count_data, na.rm=T) # 1.6 =-0.84
mean(jm_data[jm_data$Location=="Oki",]$count_data, na.rm=T) # count data only


jm_data<-read.csv("input_data/dens_abs_1km_extract.csv", h=T)

nrow(jm_data[jm_data$Location=="Izu"&jm_data$Density>0,])/
  nrow(jm_data[jm_data$Location=="Izu",])*100 # 16.8
nrow(jm_data[jm_data$Location=="Birojima"&jm_data$Density>0,])/
  nrow(jm_data[jm_data$Location=="Birojima",])*100 # 21.1
nrow(jm_data[jm_data$Location=="Hchijojima"&jm_data$Density>0,])/
  nrow(jm_data[jm_data$Location=="Hchijojima",])*100 # 3.3
nrow(jm_data[jm_data$Location=="Oki"&jm_data$Density>0,])/
  nrow(jm_data[jm_data$Location=="Oki",])*100 # 6.1

$count_data, na.rm=T) # 2.2
mean(jm_data[jm_data$Location=="Hchijojima",]$count_data, na.rm=T) # 1.6
mean(jm_data[jm_data$Location=="Oki",]$count_data, na.rm=T) # count data only

# remodelling reporting
jm_data<-read.csv("input_data/dens_abs_RAW_extract.csv", h=T)

jm_data$count_data<-as.numeric(as.character(jm_data$count_data))

# number of encounters with JM (individuals or groups)
sum(jm_data[jm_data$Location=="Izu",]$pres_abs, na.rm=T) # 867
sum(jm_data[jm_data$Location=="Birojima",]$pres_abs, na.rm=T) # 330
sum(jm_data[jm_data$Location=="Hchijojima",]$pres_abs, na.rm=T) # 10
sum(jm_data[jm_data$Location=="Oki",]$pres_abs, na.rm=T) # 14

out<-NULL
for(i in c("Izu", "Birojima", "Hchijojima","Oki"))
{
int1<-data.frame(Col=i, 
                SST_ave=mean(jm_data[jm_data$Location==i,]$SST, na.rm=T),
                SST_std=sd(jm_data[jm_data$Location==i,]$SST, na.rm=T), 
                CHL_ave=mean(jm_data[jm_data$Location==i,]$CHLA, na.rm=T),
                CHL_std=sd(jm_data[jm_data$Location==i,]$CHLA, na.rm=T), 
                BTY_ave=mean(jm_data[jm_data$Location==i,]$BATHY, na.rm=T),
                BTY_std=sd(jm_data[jm_data$Location==i,]$BATHY, na.rm=T), 
                DLD_ave=mean(jm_data[jm_data$Location==i,]$D_LAND, na.rm=T),
                DLD_std=sd(jm_data[jm_data$Location==i,]$D_LAND, na.rm=T))
out<-rbind(out, int1)
}

out$n_encounter=0                
out[1,]$n_encounter<-sum(jm_data[jm_data$Location=="Izu",]$pres_abs, na.rm=T) # 867
out[2,]$n_encounter<-sum(jm_data[jm_data$Location=="Birojima",]$pres_abs, na.rm=T) # 330
out[3,]$n_encounter<-sum(jm_data[jm_data$Location=="Hchijojima",]$pres_abs, na.rm=T) # 10
out[4,]$n_encounter<-sum(jm_data[jm_data$Location=="Oki",]$pres_abs, na.rm=T) # 14

write.csv(out, "~/research/miller_et_al/remodelling_2017/Oceanographic_cols.csv", quote=F, row.names=F)
                




