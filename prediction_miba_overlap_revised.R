
## japanese marine iba network coverage of jm spatial prediction ##

rm(list=ls())
setwd("~/research/miller_et_al/")


#libraries
library(raster)
library(rgeos)
library(rgdal)
library(maptools)

#crop out Korean and conflict zone jm predictions

cols<-read.csv("input_data/JM_colony_Japan_12_12_12_japan_comments.csv", h=T)

cols<-cols[cols$Location!="Daegugul Island" & cols$Location!="Dok Island/Takeshima Island", ] #removes korean and conflict zone colonies

sp_cols<-SpatialPoints(data.frame(cols$Longitude, cols$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))

jap_miba<-readOGR( layer="japan_mIBA_network", dsn="GIS", verbose=TRUE)

NAGASHIMA<-jap_miba[grep("Nagashima", jap_miba$NatName),]
NAGASHIMA@polygons[[1]]@ID<-"99"
row.names(NAGASHIMA@data)<-"99"

jap_miba<-jap_miba[which(jap_miba$NatName!="NA"),] # selects only marine ibas, not islands

jap_miba<-spRbind(jap_miba, NAGASHIMA)

DgProj <- CRS("+proj=laea +lon_0=134 +lat_0=34")

TshapeProj <- spTransform(sp_cols, CRS=DgProj)
IBAProj <- spTransform(jap_miba, CRS=DgProj)

plot(IBAProj)
plot(TshapeProj, add=T, border=2)


for(i in 1:length(sp_cols)) # loop yoinked from extract_raster_making.r
{
  
  

  TBuffProj <- gBuffer(TshapeProj,byid=T, width=8000, quadsegs=50) #note the new buffer = max range from data
  TBuffProj <-SpatialPolygonsDataFrame(TBuffProj, data=cols)
  TBuffProj$area<-gArea(TBuffProj, byid=T )
  
  int1<-gIntersection(TBuffProj,IBAProj )
  
  # need to figure out differences due to land being cut out of
  # miba layer but not buffer layer. Raster might be the best approach
  # what did i do originally,maybe some code to jump on there
  # doesnt have to be exact just approximate
  
  plot(TBuffProj)
  plot(TshapeProj, add=T, col="grey")
  TBuffProj@polygons[[1]]@ID <- as.character(i)
  TBuffWgs <- spTransform(TBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))#from ras
  col_inf<-rasterize(TBuffWgs, jm_pred, field=1, background=0)
  if(i == 1) {col_fin <- col_inf} else {col_fin <- col_fin + col_inf}
  print(i)
} 

col_fin[values(col_fin>0)]<-1
col_fin[values(col_fin==0)]<-NA


pred_clip<-jm_pred*col_fin #clips prediction to max foraging range (48km!!!!)

#normalization as per Lavers paper
pred_norm<- (pred_clip - (min(values(pred_clip), na.rm=T)))/(max(values(pred_clip), na.rm=T)-min(values(pred_clip), na.rm=T))

writeRaster(pred_norm, "results/Ensemble_pred_clip_norm.tif", overwrite=TRUE)



jap_miba<-readOGR( layer="japan_mIBA_network", dsn="GIS", verbose=TRUE)

NAGASHIMA<-jap_miba[grep("Nagashima", jap_miba$NatName),]

jap_miba<-jap_miba[grep("- Marine", jap_miba$NatName),] # selects only marine ibas, not islands

jap_miba<-spRbind(jap_miba, NAGASHIMA)

ras<-raster("extract_rasters/dist_land_1km.tif")

plot(pred_norm);plot(jap_miba, add=T)

ovrlp_out<-NULL

for(j in c(0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
{
  
   jm_class_pred<-pred_norm
   jm_class_pred[values(jm_class_pred>=j)]<-1
   jm_class_pred[values(jm_class_pred<j)]<-0 
   
   miba_ras<-rasterize(jap_miba, jm_class_pred, field=1, background=0)
   
   miba_overlap<-jm_class_pred+miba_ras
    
  #length(which((values(jm_class_pred)==1))) # area of avail JM breeding habitat at pred level
  #length(which((values(miba_overlap)==2))) # area of avail JM breeding habitat at pred level within miba network
  
   ovrlp<-round(length(which((values(miba_overlap)==2)))/length(which((values(jm_class_pred)==1)))*100)
   area_miba_olp<-length(which(values(miba_overlap)==2))
   area_jm_pred_class<-length(which(values(jm_class_pred)==1))
   
   ovrlp_out<-rbind(ovrlp_out, data.frame(level=j, overlap_perc=ovrlp, area_jm_pred_class=area_jm_pred_class, area_miba_olp=area_miba_olp))
  
   print(j)
}


#current japanese miba network covers..


write.csv(ovrlp_out , "results/overlap_all_miba_CORRECT.csv", quote=F, row.names=F)

olp<-ggplot(ovrlp_out[ovrlp_out$level>0,],
            aes(x=level, y=overlap_perc)) +geom_point()+geom_line()+ scale_x_continuous(
          breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))  + scale_y_continuous(
            limits=c(0, 100)) + xlab("Core level of JM probability of occurrence")  + ylab("% overlap with Japanese mIBA network")+theme_bw()

png("D:/BIRDLIFE/miller_et_al/results/pred_miba_overlap.png", width = 8, height =6 , units ="in", res =600)

olp

dev.off()



## Appendix analyses of % miba contribution to JM protection

jm_pred<-raster("results/Ensemble_pred_clip_norm.tif")

pred_mask<-jm_pred
pred_mask[values(pred_mask>=0)]<-1 # created to remove land from both rasters in loop

jap_miba<-readOGR( layer="japan_mIBA_network", dsn="GIS", verbose=TRUE)

NAGASHIMA<-jap_miba[grep("Nagashima", jap_miba$NatName),]

jap_miba<-jap_miba[grep("- Marine", jap_miba$NatName),] # selects only marine ibas, not islands

jap_miba<-spRbind(jap_miba, NAGASHIMA)

ras<-raster("extract_rasters/dist_land_1km.tif")

plot(jm_pred);plot(jap_miba, add=T)

overlap_all<-data.frame(core_levels=c("all", 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
for( i in 1: length(jap_miba))
{
  j1<-jap_miba[i,]
  jr1<-rasterize(j1, ras, field=1, background=0)
  jr1<-jr1*pred_mask
  #miba_area<-length(which((values(jr1)>0))) # area of miba km^2. This is not absulte correct as it has been altered to
  # match that of the prediction, which relies on satellite imagery that doent not cover all inshore areas and match coastline
  if(is.na(max(values(jr1), na.rm=T))){print(i);next}
  
  jr1_pred<-jr1+jm_pred
  
 
  miba_overlap_levels=c(  
                          length(which(values(jr1_pred)>=1)),
                          length(which(values(jr1_pred)>=1.1)),
                          length(which(values(jr1_pred)>=1.2)),
                          length(which(values(jr1_pred)>=1.3)),
                          length(which(values(jr1_pred)>=1.4)),
                          length(which(values(jr1_pred)>=1.5)),
                          length(which(values(jr1_pred)>=1.6)),
                          length(which(values(jr1_pred)>=1.7)),
                          length(which(values(jr1_pred)>=1.8)),
                          length(which(values(jr1_pred)>=1.9)))
  
  overlap_all<-cbind(overlap_all, name=miba_overlap_levels)  
  
  names(overlap_all)[names(overlap_all)=="name"]<-as.character(j1$SitRecID)
  print(overlap_all)
  }

write.csv(overlap_all , "results/overlap_detail_miba_CORRECT.csv", quote=F, row.names=F)

write.csv(data.frame(t(overlap_all)), "results/overlap_detail_miba_trans_CORRECT.csv", quote=F, row.names=T)

## Appendix analyses of the % coverage of increasing important habitat
## that each miba contains i.e. How much of the miba's area is covered
## by habitat quality X

jm_pred<-raster("results/Ensemble_pred_clip_norm.tif")

pred_mask<-jm_pred
pred_mask[values(pred_mask>=0)]<-1 # created to remove land from both rasters in loop

jap_miba<-readOGR( layer="japan_mIBA_network", dsn="GIS", verbose=TRUE)

jap_miba<-jap_miba[grep("- Marine", jap_miba$NatName),] # selects only marine ibas, not islands

ras<-raster("extract_rasters/dist_land_1km.tif")

plot(jm_pred);plot(jap_miba, add=T)

overlap_all<-NULL
for( i in 1: length(jap_miba))
    {
      j1<-jap_miba[i,]
      jr1<-rasterize(j1, ras, field=1, background=0)
      jr1<-jr1*pred_mask
      #miba_area<-length(which((values(jr1)>0))) # area of miba km^2. This is not absulte correct as it has been altered to
      # match that of the prediction, which relies on satellite imagery that doent not cover all inshore areas and match coastline
      if(is.na(max(values(jr1), na.rm=T))){print(i);next}
      
      jr1_pred<-jr1+jm_pred
      
      overlap_out<-data.frame(Name=jap_miba[i,]$NatName, SitRecID=jap_miba[i,]$SitRecID, abs_area=jap_miba[i,]$SitArea_1,
                              miba_rast_area=length(which((values(jr1)>0))), 
                              overlap_all=length(which((values(jr1_pred)>1)))/length(which((values(jr1)>0)))*100,
                              overlap_0.1_core=round(length(which((values(jr1_pred)>=1.1)))/length(which((values(jr1)>0)))*100),
                              overlap_0.2_core=round(length(which((values(jr1_pred)>=1.2)))/length(which((values(jr1)>0)))*100),
                              overlap_0.3_core=round(length(which((values(jr1_pred)>=1.3)))/length(which((values(jr1)>0)))*100),
                              overlap_0.4_core=round(length(which((values(jr1_pred)>=1.4)))/length(which((values(jr1)>0)))*100),
                              overlap_0.5_core=round(length(which((values(jr1_pred)>=1.5)))/length(which((values(jr1)>0)))*100),
                              overlap_0.6_core=round(length(which((values(jr1_pred)>=1.6)))/length(which((values(jr1)>0)))*100),
                              overlap_0.7_core=round(length(which((values(jr1_pred)>=1.7)))/length(which((values(jr1)>0)))*100),
                              overlap_0.8_core=round(length(which((values(jr1_pred)>=1.8)))/length(which((values(jr1)>0)))*100),
                              overlap_0.9_core=round(length(which((values(jr1_pred)>=1.9)))/length(which((values(jr1)>0)))*100))
      
      overlap_all<-rbind(overlap_out, overlap_all)                        
      
    }

write.csv(overlap_all , "results/overlap_pred_miba_CORRECT.csv", quote=F, row.names=F)

library(ggplot2)

qplot(y=as.vector(colSums(overlap_all[,5:length(overlap_all)], na.rm=T)), 
                  x=seq(1,10), geom="point")

no_na_overlap<-t(na.omit(overlap_all))





