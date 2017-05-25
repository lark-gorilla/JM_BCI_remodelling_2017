
## japanese marine iba network coverage of jm spatial prediction ##

rm(list=ls())
setwd("~/research/miller_et_al/")


#libraries
library(raster)
library(rgeos)
library(rgdal)
library(maptools)

land_mask<-raster("GIS/land_mask_wgs.tif")
land_mask[is.na(values(land_mask)),]<-0
land_mask[values(land_mask)==1,]<-NA


#crop out Korean and conflict zone jm predictions

cols<-read.csv("input_data/JM_colony_Japan_12_12_12_japan_comments.csv", h=T)

cols<-cols[cols$Location!="Daegugul Island" & cols$Location!="Dok Island/Takeshima Island", ] #removes korean and conflict zone colonies

sp_cols<-SpatialPoints(data.frame(cols$Longitude, cols$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

jap_miba<-readOGR( layer="japan_mIBA_network_torishima_northwest_joined", dsn="GIS", verbose=TRUE)

jap_miba<-jap_miba[which(jap_miba$NatName!="NA"),] # selects only marine ibas, not islands

# drop "Kamogawa estuary" as not around a JM colony
jap_miba<-jap_miba[jap_miba$NatName!="Kamogawa estuary",] 


plot(land_mask)
plot(jap_miba, add=T)
plot(sp_cols, add=T)

miba_overlap<-NULL
for ( i in unique(jap_miba$NatName))
{
  miba<-jap_miba[jap_miba$NatName==i,]
  
  ## remove holes (islands from polygon)
  #https://gis.stackexchange.com/questions/224048/deleting-inner-holes-rings-borders-of-a-spatial-polygon-in-r
  ring = SpatialPolygons(
    list(Polygons(list(miba@polygons[[1]]@Polygons[[1]]),ID=1)),
    proj4string=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  miba<-SpatialPolygonsDataFrame(ring,data=miba@data, match.ID=FALSE)
  ####
  
  in_cols<-sp_cols[miba,] # select cols inside miba poly

  miba_ras<-rasterize(miba, land_mask, field=1, background=0)

  Tshape <- in_cols
  DgProj <- CRS(paste("+proj=laea +lon_0=", Tshape@coords[1], " +lat_0=", Tshape@coords[2], sep=""))
  TshapeProj <- spTransform(Tshape, CRS=DgProj)
  TBuffProj8 <- gBuffer(TshapeProj, width=8000, quadsegs=50)
  TBuffProj11 <- gBuffer(TshapeProj, width=11000, quadsegs=50)
  TBuffProj30 <- gBuffer(TshapeProj, width=30000, quadsegs=50)
  TBuffProj37 <- gBuffer(TshapeProj, width=37000, quadsegs=50)
 
  TBuffWgs8 <- spTransform(TBuffProj8, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  TBuffWgs11 <- spTransform(TBuffProj11, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  TBuffWgs30 <- spTransform(TBuffProj30, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#from ras
  TBuffWgs37 <- spTransform(TBuffProj37, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#from ras
  
  col_8km<-rasterize(TBuffWgs8, land_mask, field=1, background=0)
  col_11km<-rasterize(TBuffWgs11, land_mask, field=1, background=0)
  col_30km<-rasterize(TBuffWgs30, land_mask, field=1, background=0)
  col_37km<-rasterize(TBuffWgs37, land_mask, field=1, background=0)
  
  cov8km<-col_8km+land_mask
  cov11km<-col_11km+land_mask
  cov30km<-col_30km+land_mask
  cov37km<-col_37km+land_mask
  miba_masked<-miba_ras+land_mask
  
  out<-data.frame(MIBA=i, miba_size=length(miba_masked[values(miba_masked)==1]),
                  buf8_size= length(cov8km[values(cov8km)==1]),
                  buf11_size= length(cov11km[values(cov11km)==1]),
                  buf30_size= length(cov30km[values(cov30km)==1]),
                  buf37_size= length(cov37km[values(cov37km)==1]))
  
  miba_overlap<-rbind(out, miba_overlap)
  print(i)
  }
 
write.csv(miba_overlap,"~/research/miller_et_al/remodelling_2017/miba_overlap.csv",
          quote=F, row.names=F)                 
   
  
## make buffer spatial layer for map
out<-NULL
for ( i in 1:length(sp_cols))
{

  Tshape <- sp_cols[i,]
  DgProj <- CRS(paste("+proj=laea +lon_0=", Tshape@coords[1], " +lat_0=", Tshape@coords[2], sep=""))
  TshapeProj <- spTransform(Tshape, CRS=DgProj)
  
  TBuffProj37 <- gBuffer(TshapeProj, width=37000, quadsegs=50)
  TBuffProj37@polygons[[1]]@ID <- as.character(paste(i, "37"))
  TBuffProj30 <- gBuffer(TshapeProj, width=30000, quadsegs=50)
  TBuffProj30@polygons[[1]]@ID <- as.character(paste(i, "30"))
  b1<-spRbind(TBuffProj37, TBuffProj30)
  TBuffProj11 <- gBuffer(TshapeProj, width=11000, quadsegs=50)
  TBuffProj11@polygons[[1]]@ID <- as.character(paste(i, "11"))
  b1<-spRbind(b1, TBuffProj11)
  TBuffProj8 <- gBuffer(TshapeProj, width=8000, quadsegs=50)
  TBuffProj8@polygons[[1]]@ID <- as.character(paste(i, "8"))
  b1<-spRbind(b1, TBuffProj8)
  
  b1Wgs <- spTransform(b1, CRS=CRS( "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  if(i==1){out<-b1Wgs}else{out<-spRbind(out, b1Wgs)}
  plot(out)
  print(i)
}

df1<-data.frame(ID=sort(rep(seq(1,22), 4)), size=rep(c("37", "30", "11", "8"), 22))
row.names(df1)<-paste(df1$ID, df1$size)

out<-SpatialPolygonsDataFrame(out, df1)

writeOGR(out, layer="cols_diff_radii_buff", dsn="GIS", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)


