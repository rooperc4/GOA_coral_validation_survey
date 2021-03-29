####IMPORT AND MERGE GOA DATA######

#Loop to make spatial lines out of the survey positions
line.x<-c(GOA_data[1,6],GOA_data[1,7])
line.y<-c(GOA_data[1,4],GOA_data[1,5])
line1<-Line(cbind(line.x,line.y))
line1<-Lines(list(line1),paste(GOA_data[1,1],GOA_data[1,2],GOA_data[1,3]))
Tows<-SpatialLines(list(line1),proj4string=CRS("+proj=longlat +datum=WGS84"))

for(i in 2:length(GOA_data[,1])){
line.x<-c(GOA_data[i,6],GOA_data[i,7])
line.y<-c(GOA_data[i,4],GOA_data[i,5])
line2<-Line(cbind(line.x,line.y))
line2<-Lines(list(line2),paste(GOA_data[i,1],GOA_data[i,2],GOA_data[i,3]))
Hauls1<-SpatialLines(list(line2),proj4string=CRS("+proj=longlat +datum=WGS84"))
Tows<-rbind(Tows,Hauls1)
}

Tows.project<-spTransform(Tows,CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

##Extract Data From Raster Stack
stack.data<-extract(GOA.stack,Tows.project,fun=mean)
GOA_data<-cbind(GOA_data,stack.data)
GOA_data<-subset(GOA_data,GOA_data$depth>0)
GOA_data[,11:28] <-apply(GOA_data[,11:28], 2, function(x){replace(x, is.na(x), 0)})
colnames(GOA_data)[31]<-"depthz"
colnames(GOA_data)[8]<-"depth"


