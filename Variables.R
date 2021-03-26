################ BASE AREA RASTER AND BATHYMETRY AND SLOPE VARIABLES ###########################

##Making data in R from ArcGIS raster layers
##specify the folder location of the spatial data, not a specific file
library(raster)
library(rgdal)
library(maptools)
library(gstat)
library(rgeos)
library(proj4)
library(sp)

#import the data from arcgis raster layer
GOAbathyz<-raster("//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/Bathy/egoa_grid")
GOAbathyz<-projectRaster(GOAbathyz, crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs",overwrite=TRUE)
breaks<-c(0,1000)
GOA1000<-cut(GOAbathyz,breaks=breaks, filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/GOA1000",progress="text",overwrite=TRUE)

#Bring in Alaska land for maps
akland<-readShapeSpatial("V:/GIS stuff/Rstuff/AKland", CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
canada_land<-readShapeSpatial("V:/GIS stuff/Rstuff/canada", CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

GOA.bathy<-mask(GOAbathyz,GOA1000,overwrite=TRUE,progress="text")
GOA.bathy<-mask(GOA.bathy,akland,overwrite=TRUE,inverse=TRUE, progress="text")

seamount.rm.x<-c(0,500000,1000000,1600000,0)
seamount.rm.y<-c(300000,800000,750000,300000,300000)
seamount.rm.pol<-Polygon(cbind(seamount.rm.x,seamount.rm.y))
seamount.rm.pol<-Polygons(list(seamount.rm.pol),"Seamount")
seamount.rm.pol<-SpatialPolygons(list(seamount.rm.pol), proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
plot(seamount.rm.pol,add=TRUE)

ebs.rm.x<-c(-500000,-500000,-150000,-150000,-500000)
ebs.rm.y<-c(620000,1200000,1200000,820000,620000)
ebs.rm.pol<-Polygon(cbind(ebs.rm.x,ebs.rm.y))
ebs.rm.pol<-Polygons(list(ebs.rm.pol),"ebs")
ebs.rm.pol<-SpatialPolygons(list(ebs.rm.pol), proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
plot(ebs.rm.pol,add=TRUE)

GOA.bathy<-mask(GOA.bathy,seamount.rm.pol,overwrite=TRUE,inverse=TRUE, progress="text")
GOA.bathy<-mask(GOA.bathy,ebs.rm.pol,overwrite=TRUE,inverse=TRUE, progress="text",file="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goabathy")

GOA.slope<-terrain(GOA.bathy, opt="slope",unit="degrees", neighbors=8, filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goaslope") 
GOA.aspect<-terrain(GOA.bathy, opt="aspect",unit="degrees", neighbors=8,overwrite=TRUE, filename="D:/GOA Coral and Sponge Model/Variables/RasterGrids/aspect") 
GOA.rugosity<-terrain(GOA.bathy, opt="TRI", neighbors=8,overwrite=TRUE, filename="D:/GOA Coral and Sponge Model/Variables/RasterGrids/rugosity") 
GOA.tpi<-terrain(GOA.bathy, opt="TPI", neighbors=8,overwrite=TRUE, filename="D:/GOA Coral and Sponge Model/Variables/RasterGrids/TPI") 


######## MEAN CURRENT VARIABLE ########
#IDW to interpolate currents at points and make currents into interpolated raster
#Get NEP current data and project
currents<-read.table("//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/Deep Sea Coral Research/Distribution model/Variables/NEP Current/speedbot_ave.asc",header=TRUE)
currents<-subset(currents,currents$Speed>0)
currents[,1]<-currents[,1]-360
current1<-cbind(currents[,1],currents[,2])
currents.project <- project(current1, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
currents.pos<-data.frame(cbind(currents.project[,1],currents.project[,2]))
currents.pos<-SpatialPoints(currents.pos, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 
currents.project<-data.frame(cbind(currents.project,currents[,3]))
#currents.select<-over(currents.pos,area[1])
#currents.project<-cbind(currents.project,currents.select)
#currents.project<-subset(currents.project,currents.project[,4]==1)
#currents.project<-currents.project[,-4]
currents.project<-SpatialPointsDataFrame(coords=c(currents.project["X1"],currents.project["X2"]),data=currents.project, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

#Interpolate to blank raster
GOA.speed.idw<-gstat(id = "X3", formula = X3~1, data=currents.project, nmax=7, set=list(idp = 2.5))
GOA.speed<-interpolate(GOA.bathy,GOA.speed.idw,overwrite=TRUE,xyOnly=TRUE, progress="text")
GOA.speed<-mask(GOA.speed,GOA.bathy,filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goaspeed",overwrite=TRUE)
plot(GOA.speed)

####### SEA WIFS VARIABLE #############
seawifs1<-array(0,dim=c(0,3))
month<-c("may","june","july","aug","sept")
year<-c(seq(2003,2011,1))
for(j in 1:9){
	for(i in 1:5){
path<-paste("V:/EBS Canyons Project/SeaWIFS data/",month[i],year[j],".xyz",sep="")
wifs1<-cbind(read.table(path,skip=1),0)
wifs1<-subset(wifs1,wifs1[,2]>=51&wifs1[,2]<=62&wifs1[,3]>0)
seawifs1<-rbind(seawifs1,wifs1)}

wifsmean<-aggregate(seawifs1,by=list(seawifs1[,1],seawifs1[,2]),mean)
wifscount<-aggregate(seawifs1,by=list(seawifs1[,1],seawifs1[,2]),sum)/wifsmean
seawifs2<-cbind(wifsmean[,3:5],wifscount[,5])
}

wifsmean1<-aggregate(seawifs2,by=list(seawifs2[,1],seawifs2[,2]),mean)
wifscount1<-aggregate(seawifs2,by=list(seawifs2[,1],seawifs2[,2]),sum)
seawifs<-cbind(wifsmean1[,3:5],wifscount1[,5])

seawif1<-cbind(seawifs[,1],seawifs[,2])
seawifs.project <- project(seawif1, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
seawifs.pos<-data.frame(cbind(seawifs.project[,1],seawifs.project[,2]))
seawifs.pos<-SpatialPoints(seawifs.pos, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 
seawifs.project<-data.frame(cbind(seawifs.project,seawifs[,3]))
seawifs.project<-SpatialPointsDataFrame(coords=c(seawifs.project["X1"],seawifs.project["X2"]),data=seawifs.project, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

#Interpolate to blank raster
GOA.colork.idw<-gstat(id = "X3", formula = X3~1, data=seawifs.project, nmax=7, set=list(idp = 2.5))
GOA.color<-interpolate(GOA.bathy,GOA.colork.idw,overwrite=TRUE,xyOnly=TRUE, progress="text")
GOA.color<-mask(GOA.color,GOA.bathy,filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goacolor",overwrite=TRUE)

####### LATITUDE AND LONGITUDE VARIABLE #######
## Make a grid of x and y from lat long##

latt <- init(GOA.bathy, v='y')
GOA.lat<-mask(latt,GOA.bathy,filename="D://GOA Coral and Sponge Model/Variables/RasterGrids/goalat",overwrite=TRUE)

lont <- init(GOA.bathy, v='x')
GOA.lon<-mask(lont,GOA.bathy,filename="D://GOA Coral and Sponge Model/Variables/RasterGrids/goalon",overwrite=TRUE)

######TIDAL CURRENT VARIABLE#########
##Make grid of points to calculate tides on###
#GOA.tide1<-GOA.bathy
#res(GOA.tide1)<-1000
#test1<-xyFromCell(GOA.tide1, seq(0:ncell(GOA.tide1)),spatial=TRUE)
#test2<-extract(GOA.bathy,test1)
#test1<-data.frame(test1)
#test2<-data.frame(test2)
#tide.pts<-cbind(test1,test2)
#tide.pts<-subset(tide.pts,tide.pts[,3]>0)
#tide.pts<-tide.pts[,1:2]
tide.pts<-cbind(GOA_hauls["Start_longC"],GOA_hauls["Start_latC"])
#tide.pts<-SpatialPoints(tide.pts, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 
#tide.pts<- spTransform(tide.pts, CRS("+proj=longlat +datum=WGS84"))
tide.pts<-data.frame(tide.pts)

####Connecting Matlab to R to Run Tidal Inversion Model######
#Make the connection for R and Matlab
#Open Matlab
#Set path to TMD
#type "MatlabServer" at prompt

#In R
library(R.matlab)
library(R.utils)
matlab<-Matlab()
is.open<-open(matlab)
print(matlab)

#Make time file in matlab to compute tides at one hour intervals for 1 year + 4 days (starting 1/1/2009)
evaluate(matlab,"SDtime=[floor(datenum([2009 1 1 1 0 0])):1/24:floor(datenum([2009 1 1 1 0 0]))+368]")

#Extract latitude and longitude of stations with positions to predict tides at
lat<-tide.pts[,2]
lon<-tide.pts[,1]
max_tide<-rep(0,length(lat))
sd_tide<-rep(0,length(lat))
mean_tide<-rep(0,length(lat))
median_tide<-rep(0,length(lat))


####START LOOP AT 377###############
#loop the tidal prediction over the positions and the year + 4 days and calculate the maximum
for(i in 377:length(lat)){
text_u<-paste("[TS,ConList]=tmd_tide_pred('DATA/Model_AleI2013',SDtime",lat[i],lon[i],"'u')",sep=",")
text_v<-paste("[TS,ConList]=tmd_tide_pred('DATA/Model_AleI2013',SDtime",lat[i],lon[i],"'v')",sep=",")

#RUN TIDAL MODEL FOR u and v
evaluate(matlab, text_u)
u <- getVariable(matlab, c("TS"))
evaluate(matlab, text_v)
v <- getVariable(matlab, c("TS"))

max_tide[i]<-max((u$TS^2+v$TS^2)^(1/2))
sd_tide[i]<-sd((u$TS^2+v$TS^2)^(1/2))
mean_tide[i]<-mean((u$TS^2+v$TS^2)^(1/2))
median_tide[i]<-median((u$TS^2+v$TS^2)^(1/2))
evaluate(matlab, "fclose('all')")
print(i)
}

tmax<-cbind(lon,lat,max_tide[1:7653],sd_tide[1:7653],mean_tide[1:7653],median_tide[1:7653])
colnames(tmax)<-c("lon","lat","max_tide","sd_tide","mean_tide","median_tide")
evaluate(matlab, "quit()")

tmax<-subset(tmax,tmax[,3]>0)
tide1<-cbind(tmax[,1],tmax[,2])
tmax.project <- project(tide1, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
tmax.project<-data.frame(cbind(tmax.project,tmax[,3]))
tmax.project<-SpatialPointsDataFrame(coords=c(tmax.project["X1"],tmax.project["X2"]),data=tmax.project, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

## kriging from gstat
v <- variogram(X3~1, tmax.project)
m <- fit.variogram(v, vgm(.2625, "Exp", 200000, .06))
g <- gstat(NULL, "X3", X3~1, tmax.project, model = m, nmax=50)
tmax.raster.kr <- interpolate(GOA.bathy, g, xyOnly=TRUE, progress="text",overwrite=TRUE )
GOA.tmax<-mask(tmax.raster.kr,GOA.bathy,filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goatmax",overwrite=TRUE)

x <- krige.cv(X3~1, tmax.project, m, nmax = 12)
bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")
plot(x$observed,x$var1.pred)
summary(lm(x$var1.pred~x$observed))
mean(x$residual^2)

tmax<-cbind(lon,lat,max_tide[1:7653],sd_tide[1:7653],mean_tide[1:7653],median_tide[1:7653])
GOA_hauls<-cbind(GOA_hauls,tmax[,3])

######## BOTTOM TEMPERATURE VARIABLE ########
#Temperature kriging to raster
tempdata<-cbind(GOA_hauls$Start_longC,GOA_hauls$Start_latC,GOA_hauls$GEAR_TEMPERATURE)
tempdata<-subset(tempdata,tempdata[,3]>=0)
temp1<-cbind(tempdata[,1],tempdata[,2])
tempdata.project <- project(temp1, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
tempdata.pos<-data.frame(cbind(tempdata.project[,1],tempdata.project[,2]))
tempdata.pos<-SpatialPoints(tempdata.pos, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 
tempdata.project<-data.frame(cbind(tempdata.project,tempdata[,3]))
tempdata.project<-SpatialPointsDataFrame(coords=c(tempdata.project["X1"],tempdata.project["X2"]),data=tempdata.project, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

## kriging from gstat
v <- variogram(X3~1, tempdata.project)
m <- fit.variogram(v, vgm(.2625, "Exp", 200000, .06))
plot(v,model=m)
g <- gstat(NULL, "X3", X3~1, tempdata.project, model = m, nmax=50)

x <- krige.cv(X3~1, tempdata.project, m, nmax = 12)
bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")
plot(x$observed,x$var1.pred)
summary(lm(x$var1.pred~x$observed))
mean(x$residual^2)

btemp.raster.kr <- interpolate(GOA.bathy, g, xyOnly=TRUE, progress="text",overwrite=TRUE )
GOA.btemp<-mask(btemp.raster.kr,GOA.bathy,filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goabtemp",overwrite=TRUE)

##############ASPECT VARIABLE##############################
GOA.aspect.project<-terrain(GOA.bathy,opt="aspect",unit="degrees",neighbors=8,"//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goaaspect",overwrite=TRUE,progress="text")

cdirection<-read.table("//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/Deep Sea Coral Research/Distribution model/Variables/Aspect/uvbot_ave_with_label.asc",header=TRUE)
cdirection<-subset(cdirection,cdirection$speed>0)
cdirection[,1]<-cdirection[,1]-360
current2<-cbind(cdirection[,1],cdirection[,2])
cdirection.project <- project(current2, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
cdirection.pos<-data.frame(cbind(cdirection.project[,1],cdirection.project[,2]))
cdirection.pos<-SpatialPoints(cdirection.pos, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

cangle<-rep(-9999,107814)
for(i in 1:107814){
cangle[i]<-circle_angle(cdirection[i,4],cdirection[i,3])}

##FUNCTION TO DO CONVERSION TO 360 DEGREES###
#circle_angle<-function(rise,run){
#if(rise>=0&run>=0){
#	ca1<-(90-atan(rise/run)*180/pi)}
#if(rise<0&run>=0){
#	ca1<-abs(atan(rise/run)*180/pi)+90}
#if(rise<0&run<0){
#	ca1<-(360-(90+atan(rise/run)*180/pi))}
#if(rise>=0&run<0){
#	ca1<-(-(atan(rise/run)*180/pi)+270)}
#	return(ca1)}

cdirection.project<-data.frame(cbind(cdirection.project,cangle))
colnames(cdirection.project)<-c("X1","X2","X3")
#Make Polygons of GOA area
GOA.bb<-bbox(GOA.bathy)
goa.x<-c(GOA.bb[1,1],GOA.bb[1,2],GOA.bb[1,2],GOA.bb[1,1],GOA.bb[1,1])
goa.y<-c(GOA.bb[2,2],GOA.bb[2,2],GOA.bb[2,1],GOA.bb[2,1],GOA.bb[2,2])
goa.pol<-Polygon(cbind(goa.x,goa.y))
goa.pol<-Polygons(list(goa.pol),"goa")
goa.pol<-SpatialPolygons(list(goa.pol), proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

cdirection.select<-over(cdirection.pos,goa.pol)
cdirection.project<-cbind(cdirection.project,cdirection.select)
cdirection.project<-subset(cdirection.project,cdirection.project[,4]==1)
cdirection.project<-cdirection.project[,-4]
cdirection.project<-SpatialPointsDataFrame(coords=c(cdirection.project["X1"],cdirection.project["X2"]),data=cdirection.project, proj4string=CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) 

#Interpolate to blank raster
GOA.cdirect.idw<-gstat(id = "X3", formula = X3~1, data=cdirection.project, nmax=7, set=list(idp = 2.5))
GOA.cdirect.raster<-interpolate(GOA.bathy,GOA.cdirect.idw,overwrite=TRUE,xyOnly=TRUE, progress="text")
plot(GOA.cdirect.raster)

GOA.align<-abs(GOA.cdirect.raster-GOA.aspect)
fun <- function(x) { x[x>180] <- abs(x-360); return(x) }
GOA.aspecta <- calc(GOA.align, fun,filename="//nmfs.local/AKC-RACE/Users/chris.rooper/Desktop/GOA Coral and Sponge Model/Variables/RasterGrids/goaaspecta",overwrite=TRUE, progress="text")

plot(GOA.aspecta)


###### RASTER STACK #######
GOA.stack<-stack(GOA.lon,GOA.lat,GOA.bathy,GOA.slope,GOA.rugosity,GOA.btemp,GOA.color,GOA.speed,GOA.aspecta,GOA.tmax)
names(GOA.stack)<-c("lon","lat","depth","slope","rugosity","btemp","color","speed","aspect","tmax")

