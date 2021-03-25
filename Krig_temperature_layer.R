
##specify the folder location of the spatial data, not a specific file
library(raster)
library(rgdal)
library(maptools)
library(gstat)
library(rgeos)
library(proj4)
library(sp)

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