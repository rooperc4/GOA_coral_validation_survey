haul<-read.csv("C:/Users/rooperc/Desktop/RACEBASE/GOAhaul.csv",header=TRUE,stringsAsFactors=FALSE)
catch<-read.csv("C:/Users/rooperc/Desktop/RACEBASE/GOAcatch.csv",header=TRUE,stringsAsFactors=FALSE)
length<-read.csv("C:/Users/rooperc/Desktop/RACEBASE/GOAlength.csv",header=TRUE,stringsAsFactors=FALSE)


haul<-subset(haul,(haul$REGION=="AI"|haul$REGION=="BS"|haul$REGION=="GOA")&haul$PERFORMANCE>=0&haul$HAUL_TYPE==3&haul$CRUISE>198200&haul$NET_WIDTH>0&haul$DISTANCE_FISHED>0)
haul<-subset(haul,haul$START_LATITUDE>0&haul$END_LATITUDE>0&haul$BOTTOM_DEPTH>0&haul$WIRE_LENGTH>0)
haul$AreaSwept<-as.numeric(haul$NET_WIDTH)*as.numeric(haul$DISTANCE_FISHED)*1000/10000

##############################################################################################
##############INVERTEBRATES###################################################################
invert_species<-read.csv("SpeciesCodes.csv",header=TRUE)
invert_species<-unique(invert_species)
invert_species<-invert_species[-81,]
catch.data<-merge(catch,invert_species,by="SPECIES_CODE",all.x=TRUE)
catch.data<-subset(catch.data,is.na(catch.data$Family)==FALSE)
haul<-merge(haul,catch.data,by="HAULJOIN",all.x=TRUE)

haul<-data.frame(HAULJOIN=haul$HAULJOIN,VESSEL=haul$VESSEL.x,CRUISE=haul$CRUISE.x,HAUL=haul$HAUL.x,START_TIME=haul$START_TIME,DISTANCE_FISHED=haul$DISTANCE_FISHED,NET_WIDTH=haul$NET_WIDTH,START_LATITUDE=as.numeric(haul$START_LATITUDE),END_LATITUDE=as.numeric(haul$END_LATITUDE),
                 START_LONGITUDE=as.numeric(haul$START_LONGITUDE),END_LONGITUDE=as.numeric(haul$END_LONGITUDE),BOTTOM_DEPTH=as.numeric(haul$BOTTOM_DEPTH),GEAR_TEMPERATURE=as.numeric(haul$GEAR_TEMPERATURE),
                 WIRE_LENGTH=as.numeric(haul$WIRE_LENGTH),AreaSwept=haul$AreaSwept,SPECIES_CODE=haul$SPECIES_CODE,Weight=as.numeric(haul$WEIGHT),Family=haul$Family,Taxon=haul$Taxon,Genus=haul$Genus)



haul.pos<-NetPosition(haul$START_LATITUDE,haul$END_LATITUDE,haul$START_LONGITUDE,haul$END_LONGITUDE,haul$WIRE_LENGTH,haul$BOTTOM_DEPTH)
haul$LONGITUDE<-(haul.pos[,3]+haul.pos[,4])/2
haul$LATITUDE<-(haul.pos[,1]+haul.pos[,2])/2

##Remove data from 1998 and 2009 Rockfish Adaptive Sampling Studies
haul<-subset(haul,haul$CRUISE!=199801&haul$CRUISE!=200902)

##Remove data from winter survey in 2001
haul<-subset(haul,haul$CRUISE!=200102)
haul<-subset(haul,haul$CRUISE!=200101|(haul$CRUISE==200101&haul$VESSEL!=143))

haul$CPUE<-haul$Weight/haul$AreaSwept
haul$CPUE[is.na(haul$CPUE)]<-0

write.csv(haul,"GOA_Trawl_data.csv",row.names=FALSE)



#############################################################################
# FUNCTION TO CALCULATE CORRECTED POSITION FOR THE SURVEY TRAWL BEHIND THE  #
# SURVEY VESSEL USING WIRE OUT AND DEPTH RECORDED FROM THE HAUL             #
# Chris Rooper                                                              #
# 2/18/2015                                                                 #
#############################################################################

#START_LAT<-55.58621*pi/180
#END_LAT<-55.57442*pi/180
#START_LON<--168.8596*pi/180
#END_LON<--168.8262*pi/180
#BOTTOM_DEPTH<-1018
#WIRE_OUT<-1803

#NetPosition(START_LAT,END_LAT,START_LON,END_LON,WIRE_OUT,BOTTOM_DEPTH)

NetPosition<-function(START_LAT,END_LAT,START_LON,END_LON,WIRE_OUT,BOTTOM_DEPTH){
START_LAT<-START_LAT*pi/180
END_LAT<-END_LAT*pi/180
START_LON<-START_LON*pi/180
END_LON<-END_LON*pi/180

Bearing<-((atan2(sin(END_LON-START_LON)*cos(END_LAT),cos(START_LAT)*sin(END_LAT)-sin(START_LAT)*cos(END_LAT)*cos(END_LON-START_LON))*180/pi)+360)%%360
Distance_behind<-sqrt(WIRE_OUT^2-BOTTOM_DEPTH^2)*-1

START_LATC<-asin(sin(START_LAT)*cos(Distance_behind/6367449)+cos(START_LAT)*sin(Distance_behind/6367449)*cos(Bearing*pi/180))*180/pi
END_LATC<-asin(sin(END_LAT)*cos(Distance_behind/6367449)+cos(END_LAT)*sin(Distance_behind/6367449)*cos(Bearing*pi/180))*180/pi

START_LONC<-(START_LON+atan2(sin(Bearing*pi/180)*sin(Distance_behind/6367449)*cos(START_LAT),cos(Distance_behind/6367449)-sin(START_LAT)*sin(START_LATC*pi/180)))*180/pi
END_LONC<-(END_LON+atan2(sin(Bearing*pi/180)*sin(Distance_behind/6367449)*cos(END_LAT),cos(Distance_behind/6367449)-sin(END_LAT)*sin(END_LATC*pi/180)))*180/pi

return(data.frame(START_LATC,END_LATC,START_LONC,END_LONC))}
