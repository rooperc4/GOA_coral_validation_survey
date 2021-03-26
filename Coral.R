## load libraries
library(mgcv)
library(PresenceAbsence)
library(gbm)
library(dismo)
library(rJava)
library(raster)
library(rgdal)
library(maptools)
library(gstat)
library(rgeos)
library(proj4)
library(sp)
library(maptools)
library(maps)
library(xtable)
library(randomForest)
###################MAXENT MODEL TO DETERMINE SUITABLE HABITAT#########################

########STEP 1 - Set up data and directories############

results.path <- "D:/GOA Coral and Sponge Model/Results/Corals"

presence.dat<-subset(GOA_data,GOA_data$Coral_generic>0)
presence.dat<-cbind(presence.dat$lon,presence.dat$lat)
dim(presence.dat)

########STEP 2 - MAKE THE RASTER STACK ############

maxent.stack<-stack(GOA.lat,GOA.bathy,GOA.btemp)
names(maxent.stack)<-c("Latitude","Depth","Temperature")

#######STEP 3 - MAKE THE MODELS#############
#Make the presence only models for winter 
corals.model<-maxent(maxent.stack,presence.dat, args=c("-P","-J"),path=results.path)
plot(corals.model)
response(corals.model)
corals.model

#######STEP 4 - MAKE THE PREDICTION RASTER ##############################################
#Predict the model to a raster to see where suitable habitat is predicted
corals.habitat<-predict(corals.model,maxent.stack,overwrite=TRUE,progress="text",filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsMaxent")
plot(corals.habitat,main="Maxent Model - Probability of suitable habitat")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)
points(presence.dat,col="purple")

######STEP 5 - TEST THE PREDICTIONS AGAINST THE TRAINING DATA############################
#Extract the predictions at the training points
train.predicted<-extract(corals.habitat,presence.dat)

#Create a vector of 1's (presence points for training data)
train.observed<-rep(1,length(train.predicted))

#Choose a random vector of absence points
train.background<-randomPoints(corals.habitat,length(presence.dat[,1]),presence.dat,tryf=50)

#Extract the predictions at these points
train.background.predicted<-extract(corals.habitat,train.background)

#Create a vector of 0's (absence points for training data)
train.background.observed<-rep(0,length(train.background.predicted))

#Bind the observations and predictions together and create a dataframe
train.predicted<-c(train.predicted,train.background.predicted)
train.observed<-c(train.observed,train.background.observed)
train.auc_data<-data.frame(cbind(seq(1,length(train.predicted),1),train.observed,train.predicted))
train.auc_data<-subset(train.auc_data,train.auc_data$train.predicted>=0)

#Calculate the AUC
auc(train.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
optimal.thresholds(train.auc_data,opt.methods=c(seq(1:9)))
maxent.threshold<-optimal.thresholds(train.auc_data,opt.methods=2)
maxent.threshold<-maxent.threshold[,2]
auc.roc.plot(train.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(train.auc_data,N.bins=10)
error.threshold.plot(train.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(train.auc_data,threshold=maxent.threshold)
pcc(cmx(train.auc_data,threshold=maxent.threshold))
sensitivity(cmx(train.auc_data,threshold=maxent.threshold))
specificity(cmx(train.auc_data,threshold=maxent.threshold))
predicted.prevalence(train.auc_data,threshold=maxent.threshold)
presence.absence.accuracy(train.auc_data,threshold=maxent.threshold)
presence.absence.hist(train.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(train.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(train.auc_data[,2],train.auc_data[,3])

##Make a suitable habitat mask
corals.habitat.mask<-cut(corals.habitat,breaks=c(maxent.threshold,1))
writeRaster(corals.habitat.mask,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsMaxentMask",overwrite=TRUE)

######SUBSET DATA WHERE SUITABLE HABITAT IS PRESENT AND CREATE MASK##################################################
corals.pos<-cbind(GOA_data$lon,GOA_data$lat)
corals.keep<-extract(corals.habitat,corals.pos)
corals.data<-subset(GOA_data,corals.keep>maxent.threshold)

corals.data[,"corals_present"]<-corals.data$Coral_generic
corals.data[,"corals_present"]<-ifelse(corals.data[,"corals_present"]>0,1,0)

##Sample training and testing data sets####
coral_indexes = sample(1:nrow(corals.data), size=0.2*nrow(corals.data))
test = corals.data[coral_indexes,]
train = corals.data[-coral_indexes,]

trainc<-subset(train,train[,"corals_present"]>0)
testc<-subset(test,test[,"corals_present"]>0)
(dim(trainc)[1]+dim(testc)[1])/(dim(train)[1]+dim(test)[1])
(dim(train)[1]+dim(test)[1])

############METHOD 1 - DELTA GAM (GAM FOR PRESENCE-ABSENCE, SET THRESHOLD, GAM FOR POSITIVE OCCURANCES #####################################
#fullmodel
corals.DG.pa<-gam(corals_present~s(lon)+s(lat)+s(depth,k=4)+s(btemp,k=4)+s(tmax,k=4)+s(speed,k=4)+s(slope,k=4)+s(color,k=4),data=train,family=binomial)
summary(corals.DG.pa)

#reduced model
corals.DG.pa<-gam(corals_present~s(lon)+s(lat)+s(depth,k=4)+s(btemp,k=4)+s(speed,k=4)+s(slope,k=4)+s(color,k=4),data=train,family=binomial)
summary(corals.DG.pa)
gam.check(corals.DG.pa)
#plot(corals.DG.pa)

#test the predictions on the training data
train_predict<-predict.gam(corals.DG.pa,train,type="response")
train.auc_data<-data.frame(cbind(seq(1,length(train_predict),1),train[,"corals_present"],train_predict))

#Calculate the AUC
auc(train.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
optimal.thresholds(train.auc_data,opt.methods=c(seq(1:9)))
train.threshold<-optimal.thresholds(train.auc_data,opt.methods=6)
train.threshold<-train.threshold[,2]
auc.roc.plot(train.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(train.auc_data,N.bins=10)
error.threshold.plot(train.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(train.auc_data,threshold=train.threshold)
pcc(cmx(train.auc_data,threshold=train.threshold))
sensitivity(cmx(train.auc_data,threshold=train.threshold))
specificity(cmx(train.auc_data,threshold=train.threshold))
predicted.prevalence(train.auc_data,threshold=train.threshold)
presence.absence.accuracy(train.auc_data,threshold=train.threshold)
presence.absence.hist(train.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(train.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(train.auc_data[,2],train.auc_data[,3],method="spearman")

#test the predictions on the test data
test_predict<-predict.gam(corals.DG.pa,test,type="response")
test.auc_data<-data.frame(cbind(seq(1,length(test_predict),1),test[,"corals_present"],test_predict))

#Calculate the AUC
auc(test.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
auc.roc.plot(test.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(test.auc_data,N.bins=10)
error.threshold.plot(test.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(test.auc_data,threshold=train.threshold)
pcc(cmx(test.auc_data,threshold=train.threshold))
sensitivity(cmx(test.auc_data,threshold=train.threshold))
specificity(cmx(test.auc_data,threshold=train.threshold))
predicted.prevalence(test.auc_data,threshold=train.threshold)
presence.absence.accuracy(test.auc_data,threshold=train.threshold)
presence.absence.hist(test.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(test.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(test.auc_data[,2],test.auc_data[,3],method="spearman")

corals.DG.pa.raster<-predict(GOA.stack, corals.DG.pa,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGAMpa",fun=predict, na.rm=TRUE,overwrite=TRUE,progress="text",type="response",newdata.guaranteed=TRUE)
plot(corals.DG.pa.raster,main="Coral Generalized Additive Model - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)
corals.presence.mask<-mask(corals.DG.pa.raster,corals.habitat.mask)
writeRaster(corals.presence.mask,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGAMpa",overwrite=TRUE)

temp_raster2<-corals.habitat.mask
temp_raster1<-cut(corals.presence.mask,breaks=c(0,train.threshold,1))
temp_raster1[temp_raster1==1]<-0
temp_raster1[temp_raster1==2]<-1
temp_raster3<-mosaic(temp_raster1,temp_raster2,fun="sum")
temp_raster3[temp_raster3==1]<-0
temp_raster3[temp_raster3==2]<-1
temp_raster3<-mask(temp_raster3,GOA.bathy)
writeRaster(temp_raster3,overwrite=TRUE,"E://EFH Descriptions 2015/Variables/Variables_GOA_1ha/coralfactor_all")




############GAM Models for Coral abundance
#fullmodel
corals.DG.cpue<-gam((Coral_generic)^.25~s(lon)+s(lat)+s(depth,k=4)+s(btemp,k=4)+s(tmax,k=4)+s(speed,k=4)+s(slope,k=4)+s(color,k=4),data=trainc,family=gaussian)
summary(corals.DG.cpue)

#reduced model
corals.DG.cpue<-gam((Coral_generic)^.25~s(lon)+s(depth,k=4)+s(slope,k=4)+s(color,k=4),data=trainc,family=gaussian)
summary(corals.DG.cpue)
gam.check(corals.DG.cpue)
#plot(corals.DG.cpue,scale=0,residuals=TRUE)

#test the predictions on the training data
observed.train<-(trainc$Coral_generic)^.25
predicted.train<-predict(corals.DG.cpue,trainc,type="response")
pred.train<-lm(observed.train~predicted.train)
summary(pred.train)
plot(observed.train,predicted.train)
abline(pred.train)

#test the predictions on the test data
observed.test<-(testc$Coral_generic)^.25
predicted.test<-predict(corals.DG.cpue,testc,type="response")
pred.test<-lm(observed.test~predicted.test)
summary(pred.test)
plot(observed.test,predicted.test)
abline(pred.test)

temp_raster<-cut(corals.presence.mask,breaks=c(train.threshold,1))
writeRaster(temp_raster,"E://EFH Descriptions 2015/Variables/Variables_GOA_1ha/coralfactor_all")
corals.DG.cpue.raster<-predict(GOA.stack, corals.DG.cpue,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGAMcpue",fun=predict, na.rm=TRUE,overwrite=TRUE,progress="text",type="response",newdata.guaranteed=TRUE)
corals.DG.cpue.raster<-mask(corals.DG.cpue.raster,temp_raster)
writeRaster(corals.DG.cpue.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGAMcpue",overwrite=TRUE)
plot(corals.DG.cpue.raster, main="Coral GAM - PA, GAM for CPUE, Threshold = 0.05")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)

#####OVERALL R2 GAM######
train.observed<-(train$Coral_generic)^0.25
train.predicted<-predict(corals.DG.pa,train,type="response")
train.predicted[train.predicted>=train.threshold]<-1
train.predicted[train.predicted<train.threshold]<-0
train.predicted.cpue<-predict(corals.DG.cpue,train,type="response")
train.predicted.cpue<-train.predicted.cpue*train.predicted
pred.train<-lm(train.predicted.cpue~train.observed)
summary(pred.train)
plot(train.observed,train.predicted.cpue)
abline(pred.train)
mse_coral<-sum(pred.train$residuals^2)/pred.train$df.residual
mse_coral

test.observed<-(test$Coral_generic)^0.25
test.predicted<-predict(corals.DG.pa,test,type="response")
test.predicted[test.predicted>=train.threshold]<-1
test.predicted[test.predicted<train.threshold]<-0
test.predicted.cpue<-predict(corals.DG.cpue,test,type="response")
test.predicted.cpue<-test.predicted.cpue*test.predicted
pred.test<-lm(test.predicted.cpue~test.observed)
summary(pred.test)
plot(test.observed,test.predicted.cpue)
abline(pred.test)
mse_coral<-sum(pred.test$residuals^2)/pred.test$df.residual
mse_coral

############METHOD 2 - Boosted Regression Tree  #####################################
####BRT PRESENCE ABSENCE MODEL######################################

corals.BRT.pa<-gbm.step(data=train, gbm.x = c(8,29:30,32,34:36,38), gbm.y = 39, family = "bernoulli", tree.complexity = 5,learning.rate = 0.005, bag.fraction = 0.5)
gbm.plot(corals.BRT.pa)
gbm.plot.fits(corals.BRT.pa)
summary(corals.BRT.pa)

#test the predictions on the training data
train_predict<-predict.gbm(corals.BRT.pa,train, n.trees=corals.BRT.pa$gbm.call$best.trees, type="response")
train.auc_data<-data.frame(cbind(seq(1,length(train_predict),1),train[,"corals_present"],train_predict))

#Calculate the AUC
auc(train.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
optimal.thresholds(train.auc_data,opt.methods=c(seq(1:9)))
train.threshold<-optimal.thresholds(train.auc_data,opt.methods=2)
train.threshold<-train.threshold[,2]
auc.roc.plot(train.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(train.auc_data,N.bins=10)
error.threshold.plot(train.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(train.auc_data,threshold=train.threshold)
pcc(cmx(train.auc_data,threshold=train.threshold))
sensitivity(cmx(train.auc_data,threshold=train.threshold))
specificity(cmx(train.auc_data,threshold=train.threshold))
predicted.prevalence(train.auc_data,threshold=train.threshold)
presence.absence.accuracy(train.auc_data,threshold=train.threshold)
presence.absence.hist(train.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(train.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(train.auc_data[,2],train.auc_data[,3])

#test the predictions on the test data
test_predict<-predict.gbm(corals.BRT.pa,test, n.trees=corals.BRT.pa$gbm.call$best.trees, type="response")
test.auc_data<-data.frame(cbind(seq(1,length(test_predict),1),test[,"corals_present"],test_predict))

#Calculate the AUC
auc(test.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
auc.roc.plot(test.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(test.auc_data,N.bins=10)
error.threshold.plot(test.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(test.auc_data,threshold=train.threshold)
pcc(cmx(test.auc_data,threshold=train.threshold))
sensitivity(cmx(test.auc_data,threshold=train.threshold))
specificity(cmx(test.auc_data,threshold=train.threshold))
predicted.prevalence(test.auc_data,threshold=train.threshold)
presence.absence.accuracy(test.auc_data,threshold=train.threshold)
presence.absence.hist(test.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(test.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(test.auc_data[,2],test.auc_data[,3])



#Predict to Raster
corals.BRT.pa.raster<-predict(GOA.stack, corals.BRT.pa, n.trees=corals.BRT.pa$gbm.call$best.trees, type="response",filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsBRTpa", na.rm=TRUE,overwrite=TRUE,progress="text",newdata.guaranteed=TRUE)
corals.BRT.pa.raster<-mask(corals.BRT.pa.raster,corals.habitat.mask,overwrite=TRUE)
writeRaster(corals.BRT.pa.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsBRTpa",overwrite=TRUE)
plot(corals.BRT.pa.raster,main="Coral Boosted Regression Tree - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)

#############BRT Models for Coral abundance
train$corals.fourth<-(train$Coral_generic)^0.25
test$corals.fourth<-(test$Coral_generic)^0.25
corals.BRT<-gbm.step(data=train, gbm.x = c(8,29:30,32,34:36,38), gbm.y = 40, family = "gaussian", tree.complexity = 5,learning.rate = 0.01, bag.fraction = 0.5)
gbm.plot(corals.BRT)
gbm.plot.fits(corals.BRT)
summary(corals.BRT)

#Test on training data
train_predict<-predict.gbm(corals.BRT,train, n.trees=corals.BRT$gbm.call$best.trees, type="response")
plot(train_predict~train$corals.fourth)
brt.pred.train<-lm(train_predict~train$corals.fourth)
summary(brt.pred.train)
abline(brt.pred.train)
mse_coral<-sum(brt.pred.train$residuals^2)/brt.pred.train$df.residual
mse_coral

#Test on test data
test_predict<-predict.gbm(corals.BRT,test, n.trees=corals.BRT$gbm.call$best.trees, type="response")
plot(test_predict~test$corals.fourth)
brt.pred.test<-lm(test_predict~test$corals.fourth)
summary(brt.pred.test)
abline(brt.pred.test)
mse_coral<-sum(brt.pred.test$residuals^2)/brt.pred.test$df.residual
mse_coral

#Predict to Raster
corals.BRT.cpue.raster<-predict(GOA.stack, corals.BRT, n.trees=corals.BRT$gbm.call$best.trees, type="response",filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsBRTcpue", na.rm=TRUE,overwrite=TRUE,progress="text",newdata.guaranteed=TRUE)
corals.BRT.cpue.raster<-mask(corals.BRT.cpue.raster,corals.habitat.mask,overwrite=TRUE)
writeRaster(corals.BRT.cpue.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsBRTcpue",overwrite=TRUE)
plot(corals.BRT.cpue.raster,main="Coral Boosted Regression Tree - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)


###########################METHOD 3 - DELTA GLM METHOD#############################################
#fullmodel
corals.GLM.pa<-glm(corals_present~lon+lat+depth+btemp+slope+tmax+speed+color+I(lon^2)++I(lat^2)+I(depth^2)+I(slope^2)+I(btemp^2)+I(tmax^2)+I(speed^2)+I(color^2),data=train,family=binomial)
summary(corals.GLM.pa)

#reduced model
corals.GLM.pa<-glm(corals_present~lon+lat+depth+btemp+slope+tmax+speed+color+I(lon^2)++I(lat^2)+I(depth^2)+I(btemp^2)+I(tmax^2)+I(speed^2)+I(color^2),data=train,family=binomial)
summary(corals.GLM.pa)
#plot(corals.GLM.pa)

#test the predictions on the training data
train_predict<-predict.glm(corals.GLM.pa,train,type="response")
train.auc_data<-data.frame(cbind(seq(1,length(train_predict),1),train[,"corals_present"],train_predict))

#Calculate the AUC
auc(train.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
optimal.thresholds(train.auc_data,opt.methods=c(seq(1:9)))
train.threshold<-optimal.thresholds(train.auc_data,opt.methods=6)
train.threshold<-train.threshold[,2]
auc.roc.plot(train.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(train.auc_data,N.bins=10)
error.threshold.plot(train.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(train.auc_data,threshold=train.threshold)
pcc(cmx(train.auc_data,threshold=train.threshold))
sensitivity(cmx(train.auc_data,threshold=train.threshold))
specificity(cmx(train.auc_data,threshold=train.threshold))
predicted.prevalence(train.auc_data,threshold=train.threshold)
presence.absence.accuracy(train.auc_data,threshold=train.threshold)
presence.absence.hist(train.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(train.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(train.auc_data[,2],train.auc_data[,3],method="spearman")

#test the predictions on the test data
test_predict<-predict.glm(corals.GLM.pa,test,type="response")
test.auc_data<-data.frame(cbind(seq(1,length(test_predict),1),test[,"corals_present"],test_predict))

#Calculate the AUC
auc(test.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
auc.roc.plot(test.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(test.auc_data,N.bins=10)
error.threshold.plot(test.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(test.auc_data,threshold=train.threshold)
pcc(cmx(test.auc_data,threshold=train.threshold))
sensitivity(cmx(test.auc_data,threshold=train.threshold))
specificity(cmx(test.auc_data,threshold=train.threshold))
predicted.prevalence(test.auc_data,threshold=train.threshold)
presence.absence.accuracy(test.auc_data,threshold=train.threshold)
presence.absence.hist(test.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(test.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(test.auc_data[,2],test.auc_data[,3])

corals.GLM.pa.raster<-predict(GOA.stack, corals.GLM.pa,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGLMpa",fun=predict, na.rm=TRUE,overwrite=TRUE,progress="text",type="response",newdata.guaranteed=TRUE)
plot(corals.GLM.pa.raster,main="Coral Generalized Linear Model - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)
corals.GLM.pa.raster<-mask(corals.GLM.pa.raster,corals.habitat.mask)
writeRaster(corals.GLM.pa.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGLMpa",overwrite=TRUE)

#############GLM Models for Coral abundance
#fullmodel
corals.GLM.cpue<-glm((Coral_generic)^.25~lon+lat+depth+btemp+slope+tmax+speed+color+I(lon^2)+I(lat^2)+I(depth^2)+I(slope^2)+I(btemp^2)+I(tmax^2)+I(speed^2)+I(color^2),data=trainc,family=gaussian)
summary(corals.GLM.cpue)

#reduced model
#fullmodel
corals.GLM.cpue<-glm((Coral_generic)^.25~lon+lat+depth+btemp+slope+color+I(lat^2)+I(btemp^2)+I(color^2),data=trainc,family=gaussian)
summary(corals.GLM.cpue)
#plot(corals.GLM.cpue,scale=0,residuals=TRUE)

#test the predictions on the training data
observed.train.GLM<-(trainc$Coral_generic)^.25
predicted.train.GLM<-predict(corals.GLM.cpue,trainc,type="response")
pred.train.GLM<-lm(observed.train.GLM~predicted.train.GLM)
summary(pred.train.GLM)
plot(observed.train.GLM,predicted.train.GLM)
abline(pred.train.GLM)

#test the predictions on the test data
observed.test.GLM<-(testc$Coral_generic)^.25
predicted.test.GLM<-predict(corals.GLM.cpue,testc,type="response")
pred.test.GLM<-lm(observed.test.GLM~predicted.test.GLM)
summary(pred.test.GLM)
plot(observed.test.GLM,predicted.test.GLM)
abline(pred.test.GLM)

temp_raster<-cut(corals.GLM.pa.raster,breaks=c(train.threshold,1))
corals.GLM.cpue.raster<-predict(GOA.stack, corals.GLM.cpue,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGLMcpue",fun=predict, na.rm=TRUE,overwrite=TRUE,progress="text",type="response",newdata.guaranteed=TRUE)
corals.GLM.cpue.raster<-mask(corals.GLM.cpue.raster,temp_raster)
writeRaster(corals.GLM.cpue.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsGLMcpue",overwrite=TRUE)
plot(corals.GLM.cpue.raster, main="Coral GAM - PA, GAM for CPUE, Threshold = 0.05")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)

#####OVERALL R2 GLM######
train.observed<-(train$Coral_generic)^0.25
train.predicted<-predict(corals.GLM.pa,train,type="response")
train.predicted[train.predicted>=train.threshold]<-1
train.predicted[train.predicted<train.threshold]<-0
train.predicted.cpue<-predict(corals.GLM.cpue,train,type="response")
train.predicted.cpue<-train.predicted.cpue*train.predicted
pred.train<-lm(train.predicted.cpue~train.observed)
summary(pred.train)
plot(train.observed,train.predicted.cpue)
abline(pred.train)
mse_coral<-sum(pred.train$residuals^2)/pred.train$df.residual
mse_coral

test.observed<-(test$Coral_generic)^0.25
test.predicted<-predict(corals.GLM.pa,test,type="response")
test.predicted[test.predicted>=train.threshold]<-1
test.predicted[test.predicted<train.threshold]<-0
test.predicted.cpue<-predict(corals.GLM.cpue,test,type="response")
test.predicted.cpue<-test.predicted.cpue*test.predicted
pred.test<-lm(test.predicted.cpue~test.observed)
summary(pred.test)
plot(test.observed,test.predicted.cpue)
abline(pred.test)
mse_coral<-sum(pred.test$residuals^2)/pred.test$df.residual
mse_coral


###########################METHOD 4 - RANDOM FOREST METHOD #########################################
###RF PRESENCE ABSENCE MODEL####

corals.RF.pa<-randomForest(as.factor(corals_present)~lat+lon+depth+slope+btemp+speed+tmax+color,data=train, ntree=1000,importance=TRUE, na.action=na.omit)
print(corals.RF.pa)
round(importance(corals.RF.pa), 2)

train_predict_RF<-predict(corals.RF.pa,train,type="prob")[,2]
train.auc_data<-data.frame(cbind(seq(1,length(train_predict_RF),1),train[,"corals_present"],train_predict_RF))

#Calculate the AUC
auc(train.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
optimal.thresholds(train.auc_data,opt.methods=c(seq(1:9)))
train.threshold<-optimal.thresholds(train.auc_data,opt.methods=2)
train.threshold<-train.threshold[,2]
auc.roc.plot(train.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(train.auc_data,N.bins=10)
error.threshold.plot(train.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(train.auc_data,threshold=train.threshold)
pcc(cmx(train.auc_data,threshold=train.threshold))
sensitivity(cmx(train.auc_data,threshold=train.threshold))
specificity(cmx(train.auc_data,threshold=train.threshold))
predicted.prevalence(train.auc_data,threshold=train.threshold)
presence.absence.accuracy(train.auc_data,threshold=train.threshold)
presence.absence.hist(train.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(train.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(train.auc_data[,2],train.auc_data[,3])


#Test on test data
test_predict_RF<-predict(corals.RF.pa,test,type="prob")[,2]
test.auc_data<-data.frame(cbind(seq(1,length(test_predict_RF),1),test[,"corals_present"],test_predict_RF))

#Calculate the AUC
auc(test.auc_data,na.rm=TRUE)

#Estimate the thresholds and calculate diagnostics
auc.roc.plot(test.auc_data,opt.methods=c(seq(1:9)))
calibration.plot(test.auc_data,N.bins=10)
error.threshold.plot(test.auc_data,opt.methods=c(seq(1:9)),opt.thresholds=TRUE)
cmx(test.auc_data,threshold=train.threshold)
pcc(cmx(test.auc_data,threshold=train.threshold))
sensitivity(cmx(test.auc_data,threshold=train.threshold))
specificity(cmx(test.auc_data,threshold=train.threshold))
predicted.prevalence(test.auc_data,threshold=train.threshold)
presence.absence.accuracy(test.auc_data,threshold=train.threshold)
presence.absence.hist(test.auc_data,color=c("green","red"),truncate.tallest=TRUE)
presence.absence.summary(test.auc_data,opt.methods=c(2,4,5),N.bins=10,truncate.tallest=TRUE)
cor.test(test.auc_data[,2],test.auc_data[,3])


#Predict to Raster
corals.RF.pa.raster<-predict(GOA.stack, corals.RF.pa,type="prob",filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsRFpa", na.rm=TRUE,overwrite=TRUE,progress="text",newdata.guaranteed=TRUE)
corals.RF.pa.raster<-mask(corals.RF.pa.raster,corals.habitat.mask,overwrite=TRUE)
writeRaster(corals.RF.pa.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsRFpa",overwrite=TRUE)
plot(corals.RF.pa.raster,main="Coral Boosted Regression Tree - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)

#####RF CPUE MODEL

corals.RF<-randomForest(corals.fourth~lat+lon+depth+slope+btemp+speed+tmax+color,data=train, ntree=1000,importance=TRUE, na.action=na.omit)
print(corals.RF)
round(importance(corals.RF), 2)

#Test on training data
train_predict_RF<-predict(corals.RF,train,type="response")
plot(train_predict_RF~train$corals.fourth)
RF.pred.train<-lm(train_predict_RF~train$corals.fourth)
summary(RF.pred.train)
abline(RF.pred.train)
mse_coral<-sum(RF.pred.train$residuals^2)/RF.pred.train$df.residual
mse_coral

#Test on test data
test_predict_RF<-predict(corals.RF,test,type="response")
plot(test_predict_RF~test$corals.fourth)
RF.pred.test<-lm(test_predict_RF~test$corals.fourth)
summary(RF.pred.test)
abline(RF.pred.test)
mse_coral<-sum(RF.pred.test$residuals^2)/RF.pred.test$df.residual
mse_coral

#Predict to Raster
corals.RF.cpue.raster<-predict(GOA.stack, corals.RF,type="response",filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsRFcpue", na.rm=TRUE,overwrite=TRUE,progress="text",newdata.guaranteed=TRUE)
corals.RF.cpue.raster<-mask(corals.RF.cpue.raster,corals.habitat.mask,overwrite=TRUE)
writeRaster(corals.RF.cpue.raster,filename="D:/GOA Coral and Sponge Model/Results/PredictionRasters/coralsRFcpue",overwrite=TRUE)
plot(corals.RF.cpue.raster,main="Coral Boosted Regression Tree - PA")
plot(akland,col=1,add=TRUE)
plot(canada_land,col=1,add=TRUE)

removeTmpFiles(h=.25)































