
#the overall conclusion is that ANN is super good for the problem while mRF fails. Try it!


#library for PROSAIL and also Spectral Angle Mapper
library(hsdar)

#library that has some extra random generation tools
library(MCMCglmm)

#Latin hypercube comes from here
library(FME)

#for more options on generating random samples
library(MCMCglmm)

#the single target classifier
library(randomForest)

#machine learning librarys
library(randomForestSRC)
library(ANN2)


#GIS/RS
library(maptools)
library(rgeos)
library(rgdal)
library(raster)
library(sp)


#for std col stats
library(matrixStats)

#setwd
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"
dir.create(dump.fld) #it gives out a warning if the folder exits
set.seed(1000)


#N Structure parameter
#Cab chlorophyll content
#Car Carotenoid content
#Cbrown Brown pigment content
#Cw Equivalent water thickness
#Cm Dry matter content
#psoil Dry/Wet soil factor
#LAI Leaf area index

#TypeLidf Type of leaf angle distribution. See details section
#lidfa Leaf angle distribution. See details section
#lidfb Leaf angle distribution. See details section
#hspot Hotspot parameter
#tts Solar zenith angle
#tto Observer zenith angle
#psi Relative azimuth angle
#parameterList An optional object of class 'data.frame'. Function will iterate over rows of parameterList setting missing entries to default values. See examples section.
#rsoil background (soil) reflectance. N

#these parameters com from here
#https://www.sciencedirect.com/science/article/pii/S0303243415000100
#also of interest: https://www.mdpi.com/2072-4292/11/13/1597/htm
#car limits taken from: #limits taken from: https://www.mdpi.com/2072-4292/8/2/119

param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  10,80, #Cab
  5,40, #Car
  0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.01,0.02, #Cm
  0.5,7),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

#param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
#  5,100, #Cab
#  5,40, #Car
#  0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
#  0.01,0.02, #Cm
#  0.5,8),#LAI
#  #0.05,0.1), #hotstop
#  nrow=5,ncol = 2,byrow = T)


#creating a training space
train.n <- 150 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
train.LHS <- Latinhyper(param.maxmin,train.n)

#to force having the limits we will add something oo the train
train.grid <- Grid(param.maxmin,train.n)
train.LHS <- rbind(train.LHS,train.grid)

valid.n <- 2* nrow(train.LHS) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin,valid.n)


#checking the stats
summary(train.LHS)
summary(valid.LHS)

#sentinel position
sun_zenith = 25.8734
obs_zenith = 5.4748
rel_azimut = 281.0692 - 5.4748


#lets build the dataset and train the model
train.trait.df <- data.frame(#N=train.LHS[,1],
  Cab=train.LHS[,1],
  Car=train.LHS[,2],
  Cw=train.LHS[,3],
  Cm=train.LHS[,4],
  LAI=train.LHS[,5],
  #hspot=train.LHS[,7],
  #TypeLidf = 0 ,
  tts = sun_zenith,
  tto = obs_zenith,
  psi = rel_azimut)

valid.trait.df <- data.frame(#N=valid.LHS[,1],
  Cab=valid.LHS[,1],
  Car=valid.LHS[,2],
  Cw=valid.LHS[,3],
  Cm=valid.LHS[,4],
  LAI=valid.LHS[,5],
  #hspot=valid.LHS[,7],
  #TypeLidf = 0 ,
  tts = sun_zenith,
  tto = obs_zenith,
  psi = rel_azimut)

#there are no NA here 
summary(train.trait.df)
summary(valid.trait.df)

#creating the spectral libraries from PROSAIL RTM
train.spclib <- PROSAIL(parameterList = train.trait.df)
train.spectr <- spectra(train.spclib)
train.spectr.df <- as.data.frame(train.spectr)

valid.spclib <- PROSAIL(parameterList = valid.trait.df)
valid.spectr <- spectra(valid.spclib)
valid.spectr.df <- as.data.frame(valid.spectr)


#ressampling to sentinel
s2.train.spclib <- spectralResampling(train.spclib,"Sentinel2",response_function = T)
s2.train.spectr <- spectra(s2.train.spclib)

s2.train.spclib.20m <- s2.train.spclib[,c(2,3,4,5,6,7,9,12,13)] #only the 20m bands
s2.train.spectr.20m <- spectra(s2.train.spclib.20m)

s2.valid.spclib <- spectralResampling(valid.spclib,"Sentinel2",response_function = TRUE)
s2.valid.spectr <- spectra(s2.valid.spclib)

s2.valid.spclib.20m <- s2.valid.spclib[,c(2,3,4,5,6,7,9,12,13)]
s2.valid.spectr.20m <- spectra(s2.valid.spclib.20m)

par(mfrow=c(1,2))
plot(s2.train.spclib)
plot(s2.valid.spclib)

#check for NA's.. sometimes it has been happening
summary(s2.train.spectr.20m)
summary(s2.valid.spectr.20m)


#some final touches
s2.train.spectr.20m.df <- as.data.frame(s2.train.spectr.20m)
names(s2.train.spectr.20m.df) <- c("B02","B03","B04",
                                   "B05","B06","B07",
                                   "B8A","B11","B12")
s2.valid.spectr.20m.df <- as.data.frame(s2.valid.spectr.20m)
names(s2.valid.spectr.20m.df) <- c("B02","B03","B04",
                                   "B05","B06","B07",
                                   "B8A","B11","B12")

head(s2.train.spectr.20m.df)
head(s2.valid.spectr.20m.df)


#now, lets do multi target ML -> single target is for later and in another file
#now we bring it into a DF
names(train.trait.df)
s2.mRF.train.df <- cbind(train.trait.df[,-c(6:8)],s2.train.spectr.20m.df) #we remove the fixed parameters
s2.mRF.valid.df <- cbind(valid.trait.df[,-c(6:8)],s2.valid.spectr.20m.df) #we remove the fixed parameters
names(s2.mRF.train.df)
names(s2.mRF.valid.df)

#training and predicting the mRF
mRF.Block50 <- rfsrc(Multivar(Cab,Car,
                              Cw,Cm,LAI)~.,data=s2.mRF.train.df,block.size = 50)

mRF.test.pred.df <- as.data.frame(get.mv.predicted(predict(mRF.Block50,
                                                           newdata=s2.valid.spectr.20m.df)))



#lets train and predict a ANN: 
#Architecture: 10,6 (1 more for each param and 1 more than the outputs)
#the first layer activation is a relu function while the second is a hyperbolic tangent
#first we need to create matrices from our data
head(s2.mRF.train.df)
Y.mat <- as.matrix(s2.mRF.train.df[,c(1:5)])
X.mat <- as.matrix(s2.mRF.train.df[,c(6:14)])
X.mat.valid <- as.matrix(s2.mRF.valid.df[,c(6:14)])

nn.5000 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,hidden.layers =c(10,6),loss.type = "squared",
                         activ.functions = c("relu",'tanh'),n.epochs = 10000,standardize = T)

nn.pred.5000 <- predict(nn.5000,X.mat.valid)

#and a df for the fit summary
mRF.temp.df <- data.frame(trait=names(valid.trait.df[,-c(6:8)]),
                          slope=NA,inter=NA,Rsqur=NA,RMSE=NA,MAE=NA,MAPE=NA)

ann.temp.df <- data.frame(trait=names(valid.trait.df[,-c(6:8)]),
                          slope=NA,inter=NA,Rsqur=NA,RMSE=NA,MAE=NA,MAPE=NA)

rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }


#plotting the MRF
par(mfrow=c(2,3))
for (i in 1:nrow(mRF.temp.df)){
  
  #linear model stats
  x <- valid.trait.df[,i]
  y <- mRF.test.pred.df[,i]
  
  plot(y,x,
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=mRF.temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  mRF.trait.fit <- lm(y~x)
  abline(mRF.trait.fit,lty=2,lwd=2)
  
  mRF.temp.df$slope[i] <- mRF.trait.fit$coefficients[[2]]
  mRF.temp.df$inter[i] <- mRF.trait.fit$coefficients[[1]]
  mRF.temp.df$Rsqur[i] <- summary(mRF.trait.fit)$r.square
  mRF.temp.df$RMSE[i]  <- RMSE.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAE[i]   <- mae.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

#these are the same for all datasets
mRF.temp.df$Mean <- colMeans(s2.mRF.valid.df[,1:5])
mRF.temp.df$std  <- colSds(as.matrix(s2.mRF.valid.df[,1:5]))
mRF.temp.df

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
mtext("Multi objective RF",cex=1.5)


### ploting the ANN
par(mfrow=c(2,3))
for (i in 1:nrow(mRF.temp.df)){
  
  #linear model stats
  x <- valid.trait.df[,i]
  y <- nn.pred.5000$predictions[,i]
  
  plot(y,x,
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=mRF.temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  ann.trait.fit <- lm(y~x)
  abline(ann.trait.fit,lty=2,lwd=2)
  
  ann.temp.df$slope[i] <- ann.trait.fit$coefficients[[2]]
  ann.temp.df$inter[i] <- ann.trait.fit$coefficients[[1]]
  ann.temp.df$Rsqur[i] <- summary(ann.trait.fit)$r.square
  ann.temp.df$RMSE[i]  <- RMSE.custom(ann.trait.fit$residuals)
  ann.temp.df$MAE[i]   <- mae.custom(ann.trait.fit$residuals)
  ann.temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

#these are the same for all datasets
ann.temp.df$Mean <- colMeans(s2.mRF.valid.df[,1:5])
ann.temp.df$std  <- colSds(as.matrix(s2.mRF.valid.df[,1:5]))
ann.temp.df

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
mtext("Multi objective ANN",cex=1.5)


nn.5000$Rcpp_ANN
summary(nn.5000)
plot(nn.5000)

ann.temp.df
mRF.temp.df

## can we plot our neural network - probably yes but i can just draw a scheme in powerpoint faster
write.csv(ann.temp.df,paste(dump.fld,
                       "Inverse_mANN_5Params_Summary.csv",sep="/"))

write.csv(mRF.temp.df,paste(dump.fld,
                            "Inverse_mRF_5Params_Summary.csv",sep="/"))
