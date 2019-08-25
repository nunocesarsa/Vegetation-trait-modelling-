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

#investigate
library(randomForestSRC)

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

#BROAD LIMITS (NOT THE SAME AS THE CITATIONS) - i am avoiding absence of any trait
param.maxmin <- matrix(c(0.2, 1.5, #leaf layers or leaf structure
                         1,80, #Cab
                         1,40, #Car
                         0.005,0.04, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.04, #Cm
                         0.1,8), #LAI
                         #0.05,0.2), #hotspot
                       nrow=6,ncol = 2,byrow = T)

#to ensure that no traits are outsides the test space, we will reduce the bit the limits of the surface
param.maxmin.valid <- param.maxmin
param.maxmin.valid[,1] <- param.maxmin.valid[,1]+param.maxmin.valid[,1]*.1
param.maxmin.valid[,2] <- param.maxmin.valid[,2]-param.maxmin.valid[,2]*.1

#first we look at our ability of predicting the traits in sentinel data and then also each of the indices
#creating a training space
train.n <- 1000 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
train.LHS <- Latinhyper(param.maxmin,train.n)

valid.n <- 1000 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin.valid,valid.n)

#checking the stats
summary(train.LHS)
summary(valid.LHS)

#lets build the dataset and train the model
train.trait.df <- data.frame(N=train.LHS[,1],
                             Cab=train.LHS[,2],
                             Car=train.LHS[,3],
                             Cw=train.LHS[,4],
                             Cm=train.LHS[,5],
                             LAI=train.LHS[,6],
                             #lidfa=train.LHS[,7],
                             #TypeLidf = 0 ,
                             tts = 0,
                             tto = 30)

valid.trait.df <- data.frame(N=valid.LHS[,1],
                             Cab=valid.LHS[,2],
                             Car=valid.LHS[,3],
                             Cw=valid.LHS[,4],
                             Cm=valid.LHS[,5],
                             LAI=valid.LHS[,6],
                             #lidfa=valid.LHS[,7],
                             #TypeLidf = 0 ,
                             tts = 0,
                             tto = 30)

#there are no NA here 
summary(train.trait.df)
summary(valid.trait.df)

head(train.trait.df)

#creating the spectral libraries from PROSAIL RTM
train.spclib <- PROSAIL(parameterList = train.trait.df)
train.spectr <- spectra(train.spclib)
train.spectr.df <- as.data.frame(train.spectr)

valid.spclib <- PROSAIL(parameterList = valid.trait.df)
valid.spectr <- spectra(valid.spclib)
valid.spectr.df <- as.data.frame(valid.spectr)


par(mfrow=c(1,2))
plot(train.spclib)
plot(valid.spclib)
#now we create a set for SENTINEL spectral libraries and select only the bands measured at 20m
#something is going wrong in this transformation, NA are being created, so they have to be removed

s2.train.spclib <- spectralResampling(train.spclib,"Sentinel2",response_function = T)
s2.train.spectr <- spectra(s2.train.spclib)

s2.train.spclib.20m <- s2.train.spclib[,c(2,3,4,5,6,7,9,12,13)] #only the 20m bands
s2.train.spectr.20m <- spectra(s2.train.spclib.20m)

s2.valid.spclib <- spectralResampling(valid.spclib,"Sentinel2",response_function = TRUE)
s2.valid.spectr <- spectra(s2.valid.spclib)

s2.valid.spclib.20m <- s2.valid.spclib[,c(2,3,4,5,6,7,9,12,13)]
s2.valid.spectr.20m <- spectra(s2.valid.spclib.20m)

#some final touches
s2.train.spectr.20m.df <- as.data.frame(s2.train.spectr.20m)
names(s2.train.spectr.20m.df) <- c("B02","B03","B04",
                                   "B05","B06","B07",
                                   "B8A","B11","B12")
s2.valid.spectr.20m.df <- as.data.frame(s2.valid.spectr.20m)
names(s2.valid.spectr.20m.df) <- c("B02","B03","B04",
                                   "B05","B06","B07",
                                   "B8A","B11","B12")

summary(s2.train.spectr.20m.df)
summary(s2.valid.spectr.20m.df)
#check if there are NA.. this is a prboelm, i think NA are being given when the values go under 0 in reflectance
#ATTENTION the changes in lines here have to be reflected in the validation data
names(s2.train.spectr.20m.df)
#s2.train.spectr.20m.df <- s2.train.spectr.20m.df[,-c(8,9)]
#s2.valid.spectr.20m.df <- s2.valid.spectr.20m.df[,-c(8,9)]

summary(s2.train.spectr.20m.df)
summary(s2.valid.spectr.20m.df)

lin.index.train <- which(is.na(s2.train.spectr.20m.df))
lin.index.valid <- which(is.na(s2.valid.spectr.20m.df))

#and now we clean all data
s2.train.spectr.20m.df.clean <- s2.valid.spectr.20m.df[-lin.index.train,]
s2.valid.spectr.20m.df.clean <- s2.valid.spectr.20m.df[-lin.index.valid,]
dim(s2.train.spectr.20m.df.clean)
dim(s2.valid.spectr.20m.df.clean)

train.trait.df.clean <- train.trait.df[-lin.index.train,]
valid.trait.df.clean <- valid.trait.df[-lin.index.valid,]

dim(train.trait.df.clean)
names(train.trait.df.clean)

#now we bring it into a DF
names(train.trait.df)
s2.mRF.train.df <- cbind(train.trait.df.clean[,-c(7,8)],s2.train.spectr.20m.df.clean) #we remove the fixed parameters
s2.mRF.valid.df <- cbind(valid.trait.df.clean[,-c(7,8)],s2.valid.spectr.20m.df.clean) #we remove the fixed parameters

names(s2.mRF.train.df)
head(s2.mRF.train.df)
#let's go first the single target RF 
#case 7 params
N.train.df   <- s2.mRF.train.df[,-c(2,3,4,5,6,7)]
Cab.train.df <- s2.mRF.train.df[,-c(1,3,4,5,6,7)]
Car.train.df <- s2.mRF.train.df[,-c(1,2,4,5,6,7)]
Cw.train.df  <- s2.mRF.train.df[,-c(1,2,3,5,6,7)]
Cm.train.df  <- s2.mRF.train.df[,-c(1,2,3,4,6,7)]
LAI.train.df <- s2.mRF.train.df[,-c(1,2,3,4,5,7)]
#ALA.train.df <- s2.mRF.train.df[,-c(1,2,3,4,5,6)] #lidfa is equivalent to ALA so we will keep it like that, helps organizing the code

#case of 6 params
N.train.df   <- s2.mRF.train.df[,-c(2:6)]
Cab.train.df <- s2.mRF.train.df[,-c(1,3:6)]
Car.train.df <- s2.mRF.train.df[,-c(1:2,4:6)]
Cw.train.df  <- s2.mRF.train.df[,-c(1:3,5:6)]
Cm.train.df  <- s2.mRF.train.df[,-c(1:4,6)]
LAI.train.df <- s2.mRF.train.df[,-c(1:5)]

#verify this
names(N.train.df)
names(Cab.train.df)
names(Car.train.df)
names(Cw.train.df)
names(Cm.train.df)
names(LAI.train.df)


st.mRF.N    <- rfsrc(N~.,data=N.train.df,block.size = 50)
st.mRF.Cab  <- rfsrc(Cab~.,data=Cab.train.df,block.size = 50)
st.mRF.Car  <- rfsrc(Car~.,data=Car.train.df,block.size = 50)
st.mRF.Cw   <- rfsrc(Cw~.,data=Cw.train.df,block.size = 50)
st.mRF.Cm   <- rfsrc(Cm~.,data=Cm.train.df,block.size = 50)
st.mRF.LAI  <- rfsrc(LAI~.,data=LAI.train.df,block.size = 50)
#st.mRF.ALA  <- rfsrc(lidfa~.,data=ALA.train.df,block.size = 50)

#let's see how accurate are these single trait prediction models
names(valid.trait.df.clean)
st.test.df <- valid.trait.df.clean[,-c(7,8)]
dim(s2.valid.spectr.20m.df.clean)

st.test.df$stPred_N   <- get.mv.predicted(predict(st.mRF.N,newdata=s2.valid.spectr.20m.df.clean))
st.test.df$stPred_Cab <- get.mv.predicted(predict(st.mRF.Cab,newdata=s2.valid.spectr.20m.df.clean))
st.test.df$stPred_Car <- get.mv.predicted(predict(st.mRF.Car,newdata=s2.valid.spectr.20m.df.clean))
st.test.df$stPred_Cw  <- get.mv.predicted(predict(st.mRF.Cw,newdata=s2.valid.spectr.20m.df.clean))
st.test.df$stPred_Cm  <- get.mv.predicted(predict(st.mRF.Cm,newdata=s2.valid.spectr.20m.df.clean))
st.test.df$stPred_LAI <- get.mv.predicted(predict(st.mRF.LAI,newdata=s2.valid.spectr.20m.df.clean))
#st.test.df$stPred_ALA <- get.mv.predicted(predict(st.mRF.ALA,newdata=s2.valid.spectr.20m.df))

head(st.test.df)

#lets plot these
names(valid.trait.df.clean)
temp.df <- data.frame(trait=names(valid.trait.df.clean[,-c(7,8)]),
                      slope=NA,inter=NA,Rsqur=NA,
                      RMSE=NA,MAE=NA,MAPE=NA)
temp.df
#separting the prediction from the table is actually useful now..
#names(st.test.df)
st.test.pred.df <- st.test.df[,7:ncol(st.test.df)]

rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }

head(valid.trait.df.clean)
head(st.test.pred.df)

par(mfrow=c(2,3))
for (i in 1:nrow(temp.df)){
  
  plot(valid.trait.df.clean[,i], st.test.pred.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  #linear model stats
  x <- valid.trait.df.clean[,i]
  y <- st.test.pred.df[,i]
  
  trait.fit <- lm(y~x)
  abline(trait.fit,lty=2,lwd=2)
  
  temp.df$slope[i] <- trait.fit$coefficients[[2]]
  temp.df$inter[i] <- trait.fit$coefficients[[1]]
  temp.df$Rsqur[i] <- summary(trait.fit)$r.square
  temp.df$RMSE[i]  <- RMSE.custom(trait.fit$residuals)
  temp.df$MAE[i]   <- mae.custom(trait.fit$residuals)
  temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

temp.df$Mean <- colMeans(s2.mRF.valid.df[,1:6])
temp.df$std  <- colSds(as.matrix(s2.mRF.valid.df[,1:6]))
temp.df

#plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#mtext("Single target RF", at=0.2, cex=1)

head(s2.mRF.train.df)

#MULTIOUTPUT Random Forest
names(s2.mRF.train.df)
names(s2.mRF.train.df)
dim(s2.mRF.train.df)
mRF.Block50 <- rfsrc(Multivar(N,Cab,Car,
                              Cw,Cm,LAI)~.,data=s2.mRF.train.df,block.size = 50)
                     #ensemble = "oob",
                     #bootstrap="by.root",samptype = "swr",
                     #xvar.wt = c(N.5,.25,.25,.5,.40,.45),
                     #importance="none")
dim(s2.mRF.train.df)
names(s2.mRF.train.df)

mRF.test.pred.df <- as.data.frame(get.mv.predicted(predict(mRF.Block50,
                                                           newdata=s2.valid.spectr.20m.df.clean)))

#we can now plot the outputs

mRF.temp.df <- data.frame(trait=names(valid.trait.df[,-c(7,8)]),
                      slope=NA,inter=NA,Rsqur=NA,RMSE=NA,MAE=NA,MAPE=NA)

summary(mRF.Block50)
#mRF.Block50$ensemble

dim(s2.valid.spectr.20m.df.clean)
dim(mRF.test.pred.df)
dim(mRF.test.pred.df)

names(mRF.test.pred.df)
names(s2.valid.spectr.20m.df.clean)

par(mfrow=c(2,3))
for (i in 1:nrow(temp.df)){
  
  plot(mRF.test.pred.df[,i],s2.mRF.valid.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  #linear model stats
  x <- s2.mRF.valid.df[,i]
  y <- mRF.test.pred.df[,i]
  
  mRF.trait.fit <- lm(y~x)
  abline(mRF.trait.fit,lty=2,lwd=2)
  
  mRF.temp.df$slope[i] <- mRF.trait.fit$coefficients[[2]]
  mRF.temp.df$inter[i] <- mRF.trait.fit$coefficients[[1]]
  mRF.temp.df$Rsqur[i] <- summary(mRF.trait.fit)$r.square
  mRF.temp.df$RMSE[i]  <- RMSE.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAE[i]   <- mae.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
mtext("Multi target RF", at=0.2, cex=1)