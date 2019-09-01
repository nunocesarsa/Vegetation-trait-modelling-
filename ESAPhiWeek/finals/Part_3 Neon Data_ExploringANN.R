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
library(neuralnet)

#library for optimization procedure - you still should read what this actually does lol
library(SoilHyP)

#GIS/RS
library(maptools)
library(rgeos)
library(rgdal)
library(raster)
library(sp)


#for std col stats
library(matrixStats)

#for kfolding
library(dismo)

#for the rmse
library(Metrics)

#for handling data frame more efficiently
library(reshape2)

#group summary
library(nlme)

#for raster stuff
library(raster)
library(maptools)

#for mean absolute percentage error
library(DescTools)

#useful functions for later
rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }


#setwd
setwd("D:/ESAPhiWeek/")
gc()
dump.neon.results <- "./Out_Neon"
dump.ml.opti.fld <- "./Out_OptimVsML"
dump.ml.fld <- "./Out_MachL"
dump.fld <- "./Out_Optim"

dir.create(dump.neon.results) #it gives out a warning if the folder exits
#setting seed
set.seed(2000)

param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  1,50, #Cab
  1,50, #Car
  0.01,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.01,0.03, #Cm
  0.5,2.5),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

param.maxmin.train <- param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  .5,60, #Cab
  .5,60, #Car
  0.005,0.035, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.005,0.035, #Cm
  0.25,3),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)



#creating a training space
train.n <- 10000
#train.grid <- Grid(param.maxmin,train.n)
#train.LHS <- Latinhyper(param.maxmin,train.n-nrow(train.grid))
train.LHS <- Latinhyper(param.maxmin.train,train.n)

#to force having the limits we will add something oo the train
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

#now, lets do multi target ML
#now we bring it into a DF
names(train.trait.df)
s2.mRF.train.df <- cbind(train.trait.df[,-c(6:8)],s2.train.spectr.20m.df) #we remove the fixed parameters
s2.mRF.valid.df <- cbind(valid.trait.df[,-c(6:8)],s2.valid.spectr.20m.df) #we remove the fixed parameters
names(s2.mRF.train.df)
names(s2.mRF.valid.df)


#now we train our model and predict against our simulated validation data
#from there we will have some inference of our model ability
#once that is done, we proceed to work on NEON data

#lets train and predict a ANN: 
#Architecture: 10,6 (1 more for each param and 1 more than the outputs)
#the first layer activation is a relu function while the second is a hyperbolic tangent
#first we need to create matrices from our data
head(s2.mRF.train.df)
Y.mat <- as.matrix(s2.mRF.train.df[,c(1:5)])
X.mat <- as.matrix(s2.mRF.train.df[,c(6:14)])
X.mat.valid <- as.matrix(s2.mRF.valid.df[,c(6:14)])

#lets make all the values between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01.c <- function(x,xmin,xmax){(x-xmin)/(xmax-xmin)}
range01.inv <- function(x,xmin,xmax){(x*(xmax-xmin)+xmin)}

Y.mat.n <- Y.mat

param.maxmin
head(Y.mat.n)
Y.mat.n[,1] <- range01.c(Y.mat.n[,1],param.maxmin.train[1,1]-.5,param.maxmin.train[1,2]+.5)
Y.mat.n[,2] <- range01.c(Y.mat.n[,2],param.maxmin.train[2,1]-.5,param.maxmin.train[2,2]+.5)
Y.mat.n[,3] <- range01.c(Y.mat.n[,3],param.maxmin.train[3,1]-.005,param.maxmin.train[3,2]+.005)
Y.mat.n[,4] <- range01.c(Y.mat.n[,4],param.maxmin.train[4,1]-.005,param.maxmin.train[4,2]+.005)
Y.mat.n[,5] <- range01.c(Y.mat.n[,5],param.maxmin.train[5,1]-.5,param.maxmin.train[5,2]+.5)
summary(Y.mat.n)

library(scales)
X.mat.s <- as.matrix(apply(X.mat,2,rescale))
X.mat.valid.s <- as.matrix(apply(X.mat.valid,2,rescale))

summary(X.mat.s)
summary(X.mat.valid.s)


#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c("tanh","relu")
lrates <- 0.005
epochs <- 500
optype <- "adam"

ann.3000 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,
                            hidden.layers =net.st,
                            loss.type = "squared",
                            activ.functions = act.fn,
                            n.epochs = epochs,
                            standardize = T,
                            learn.rates = lrates,
                            optim.type = optype)



ann.3000.n <- neuralnetwork(X=X.mat,y=Y.mat.n,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype)
#plot(ann.3000)
#plot(ann.3000.n)
#with grid  0.172064
#without grid < 0.1



#lets check how our model looks likes
ann.3000.pred <- predict(ann.3000,X.mat.valid)
ann.3000.pred.df <- as.data.frame(ann.3000.pred)

#lets check how our model looks likes
ann.3000.n.pred <- predict(ann.3000.n,X.mat.valid.s)
ann.3000.n.pred.df <- as.data.frame(ann.3000.n.pred)


names(ann.3000.pred.df) <- paste("ann",names(valid.trait.df[,-c(6:8)]),sep="_")
names(ann.3000.n.pred.df) <- names(ann.3000.pred.df) 

summary(ann.3000.n.pred.df)
summary(ann.3000.pred$predictions)

#visual model check up
par(mfrow=c(1,5))
model.df <- cbind(valid.trait.df[,-c(6:8)],ann.3000.pred.df)
plot(model.df$Cab,model.df$ann_Cab)
plot(model.df$Car,model.df$ann_Car)
plot(model.df$Cw,model.df$ann_Cw)
plot(model.df$Cm,model.df$ann_Cm)
plot(model.df$LAI,model.df$ann_LAI)


ann.3000.n.pred.df$ann_Cab <- range01.inv(ann.3000.n.pred.df$ann_Cab ,param.maxmin.train[1,1]-.5,param.maxmin.train[1,2]+.5)
ann.3000.n.pred.df$ann_Car <- range01.inv(ann.3000.n.pred.df$ann_Car ,param.maxmin.train[2,1]-.5,param.maxmin.train[2,2]+.5)
ann.3000.n.pred.df$ann_Cw <- range01.inv(ann.3000.n.pred.df$ann_Cw ,param.maxmin.train[3,1]-.005,param.maxmin.train[3,2]+.005)
ann.3000.n.pred.df$ann_Cm <- range01.inv(ann.3000.n.pred.df$ann_Cm ,param.maxmin.train[4,1]-.005,param.maxmin.train[4,2]+.005)
ann.3000.n.pred.df$ann_LAI <- range01.inv(ann.3000.n.pred.df$ann_LAI ,param.maxmin.train[5,1]-.5,param.maxmin.train[5,2]+.5)


par(mfrow=c(1,5))
plot(model.df$Cab,ann.3000.n.pred.df$ann_Cab)
plot(model.df$Car,ann.3000.n.pred.df$ann_Car)
plot(model.df$Cw,ann.3000.n.pred.df$ann_Cw)
plot(model.df$Cm,ann.3000.n.pred.df$ann_Cm)
plot(model.df$LAI,ann.3000.n.pred.df$ann_LAI)


#final.shp <- readShapePoints(paste(dump.neon.results,"NEON_SCE_and_ANN_and_ESA_xy.shp",sep="/"))
final.csv <- read.csv(paste(dump.neon.results,"NEON_SCE_and_ANN_and_ESA_xy.csv",sep="/"))


#bands
names(final.csv)
final.bands <- final.csv[,c(12:20)]
head(final.bands)
ann.3000.pred.bands <- predict(ann.3000,as.matrix(final.bands))
ann.3000.pred.bands.n <- predict(ann.3000.n,as.matrix(final.bands))

par(mfrow=c(1,1))
plot(final.csv$LAI,ann.3000.pred.bands$predictions[,5])
plot(final.csv$LAI,range01.inv(ann.3000.pred.bands.n$predictions[,5],param.maxmin[5,1]-.5,param.maxmin[5,2]+.5))

plot(final.csv$ESA_LAI_cab,ann.3000.pred.bands$predictions[,1]*ann.3000.pred.bands$predictions[,5])

wavelength(s2.valid.spclib.20m)

ann.3000.pred.bands.df <- as.data.frame(ann.3000.pred.bands)
names(ann.3000.pred.bands.df )
names(ann.3000.pred.bands.df ) <- paste("ann2",names(valid.trait.df[,-c(6:8)]),sep="_")


final.csv.new <- cbind(final.csv,ann.3000.pred.bands.df)

summary(final.csv.new)
head(final.csv.new)

write.csv(final.csv.new,paste(dump.neon.results,"NEON_SCE_and_ANN_and_ESA_xy_new.csv",sep="/"))

library(ggplot2)

p.all <- ggplot(final.csv.new, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "orange",size=2) + 
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "orange",linetype = "dashed",size=.5)+
  geom_point(aes(y = ESA_LAI), color="steelblue",size=2)+
  geom_smooth(aes(x=LAI,y = ESA_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  geom_point(aes(y = ann2_LAI), color="darkred",size=2) +
  geom_smooth(aes(x=LAI,y = ann2_LAI),color = "darkred",linetype = "dashed",size=.5)+
  xlab("NEON LAI") + ylab("LAI (models)")+
  theme_bw()

p.ESA <- ggplot(final.csv.new, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = ESA_LAI), color="steelblue",size=2)+
  geom_smooth(aes(x=LAI,y = ESA_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  xlab("NEON LAI") + ylab("ESA LAI")+
  theme_bw()

p.SCE <- ggplot(final.csv.new, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "orange",size=2) + 
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "orange",linetype = "dashed",size=.5)+
  xlab("NEON LAI") + ylab("SCE LAI")+
  theme_bw()

p.me <- ggplot(final.csv.new, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = ann2_LAI), color="darkred",size=2) +
  geom_smooth(aes(x=LAI,y = ann2_LAI),color = "darkred",linetype = "dashed",size=.5)+
  xlab("NEON LAI") + ylab("ANN LAI")+
  theme_bw()

p.ESA.me <- ggplot(final.csv.new, aes(x=ESA_LAI))+
  geom_point(aes(y = ann2_LAI), color="darkred",size=2) +
  geom_smooth(aes(x=ESA_LAI,y = ann2_LAI),color = "darkred",linetype = "dashed",size=.5)+
  xlab("ESA LAI") + ylab("Multi-output LAI")+
  theme_bw()

p.ESA.SCE <- ggplot(final.csv.new, aes(x=ESA_LAI))+
  geom_point(aes(y = SCE_LAI), color = "orange",size=2) + 
  geom_smooth(aes(x=ESA_LAI,y = SCE_LAI),method = "lm",color = "orange",linetype = "dashed",size=.5)+
  xlab("ESA LAI") + ylab("SCE LAI")+
  theme_bw()

p.SCE.me <- ggplot(final.csv.new, aes(x=SCE_LAI))+
  geom_point(aes(y = ann2_LAI), color="darkred",size=2) +
  geom_smooth(aes(x=SCE_LAI,y = ann2_LAI),color = "darkred",linetype = "dashed",size=.5)+
  xlab("SCE LAI") + ylab("Multi-output LAI")+
  theme_bw()


p.all



library(gridExtra)
grid.arrange(p.ESA ,p.SCE,p.me,
             p.ESA.SCE,p.ESA.me,p.SCE.me,  
             nrow = 2,ncol=3,
             top="Vs NEON Data")



