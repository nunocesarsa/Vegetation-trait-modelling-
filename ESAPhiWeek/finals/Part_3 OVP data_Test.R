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

#setwd
setwd("D:/ESAPhiWeek/")
gc()
dump.esa.results <- "./Out_ESA"
dump.neon.results <- "./Out_Neon"
dump.ml.opti.fld <- "./Out_OptimVsML"
dump.ml.fld <- "./Out_MachL"
dump.fld <- "./Out_Optim"

dir.create(dump.esa.results)


#useful functions for later
rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }

gc()

#first lets extract the data from OVP 
grass.shp <- readShapePoints("D:/ESAPhiWeek/ESA_data/OVP_GrassSamples_PHIWEEK.shp")
S2.20m.img <- stack("D:/ESAPhiWeek/ESA_data/S2B_MSIL2A_20190717_T31UFU_Small.tif")/10000

plot(S2.20m.img$S2B_MSIL2A_20190717_T31UFU_Small.9)
points(grass.shp)

names(S2.20m.img)<-c("B02","B03","B04",
                     "B05","B06","B07",
                     "B8A","B11","B12")
#loading bio
ESA.bio.20m <- stack("D:/ESAPhiWeek/ESA_Bio_Snap/S2B_MSIL2A_20190717T104029_resampled_biophysical.tif")
ESA.bio.20m.LAI <- raster(ESA.bio.20m,1)
names(ESA.bio.20m.LAI)<-c("ESA_LAI")

s2.OVP.shp <- raster::extract(S2.20m.img,
                              grass.shp,
                                   sp=T)

s2.OVP.shp <- raster::extract(ESA.bio.20m.LAI,
                              s2.OVP.shp,
                              sp=T)

s2.OVP.df <- as.data.frame(s2.OVP.shp)
head(s2.OVP.df)

#lets see a summary of ESA LAI
summary(s2.OVP.df)


##################################################
############ herein ANN training #################
#################################################

#lets train the ANN
param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  1,50, #Cab
  1,50, #Car
  0.005,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.005,0.02, #Cm
  0.2,8.5),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)


#creating a training space
train.n <- 10000
#train.grid <- Grid(param.maxmin,train.n)
#train.LHS <- Latinhyper(param.maxmin,train.n-nrow(train.grid))
train.LHS <- Latinhyper(param.maxmin,train.n)

#to force having the limits we will add something oo the train
#train.LHS <- rbind(train.LHS,train.grid)

#valid.n <- 2* nrow(train.LHS) #this represents the number of runs that prosail will be, x pts per Trait
valid.n <- 300
valid.LHS <- Latinhyper(param.maxmin,valid.n)


#checking the stats
summary(train.LHS)
summary(valid.LHS)

#sentinel position very randomly
#sun zen and rel azi are stored in the metadata but the obs zenith is per band 
sun_zenith = 33.4788548875136
obs_zenith = 9.5
rel_azimut = 154.782326698598 - 107 #this is a total guess from the histogram.. esa is so stupid lol



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


#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c("tanh","relu")
lrates <- 0.005
epochs <- 3000
optype <- "adam"

set.seed(2000)
ann.3000 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype)
head(Y.mat)
Y.mat.LAI <- Y.mat[,5]
#ann.3000.LAI <- neuralnetwork(X=X.mat,y=Y.mat.LAI,regression=T,
#                          hidden.layers =net.st,
#                          loss.type = "squared",
#                          activ.functions = act.fn,
#                          n.epochs = epochs,
#                          standardize = T,
#                          learn.rates = lrates,
#                          optim.type = optype)



#lets check how our model looks likes
ann.3000.pred <- predict(ann.3000,X.mat.valid)
ann.3000.pred.df <- as.data.frame(ann.3000.pred)
names(ann.3000.pred.df) <- paste("ann",names(valid.trait.df[,-c(6:8)]),sep="_")

ann.3000.pred.LAI <- predict(ann.3000.LAI,X.mat.valid)
ann.3000.pred.df.LAI <- as.data.frame(ann.3000.pred.LAI)
#names(ann.3000.pred.df.LAI) <- paste("ann",names(valid.trait.df[,-c(6:8)]),sep="_")


#visual model check up
model.df <- cbind(valid.trait.df[,-c(6:8)],ann.3000.pred.df)
plot(model.df$Cab,model.df$ann_Cab)
plot(model.df$Car,model.df$ann_Car)
plot(model.df$Cw,model.df$ann_Cw)
plot(model.df$Cm,model.df$ann_Cm)
plot(model.df$LAI,model.df$ann_LAI)


#we can predict now
head(s2.OVP.df)
s2.OVP.df.bands <- s2.OVP.df[,c(2:10)]
head(s2.OVP.df.bands)
s2.OVP.df.bands.mat <- as.matrix(s2.OVP.df.bands)

ann.s2.pred <- predict(ann.3000,s2.OVP.df.bands.mat)
ann.s2.pred.df <- as.data.frame(ann.s2.pred )

ann.s2.pred.LAI <- predict(ann.3000.LAI,s2.OVP.df.bands.mat)
ann.s2.pred.df.LAI <- as.data.frame(ann.s2.pred.LAI )

#quick check
head(ann.s2.pred.df$predictions.LAI)
par(mfrow=c(1,1))
plot(s2.OVP.df$ESA_LAI,ann.s2.pred.df$predictions.LAI)
plot(s2.OVP.df$ESA_LAI,ann.s2.pred.df.LAI$y_1)

temp.df <- ann.s2.pred.df
temp.df$ESALai <- s2.OVP.df$ESA_LAI

write.csv(temp.df,paste(dump.esa.results,"ESA_Preds_10k.csv",sep="/"))

################################################
########### herein optimization ################
################################################


sam.fun.spclib.5 <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  Cm_init <- list.params[4]
  Car_init <- list.params[5]
  
  #sun zen and rel azi are stored in the metadata but the obs zenith is per band 
  sun_zenith = 33.4788548875136
  obs_zenith = 9.5
  rel_azimut = 154.782326698598 - 107
  
  init.param <- data.frame(Cab=Cab_init,
                           Cw=Cw_init,
                           LAI=LAI_init,
                           Cm=Cm_init,
                           Car=Car_init,
                           tts = sun_zenith,
                           tto = obs_zenith,
                           psi = rel_azimut)
  
  #this fetches new prosail parameters.. if you want to minimize towards an observation, then this, should be the observation
  
  if (is.null(target.spectra)==T){
    print("Warning: its minimizing to the default PROSAIL parameters") #uncomment if you want to get a spammy console
    
    target.spectra <- spectralResampling(target.spectra,
                                         "Sentinel2",response_function = TRUE)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    
    my.spectra <- spectralResampling(PROSAIL(parameterList = init.param),
                                     "Sentinel2",response_function = TRUE)
    
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,def.prosail)
    fin.val <- spec.dist[1]
  }
  
  if (is.null(target.spectra)==F) {
    
    #print("Attempting to find the vegetation traits") #uncomment if you want to get a spammy console
    #print(prosail.param)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    my.spectra <- spectralResampling(PROSAIL(parameterList = init.param),
                                     "Sentinel2",response_function = TRUE)
    my.spectra <- my.spectra[,c(2,3,4,5,6,7,9,12,13)]
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,target.spectra)
    fin.val <- spec.dist[1]
    
  }
  
  
  return(fin.val)
}


#now we must create a spectral lib from our reflectance bands
bands.optim <- s2.OVP.df.bands 
wavelength(s2.valid.spclib.20m)

#as an example
target.spectra <- speclib(spectra = as.matrix(bands.optim[1,]),wavelength(s2.valid.spclib.20m))
plot(target.spectra)


#the limits of our convergence - they are are fitting the TRAIT order inside the COST FUNCTION
llim.5 <- c(param.maxmin[1,1],
            param.maxmin[3,1],
            param.maxmin[5,1],
            param.maxmin[4,1],
            param.maxmin[2,1])

ulim.5 <- c(param.maxmin[1,2],
            param.maxmin[3,2],
            param.maxmin[5,2],
            param.maxmin[4,2],
            param.maxmin[2,2])

#now a random set of initial conditions
min.max.Cab <- param.maxmin[1,]
min.max.Car <- param.maxmin[2,]
min.max.Cw <- param.maxmin[3,]
min.max.Cm <- param.maxmin[4,]
min.max.LAI <- param.maxmin[5,]
init.par.df <- data.frame(Cab=rtnorm(1,mean=mean(min.max.Cab),sd=mean(min.max.Cab)/4,lower=min.max.Cab[1],upper=min.max.Cab[2]), 
                          Cw= rtnorm(1,mean=mean(min.max.Cw),sd=mean(min.max.Cw)/4,lower=min.max.Cw[1],upper=min.max.Cw[2]),
                          LAI=rtnorm(1,mean=mean(min.max.LAI),sd=mean(min.max.LAI )/4,lower=min.max.LAI[1],upper=min.max.LAI[2]),
                          Cm= rtnorm(1,mean=mean(min.max.Cm),sd=mean(min.max.Cm)/4,lower=min.max.Cm[1],upper=min.max.Cm[2]),
                          Car=rtnorm(1,mean=mean(min.max.Car),sd=mean(min.max.Car)/4,lower=min.max.Car[1],upper=min.max.Car[2])) 
init.par.5 <- as.numeric(init.par.df[1,])

#a place to hold the ouptus


out.df.5 <- data.frame(ESA_LAI=s2.OVP.df[,c(11)])
out.df.5$SCE_Cab <- NA
out.df.5$SCE_Cw  <- NA
out.df.5$SCE_LAI <- NA
out.df.5$SCE_Cm  <- NA
out.df.5$SCE_Car <- NA
out.df.5$SCE_iter <- NA

#sometimes this loop breaks down, seems to me the optimization calls for parameters completely outside PROSAIL range, generates NaN and thus fails SAM function
#for( i in 14:17){
#i <- 1
############### THIS TAKES LIKE 6H #####################
for( i in 1:nrow(out.df.5)){
  print(paste("Estimating sample:", i))
  
  #print(df.param.list[i,])
  #tgt.spec <- PROSAIL(parameterList = df.param.list[i,])
  tgt.spec  <- speclib(spectra = as.matrix(bands.optim[i,]),wavelength(s2.valid.spclib.20m))
  
  SCEoptim.pred <- SCEoptim(FUN = sam.fun.spclib.5,
                            par = init.par.5,
                            #method="L-BFGS-B",
                            lower=llim.5,
                            upper=ulim.5,
                            target.spectra=tgt.spec)
  
  #storing the outputs
  out.df.5$SCE_Cab[i] <- SCEoptim.pred$par[1]
  out.df.5$SCE_Cw[i]  <- SCEoptim.pred$par[2]
  out.df.5$SCE_LAI[i] <- SCEoptim.pred$par[3]
  out.df.5$SCE_Cm[i]  <- SCEoptim.pred$par[4]
  out.df.5$SCE_Car[i] <- SCEoptim.pred$par[5]
  
  
  out.df.5$SCE_iter[i] <- SCEoptim.pred$iterations
  
}

SCEoptim.pred$par
init.par.5
llim.5
ulim.5

plot(out.df.5$ESA_LAI,out.df.5$SCE_LAI)
plot(out.df.5$ESA_LAI,ann.s2.pred.df$predictions.LAI,pch=14)
points(out.df.5$ESA_LAI,out.df.5$SCE_LAI,pch=19)

final.df <- cbind(out.df.5,ann.s2.pred.df)
head(final.df)
summary(final.df)

write.csv(final.df,paste(dump.esa.results,"ESA_Preds_new.csv",sep="/"))
