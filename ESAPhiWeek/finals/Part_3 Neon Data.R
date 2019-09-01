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
  0.5,4.5),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

#creating a training space
train.n <- 3000
#train.grid <- Grid(param.maxmin,train.n)
#train.LHS <- Latinhyper(param.maxmin,train.n-nrow(train.grid))
train.LHS <- Latinhyper(param.maxmin,train.n)

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


#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c("tanh","relu")
lrates <- 0.005
epochs <- 3000
optype <- "adam"

ann.3000 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype)
#with grid  0.172064
#without grid < 0.1

#lets check how our model looks likes
ann.3000.pred <- predict(ann.3000,X.mat.valid)
ann.3000.pred.df <- as.data.frame(ann.3000.pred)
names(ann.3000.pred.df) <- paste("ann",names(valid.trait.df[,-c(6:8)]),sep="_")

#visual model check up
model.df <- cbind(valid.trait.df[,-c(6:8)],ann.3000.pred.df)
plot(model.df$Cab,model.df$ann_Cab)
plot(model.df$Car,model.df$ann_Car)
plot(model.df$Cw,model.df$ann_Cw)
plot(model.df$Cm,model.df$ann_Cm)
plot(model.df$LAI,model.df$ann_LAI)


#############################################################
############## Working on neon data stars here ##############
############### But for now only for ML #####################
#############################################################

#now, we will get the samples from the point shapefile we created elsewhere - this case only grasslands
#neon.test.shp <- readShapePoints("D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_small_pts_5kSample.shp")
neon.test.shp <- readShapePoints("D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_pts_SmallSample.shp")

#now lets extract the values from the sentinel image
s2.AOI.stack <- stack("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif")/10000 #to convert to reflectance

#the names of the layers come in different order, for simplicity lets make it all the same
names(s2.AOI.stack)<-c("B02","B03","B04",
                       "B05","B06","B07",
                       "B8A","B11","B12")

s2.AOI.refl.shp <- raster::extract(s2.AOI.stack,
                                   neon.test.shp,
                                   sp=T)

#now lets load all the data from the traits of NEON and have that also in our shape
#loading the raster
path2neon <- list.files("D:/NEON_Data/CLBJ/Processed/",
                        pattern=".tif",
                        full.names=T)
name2neon <- list.files("D:/NEON_Data/CLBJ/Processed/",
                        pattern=".tif")

#now we remove the one that doesn't matter.. BEWARE - if you change anything .tif in th folder, this might change
name2neon
path2neon <- path2neon[-1]
name2neon <- name2neon[-1]
path2neon
name2neon

#they also need to be cropped to the same extents or we cant stack them
AOI_small.shp <- readShapePoly("D:/NEON_Data/CLBJ/qgis_shape_bound/qgis_shape_bound_small_Area_diss.shp")
for (i in 1:length(path2neon)){
  
  out.file <- paste("D:/NEON_Data/CLBJ/Processed_cropped/Small_",name2neon[i],sep="")
  print(out.file)
  crop(x=raster(path2neon[i]), y=AOI_small.shp, 
       filename=out.file,
       options=c("COMPRESS=LZW"),
       overwrite=TRUE)
}

#loading the raster
path2neon.small <- list.files("D:/NEON_Data/CLBJ/Processed_cropped/",
                              pattern=".tif",
                              full.names=T)
name2neon.small <- list.files("D:/NEON_Data/CLBJ/Processed_cropped/",
                              pattern=".tif")

#now we can stack and extract everything in one go
valid.shp <- raster::extract(stack(path2neon.small ),
                             s2.AOI.refl.shp,
                             sp=T)

dim(valid.shp)
head(valid.shp)
names(valid.shp)
summary(valid.shp)
#remoing the sel tag 
valid.shp <- valid.shp[,-4]
head(valid.shp)
#now we want proper names... these names are so freaking bad
names(valid.shp) <- c(names(valid.shp)[1:12],
                      c("ARVI","BIOM","CHM",
                        "EVI","fPAR","LAI",
                        "MSI","NDII","NDLI",
                        "NDNI","NDVI","NDWI",
                        "NDMI","PRI","SAVI",
                        "WBI"))

#no we have everything to make some prediction, we can go use our trained mRF to predict some outputs
neon.valid.df <- as.data.frame(valid.shp)
head(neon.valid.df)


#only the bands to avoid confusion
neon.s2.Obs.df.bandsOnly <- neon.valid.df[,c(4:12)]
head(neon.s2.Obs.df.bandsOnly)

neon.s2.obs.mat <- as.matrix(neon.s2.Obs.df.bandsOnly) #notice, these are just 


###################################################
############ Predicting the ANN on NEON data ######
###################################################


neon.pred.3000 <- predict(ann.3000,neon.s2.obs.mat)

#lets wrap it up in a dataframe
test.df <- as.data.frame(neon.pred.3000)
test.df.xy <- as.data.frame(neon.pred.3000)
names(test.df) <- paste("ann_neon",names(valid.trait.df[,-c(6:8)]),sep="_")
names(test.df.xy) <- paste("ann_neon",names(valid.trait.df[,-c(6:8)]),sep="_")


names(neon.valid.df)
neon.valid.df.neonOnly <- neon.valid.df[,c(4:28)]
neon.valid.df.neonOnly.xy <- neon.valid.df[,c(4:30)]
names(neon.valid.df.neonOnly.xy)

test.df <- cbind(test.df,neon.valid.df.neonOnly)
test.df.xy <- cbind(test.df.xy,neon.valid.df.neonOnly.xy)
names(test.df)
names(test.df.xy)


#The way i did my sampling (using a lwer res map) results in some outliers in the validation data
#lets remove outliers in each index and also reflectances that are offbeat
test.df.redux <- test.df
neon.valid.df.redux <- test.df.xy 

for (i in 6:ncol(test.df)){
  outlier.redux <- boxplot(test.df.redux[,i], plot=FALSE)$out
  print(length(outlier.redux))
  
  if (length(outlier.redux)!=0){
    test.df.redux<-test.df.redux[-which(test.df.redux[,i] %in% outlier.redux),]
    #neon.valid.df.redux <- neon.valid.df[-which(test.df.redux[,i] %in% outlier.redux),]
    neon.valid.df.redux <- neon.valid.df.redux[-which(neon.valid.df.redux[,i] %in% outlier.redux),]
  }
  #boxplot(test.df.redux[,i],main=names(test.df.redux)[i])
}

plot(test.df.redux$ann_neon_Car,neon.valid.df.redux$ann_neon_Car)

dim(test.df.redux)
#lets do a quick check
par(mfrow=c(1,1))
plot(test.df.redux$LAI,test.df.redux$ann_neon_LAI,xlab="NEON data",ylab="Prediction (ANN): LAI")
abline(lm(test.df.redux$ann_neon_LAI~test.df.redux$LAI))
#plot(test.df.redux$NDWI,test.df.redux$predictions.Cw,xlab="NEON data",ylab="Prediction (ANN): LAI")



#time for optimization

###############################################################
########### SCE Optimization from herein ######################
###############################################################

#first the cost function

sam.fun.spclib.5 <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  Cm_init <- list.params[4]
  Car_init <- list.params[5]
  
  #sentinel position
  sun_zenith = 25.8734
  obs_zenith = 5.4748
  rel_azimut = 281.0692 - 5.4748
  
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
bands.optim <- test.df.redux[,6:14]
wavelength(s2.valid.spclib.20m)

#as an example
target.spectra <- speclib(spectra = as.matrix(bands.optim[1,]),wavelength(s2.valid.spclib.20m))


param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  1,50, #Cab
  1,50, #Car
  0.01,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.01,0.03, #Cm
  0.5,4.5),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)


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
names(test.df.redux)
out.df.5 <- data.frame(LAI=test.df.redux[,c(20)])
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

write.csv(out.df.5 ,paste(dump.neon.results,"Optim_Iterations.csv",sep="/"))


#######################################################
head(neon.valid.df.redux)
plot(neon.valid.df.redux$LAI,out.df.5$SCE_LAI)
plot(out.df.5$LAI,out.df.5$SCE_LAI)

head(out.df.5)
#lets append the coordinates to a data frame
head(test.df.redux)

temp.final.df <- cbind(out.df.5[,-c(1,7)],neon.valid.df.redux)
head(temp.final.df)

write.csv(temp.final.df,paste(dump.neon.results,"NEON_SCE_and_ANN_xy.csv",sep="/"))
plot(temp.final.df$SCE_LAI,temp.final.df$ann_neon_LAI)

#lets make it in coordinates
temp.final.shp <- temp.final.df

coordinates(temp.final.shp) <- ~coords.x1+coords.x2
library(sp)
library(maptools)
writePointsShape(temp.final.shp,paste(dump.neon.results,"NEON_SCE_and_ANN_xy.shp",sep="/"))


#extracting values from ESA Raster
ESA.bio <- stack("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI_ESABio.tif")

par(mfrow=c(1,1))
plot(ESA.bio)
names(ESA.bio)<- c("ESA_LAI","LAI_flag",
                   "ESA_LAI_cab","ESA_LAI_cab_flag",
                   "ESA_LAI_car","ESA_LAI_car_flag",
                   "ESA_fapar","ESA_fapar_flag",
                   "ESA_fcover","ESA_fcover_flag")

final.shp <- extract(ESA.bio,temp.final.shp,sp=T)

writePointsShape(final.shp,paste(dump.neon.results,"NEON_SCE_and_ANN_and_ESA_xy.shp",sep="/"))
final.shp.df <- as.data.frame(final.shp)
write.csv(final.shp.df ,paste(dump.neon.results,"NEON_SCE_and_ANN_and_ESA_xy.csv",sep="/"))

plot(final.shp.df$ann_neon_LAI,final.shp$ESA_LAI)

#RESULTS ARE NOT GOOD - but i can retry the ANN - perhaps with different activations