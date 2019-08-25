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

#for scaling and unscaling variables
library(DMwR)

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
  15,45, #Cab
  10,50, #Car
  0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.01,0.02, #Cm
  0.1,4.5),#LAI
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
train.n <- 250 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
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

Y_mat.scale <- scale(Y.mat,center = T)
Y_mat.unscaled <- unscale(Y_mat.scale,Y_mat.scale)

head(Y_mat.unscaled)
head(Y.mat)

dim(Y.mat)
18+5
nn.5000 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,hidden.layers =c(10,6),loss.type = "squared",
                         activ.functions = c("relu",'tanh'),n.epochs = 5000,standardize = T,random.seed = 1000)

nn.pred.5000 <- predict(nn.5000,X.mat.valid)
head(nn.pred.5000$predictions)

nn.pred.5000$predictions <- unscale(nn.pred.5000$predictions,Y_mat.scale)
head(nn.pred.5000$predictions)

#which activation functions work better
library(gtools)
help(permn)
list.activ <- c('tanh',"sigmoid",'relu','linear','ramp')
list.activ <- c('tanh','relu','linear','ramp')
list.activ <- c('tanh','relu','linear')
act.func.list <- permutations(n=3,r=3,v=list.activ,repeats.allowed=T)
act.func.list
nrow(act.func.list)

for (i in 1:nrow(act.func.list)){
  ann.func <- c(act.func.list[i,1],act.func.list[i,2],act.func.list[i,3])
  print(ann.func)
  nn.temp <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,hidden.layers =c(10,9,5),loss.type = "squared",
                           activ.functions = ann.func,n.epochs = 1000,standardize = T,random.seed = 100)
}


#lets see how the plots look like
ann.temp.df <- data.frame(trait=names(valid.trait.df[,-c(6:8)]),
                          slope=NA,inter=NA,Rsqur=NA,RMSE=NA,MAE=NA,MAPE=NA)

rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }

### ploting the ANN
par(mfrow=c(2,3))
for (i in 1:nrow(ann.temp.df)){
  
  #linear model stats
  x <- valid.trait.df[,i]
  y <- nn.pred.5000$predictions[,i]
  
  plot(y,x,
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=ann.temp.df[i,1],pch=19,cex=.5,col="salmon")
  
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

#lets further train our network with the newly added data


#############################################################
############## Working on neon data stars here ##############
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

#now lets do some predictions of my ANN
nn.pred.5000 <- predict(nn.5000,neon.s2.obs.mat)

#lets wrap it up in a dataframe
test.df <- as.data.frame(nn.pred.5000)
names(neon.valid.df)
neon.valid.df.neonOnly <- neon.valid.df[,c(4:28)]
names(neon.valid.df.neonOnly)

test.df <- cbind(test.df,neon.valid.df.neonOnly)
names(test.df)



#The way i did my sampling (using a lwer res map) results in some outliers in the validation data
#lets remove outliers in each index and also reflectances that are offbeat

test.df.redux <- test.df
for (i in 6:ncol(test.df)){
  outlier.redux <- boxplot(test.df.redux[,i], plot=FALSE)$out
  print(length(outlier.redux))
  
  if (length(outlier.redux)!=0){
    test.df.redux<-test.df.redux[-which(test.df.redux[,i] %in% outlier.redux),]
  }
  boxplot(test.df.redux[,i],main=names(test.df.redux)[i])
}

dim(test.df.redux)
summary(test.df.redux)

#using the function of https://www.mdpi.com/2072-4292/11/13/1597/htm - AGB as a function of LAI*cm
par(mfrow=c(1,1))
plot(test.df.redux$LAI,test.df.redux$predictions.LAI,xlab="NEON data",ylab="Prediction (ANN): LAI")
#plot(test.df.redux$NDWI,test.df.redux$predictions.Cw,xlab="NEON data",ylab="Prediction (ANN): LAI")

#canopy water content indices
par(mfrow=c(3,2))

plot(test.df.redux$MSI ,test.df.redux$predictions.Cw,xlab="NEON data: MSI",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NDII,test.df.redux$predictions.Cw,xlab="NEON data: ndii",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NDWI,test.df.redux$predictions.Cw,xlab="NEON data: ndwi",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NDMI,test.df.redux$predictions.Cw,xlab="NEON data: ndmi",ylab="Prediction (mRF): Cw")
plot(test.df.redux$WBI ,test.df.redux$predictions.Cw,xlab="NEON data: wbi",ylab="Prediction (mRF): Cw")

##########################
########### recreating the indices

par(mfrow=c(1,1))
plot(test.df.redux$BIOM ,test.df.redux$predictions.LAI*test.df.redux$predictions.Cm,xlab="NEON data: biomass",ylab="Prediction: LAI*CM")

names(test.df.redux)
summary(test.df.redux)
new.param.list.redux <- test.df.redux[,1:5]
head(new.param.list.redux)
#names(new.param.list.redux) <- names(as.data.frame(out.mRF.test.pred ))
names(new.param.list.redux) <- c("Cab","Car","Cw","Cm","LAI")
head(new.param.list.redux)

#new.param.list.redux$TypeLidf="Erectophile"
new.param.list.redux$tts=sun_zenith #fixed zenith, azimuth and rel azimuth
new.param.list.redux$tto=obs_zenith
new.param.list.redux$psi=rel_azimut
#new.param.list.redux$psoil=1  #psoil - fixed this as the average that i saw from the paper [0.5 to 1.5]
#new.param.list.redux$lidfa=65
head(new.param.list.redux)
#names(new.param.list.redux) <- c("Cab","Car","Cw","Cm","LAI","tts","tto","psi")
summary(new.param.list.redux)
pred.prosail.scplib.redux <- PROSAIL(parameterList = new.param.list.redux)

par(mfrow=c(1,1))
plot(pred.prosail.scplib.redux)

#now we create a dataframe to receive everything
pred.spectra <- as.data.frame(spectra(pred.prosail.scplib))
names(pred.spectra) <- paste("band_",pred.prosail.scplib@wavelength,sep="")
names(pred.spectra)

#now we calculate the indices
pred.indices <- data.frame(linenr=seq(1:nrow(pred.spectra)))
#sources: https://www.harrisgeospatial.com/docs/CanopyWaterContent.html
pred.indices$MSI  <- pred.spectra$band_1599/pred.spectra$band_819
pred.indices$NDII <- (pred.spectra$band_819-pred.spectra$band_1649)/(pred.spectra$band_819+pred.spectra$band_1649)
pred.indices$NDWI <- (pred.spectra$band_857-pred.spectra$band_1241)/(pred.spectra$band_857+pred.spectra$band_1241)
pred.indices$NDMI <- (pred.spectra$band_860-(pred.spectra$band_1640-pred.spectra$band_2130))/(pred.spectra$band_860+(pred.spectra$band_1640-pred.spectra$band_2130))
pred.indices$WBI  <- pred.spectra$band_970/pred.spectra$band_900

