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
dump.ml.fld <- "./Out_MachL"
dir.create(dump.ml.fld) #it gives out a warning if the folder exits
#setting seed
set.seed(2000)


#first we generate a dataset which we will divide in 3 parts:
#now we set up the upper and lower limits (important! to ensure convergence)
param.maxmin <- matrix(c(5,100, #Cab
                         5,50, #Car
                         0.005,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.03, #Cm
                         0.5,9.5),#LAI
                       #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)

min.max.Cab <- param.maxmin[1,]
min.max.Car <- param.maxmin[2,]
min.max.Cw <- param.maxmin[3,]
min.max.Cm <- param.maxmin[4,]
min.max.LAI <- param.maxmin[5,]

#now we create a set of initial parameters to be able to start prosail
LHS.0500 <- Latinhyper(param.maxmin,0500)
LHS.1500 <- Latinhyper(param.maxmin,1500)
LHS.3000 <- Latinhyper(param.maxmin,3000)

LHS.valid <- Latinhyper(param.maxmin,6000)

#using the latin hypercube approach for sample generation
train.par.df.0500 <- data.frame(Cab=LHS.0500[,1],
                                Car=LHS.0500[,2],
                                Cw= LHS.0500[,3],
                                Cm= LHS.0500[,4],
                                LAI=LHS.0500[,5])

train.par.df.1500 <- data.frame(Cab=LHS.1500[,1],
                                Car=LHS.1500[,2],
                                Cw= LHS.1500[,3],
                                Cm= LHS.1500[,4],
                                LAI=LHS.1500[,5])

train.par.df.3000 <- data.frame(Cab=LHS.3000[,1],
                                Car=LHS.3000[,2],
                                Cw= LHS.3000[,3],
                                Cm= LHS.3000[,4],
                                LAI=LHS.3000[,5])

valid.par.df <- data.frame(Cab=LHS.valid[,1],
                           Car=LHS.valid[,2],
                           Cw= LHS.valid[,3],
                           Cm= LHS.valid[,4],
                           LAI=LHS.valid[,5])
                           

#useful for later
param.names <- names(train.par.df.0500)

############## converting to df
#we convert it to a spectral object
spclib.plist.0500 <-PROSAIL(parameterList = train.par.df.0500 )
spclib.plist.1500 <-PROSAIL(parameterList = train.par.df.1500 )
spclib.plist.3000 <-PROSAIL(parameterList = train.par.df.3000 )
spclib.plist.valid <- PROSAIL(parameterList = valid.par.df)

#we will focus on only sentinel 2 data
s2.spclib.0500 <- spectralResampling(spclib.plist.0500,"Sentinel2",response_function = T)
s2.spclib.1500 <- spectralResampling(spclib.plist.1500,"Sentinel2",response_function = T)
s2.spclib.3000 <- spectralResampling(spclib.plist.3000,"Sentinel2",response_function = T)
s2.spclib.valid <- spectralResampling(spclib.plist.valid,"Sentinel2",response_function = T)

wavelength(s2.spclib.0500)

#only the 20m bands
s2.spclib.0500.20m <- s2.spclib.0500[,c(2,3,4,5,6,7,9,12,13)]
s2.spclib.1500.20m <- s2.spclib.1500[,c(2,3,4,5,6,7,9,12,13)] 
s2.spclib.3000.20m <- s2.spclib.3000[,c(2,3,4,5,6,7,9,12,13)]
s2.spclib.valid.20m <- s2.spclib.valid[,c(2,3,4,5,6,7,9,12,13)]

s2.spectra.0500 <- as.data.frame(spectra(s2.spclib.0500.20m))
s2.spectra.1500 <- as.data.frame(spectra(s2.spclib.1500.20m))
s2.spectra.3000 <- as.data.frame(spectra(s2.spclib.3000.20m))
s2.spectra.valid <- as.data.frame(spectra(s2.spclib.valid.20m))

s2.band.names <- c("B02","B03","B04",
                   "B05","B06","B07",
                   "B8A","B11","B12")

names(s2.spectra.0500) <- s2.band.names
names(s2.spectra.1500) <- s2.band.names
names(s2.spectra.3000) <- s2.band.names
names(s2.spectra.valid) <- s2.band.names

#lets combine all the data
df.0500 <-  cbind(train.par.df.0500,s2.spectra.0500)
df.1500 <-  cbind(train.par.df.1500,s2.spectra.1500)
df.3000 <-  cbind(train.par.df.3000,s2.spectra.3000)
df.valid <- cbind(valid.par.df,s2.spectra.valid)


#lets train the ANN
Y.mat.0500 <- as.matrix(df.0500[,c(1:5)])
Y.mat.1500 <- as.matrix(df.1500[,c(1:5)])
Y.mat.3000 <- as.matrix(df.3000[,c(1:5)])

X.mat.0500 <- as.matrix(df.0500[,c(6:ncol(df.0500))])
X.mat.1500 <- as.matrix(df.1500[,c(6:ncol(df.1500))])
X.mat.3000 <- as.matrix(df.3000[,c(6:ncol(df.1500))])


rsq <- function (x, y) cor(x, y) ^ 2


#first we do the kfolding and select the best models - each set must be k folded independently
nrfolds <- 5 #a number of folds to use
fold.0500 <- kfold(df.0500,nrfolds)
fold.1500 <- kfold(df.1500,nrfolds)
fold.3000 <- kfold(df.3000,nrfolds)
table(fold.0500)
table(fold.1500)
table(fold.3000)

#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c('relu',"sigmoid")
lrates <- 0.005
epochs <- 2000
optype <- "adam"


param.names

#loop for 500 input points
#i <- 1
train.df <- df.0500
fold.selection <- fold.0500
for (i in 1:nrfolds){
  
  temp.train.df <- train.df[which(fold.selection != i),]
  temp.valid.df <- train.df[which(fold.selection == i),]
  
  
  
  #creating the mRF model
  temp.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                    data=temp.train.df,block.size = 25)
  temp.pred.mRF <- predict(temp.mRF,
                           newdata=temp.valid.df)
  temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))
  
  #creating the ANN model
  temp.X.mat <- as.matrix(temp.train.df[,c(6:ncol(temp.train.df))])
  temp.Y.mat <- as.matrix(temp.train.df[,c(1:5)])
  
  
  temp.X.mat.valid <- as.matrix(temp.valid.df[,c(6:ncol(temp.valid.df))])
  
  
  temp.ann <- neuralnetwork(X=temp.X.mat,y=temp.Y.mat,regression=T,
                           hidden.layers =net.st,
                           loss.type = "squared",
                           activ.functions = act.fn,
                           n.epochs = epochs,
                           standardize = T,
                           learn.rates = lrates,
                           optim.type = optype)
  
  #now i could iterate for each trait but im to lazy for that
  temp.ann.pred <- predict(temp.ann,temp.X.mat.valid)
  
  head(temp.mv.pred.mRF)
  head(temp.ann.pred)
  
  
  for (j in 1:length(param.names)){
    
    print(paste("Storing statistics for ",
                i,"ith fold of variable",param.names[j]))
    
    
    if ((i == 1) & (j == 1)){
      
      #the creates a data frame to store the outputs
      #j<- 1
      k.fold.mRF.df <- data.frame(kf_th=i,Model="mRF",variable=param.names[j])
      k.fold.ann.df <- data.frame(kf_th=i,Model="ann",variable=param.names[j])
      
      #this fetches the variables
      mRF.pred <- temp.mv.pred.mRF[,j]
      ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      lm.mRF <- lm(mRF.pred ~ temp.valid.df[,j])
      lm.ann <- lm(ann.pred ~ temp.valid.df[,j])
      
      #stores the data
      k.fold.mRF.df$Slpe <- lm.mRF$coefficients[[2]]
      k.fold.mRF.df$Inte <- lm.mRF$coefficients[[1]]
      k.fold.mRF.df$Rsqr <- summary(lm.mRF)$r.square
      k.fold.mRF.df$MAE  <- DescTools::MAE(lm.mRF)
      k.fold.mRF.df$MAPE <- DescTools::MAPE(lm.mRF)
      k.fold.mRF.df$SMAPE <- DescTools::SMAPE(lm.mRF)
      k.fold.mRF.df$MSE <- DescTools::MSE(lm.mRF)
      k.fold.mRF.df$RMSE <- DescTools::RMSE(lm.mRF)
      
      k.fold.ann.df$Slpe <- lm.ann$coefficients[[2]]
      k.fold.ann.df$Inte <- lm.ann$coefficients[[1]]
      k.fold.ann.df$Rsqr <- summary(lm.ann)$r.square
      k.fold.ann.df$MAE  <- DescTools::MAE(lm.ann)
      k.fold.ann.df$MAPE <- DescTools::MAPE(lm.ann)
      k.fold.ann.df$SMAPE <- DescTools::SMAPE(lm.ann)
      k.fold.ann.df$MSE <- DescTools::MSE(lm.ann)
      k.fold.ann.df$RMSE <- DescTools::RMSE(lm.ann)
      
      
      #plot checking
      #plot(temp.valid.df[,i],temp.valid.df[,i])
      #points(temp.valid.df[,i],mRF.pred,pch=19,col="blue")
      #points(temp.valid.df[,i],ann.pred,pch=19,col="red")
      
      
      
    } else {
      
      #print(j)
      #we have to repeat above and append the dataframe
      
      loop.k.fold.mRF.df <- data.frame(kf_th=i, Model="mRF", variable=param.names[j])
      loop.k.fold.ann.df <- data.frame(kf_th=i, Model="ann", variable=param.names[j])
                          
      #this fetches the variables
      loop.mRF.pred <- temp.mv.pred.mRF[,j]
      loop.ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      loop.lm.mRF <- lm(loop.mRF.pred ~ temp.valid.df[,j])
      loop.lm.ann <- lm(loop.ann.pred ~ temp.valid.df[,j])
      
      #stores the data
      loop.k.fold.mRF.df$Slpe <- loop.lm.mRF$coefficients[[2]]
      loop.k.fold.mRF.df$Inte <- loop.lm.mRF$coefficients[[1]]
      loop.k.fold.mRF.df$Rsqr <- summary(loop.lm.mRF)$r.square
      loop.k.fold.mRF.df$MAE  <- DescTools::MAE(loop.lm.mRF)
      loop.k.fold.mRF.df$MAPE <- DescTools::MAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$SMAPE <- DescTools::SMAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$MSE <- DescTools::MSE(loop.lm.mRF)
      loop.k.fold.mRF.df$RMSE <- DescTools::RMSE(loop.lm.mRF)
      
      loop.k.fold.ann.df$Slpe <- loop.lm.ann$coefficients[[2]]
      loop.k.fold.ann.df$Inte <- loop.lm.ann$coefficients[[1]]
      loop.k.fold.ann.df$Rsqr <- summary(loop.lm.ann)$r.square
      loop.k.fold.ann.df$MAE  <- DescTools::MAE(loop.lm.ann)
      loop.k.fold.ann.df$MAPE <- DescTools::MAPE(loop.lm.ann)
      loop.k.fold.ann.df$SMAPE <- DescTools::SMAPE(loop.lm.ann)
      loop.k.fold.ann.df$MSE <- DescTools::MSE(loop.lm.ann)
      loop.k.fold.ann.df$RMSE <- DescTools::RMSE(loop.lm.ann)
      
      k.fold.mRF.df <- rbind(k.fold.mRF.df,loop.k.fold.mRF.df  )
      k.fold.ann.df <- rbind(k.fold.ann.df,loop.k.fold.ann.df  )
      
      
      
  
      
    } #end if
  } #end for params
} #end for k folds

write.csv(k.fold.mRF.df,paste(dump.ml.fld,
                              "S2_mRF_Kfolds_0500_samples.csv",sep = "/"))
write.csv(k.fold.ann.df,paste(dump.ml.fld,
                              "S2_ann_Kfolds_0500_samples.csv",sep = "/"))

train.df <- df.1500
fold.selection <- fold.1500
for (i in 1:nrfolds){
  
  temp.train.df <- train.df[which(fold.selection != i),]
  temp.valid.df <- train.df[which(fold.selection == i),]
  
  
  
  #creating the mRF model
  temp.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                    data=temp.train.df,block.size = 25)
  temp.pred.mRF <- predict(temp.mRF,
                           newdata=temp.valid.df)
  temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))
  
  #creating the ANN model
  temp.X.mat <- as.matrix(temp.train.df[,c(6:ncol(temp.train.df))])
  temp.Y.mat <- as.matrix(temp.train.df[,c(1:5)])
  temp.X.mat.valid <- as.matrix(temp.valid.df[,c(6:ncol(temp.valid.df))])
  
  
  temp.ann <- neuralnetwork(X=temp.X.mat,y=temp.Y.mat,regression=T,
                            hidden.layers =net.st,
                            loss.type = "squared",
                            activ.functions = act.fn,
                            n.epochs = epochs,
                            standardize = T,
                            learn.rates = lrates,
                            optim.type = optype)
  
  #now i could iterate for each trait but im to lazy for that
  temp.ann.pred <- predict(temp.ann,temp.X.mat.valid)
  
  #head(temp.mv.pred.mRF)
  #head(temp.ann.pred)
  
  
  for (j in 1:length(param.names)){
    
    print(paste("Storing statistics for ",
                i,"ith fold of variable",param.names[j]))
    
    
    if ((i == 1) & (j == 1)){
      
      #the creates a data frame to store the outputs
      #j<- 1
      k.fold.mRF.df <- data.frame(kf_th=i,Model="mRF",variable=param.names[j])
      k.fold.ann.df <- data.frame(kf_th=i,Model="ann",variable=param.names[j])
      
      #this fetches the variables
      mRF.pred <- temp.mv.pred.mRF[,j]
      ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      lm.mRF <- lm(mRF.pred ~ temp.valid.df[,j])
      lm.ann <- lm(ann.pred ~ temp.valid.df[,j])
      
      
      temp.mv.pred.mRF[,j]
      
      #stores the data
      k.fold.mRF.df$Slpe <- lm.mRF$coefficients[[2]]
      k.fold.mRF.df$Inte <- lm.mRF$coefficients[[1]]
      k.fold.mRF.df$Rsqr <- summary(lm.mRF)$r.square
      k.fold.mRF.df$MAE  <- DescTools::MAE(lm.mRF)
      k.fold.mRF.df$MAPE <- DescTools::MAPE(lm.mRF)
      k.fold.mRF.df$SMAPE <- DescTools::SMAPE(lm.mRF)
      k.fold.mRF.df$MSE <- DescTools::MSE(lm.mRF)
      k.fold.mRF.df$RMSE <- DescTools::RMSE(lm.mRF)
      
      k.fold.ann.df$Slpe <- lm.ann$coefficients[[2]]
      k.fold.ann.df$Inte <- lm.ann$coefficients[[1]]
      k.fold.ann.df$Rsqr <- summary(lm.ann)$r.square
      k.fold.ann.df$MAE  <- DescTools::MAE(lm.ann)
      k.fold.ann.df$MAPE <- DescTools::MAPE(lm.ann)
      k.fold.ann.df$SMAPE <- DescTools::SMAPE(lm.ann)
      k.fold.ann.df$MSE <- DescTools::MSE(lm.ann)
      k.fold.ann.df$RMSE <- DescTools::RMSE(lm.ann)
      
      
      #plot checking
      #plot(temp.valid.df[,i],temp.valid.df[,i])
      #points(temp.valid.df[,i],mRF.pred,pch=19,col="blue")
      #points(temp.valid.df[,i],ann.pred,pch=19,col="red")
      
      
      
    } else {
      
      #print(j)
      #we have to repeat above and append the dataframe
      
      loop.k.fold.mRF.df <- data.frame(kf_th=i, Model="mRF", variable=param.names[j])
      loop.k.fold.ann.df <- data.frame(kf_th=i, Model="ann", variable=param.names[j])
      
      #this fetches the variables
      loop.mRF.pred <- temp.mv.pred.mRF[,j]
      loop.ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      loop.lm.mRF <- lm(loop.mRF.pred ~ temp.valid.df[,j])
      loop.lm.ann <- lm(loop.ann.pred ~ temp.valid.df[,j])
      
      #stores the data
      loop.k.fold.mRF.df$Slpe <- loop.lm.mRF$coefficients[[2]]
      loop.k.fold.mRF.df$Inte <- loop.lm.mRF$coefficients[[1]]
      loop.k.fold.mRF.df$Rsqr <- summary(loop.lm.mRF)$r.square
      loop.k.fold.mRF.df$MAE  <- DescTools::MAE(loop.lm.mRF)
      loop.k.fold.mRF.df$MAPE <- DescTools::MAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$SMAPE <- DescTools::SMAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$MSE <- DescTools::MSE(loop.lm.mRF)
      loop.k.fold.mRF.df$RMSE <- DescTools::RMSE(loop.lm.mRF)
      
      loop.k.fold.ann.df$Slpe <- loop.lm.ann$coefficients[[2]]
      loop.k.fold.ann.df$Inte <- loop.lm.ann$coefficients[[1]]
      loop.k.fold.ann.df$Rsqr <- summary(loop.lm.ann)$r.square
      loop.k.fold.ann.df$MAE  <- DescTools::MAE(loop.lm.ann)
      loop.k.fold.ann.df$MAPE <- DescTools::MAPE(loop.lm.ann)
      loop.k.fold.ann.df$SMAPE <- DescTools::SMAPE(loop.lm.ann)
      loop.k.fold.ann.df$MSE <- DescTools::MSE(loop.lm.ann)
      loop.k.fold.ann.df$RMSE <- DescTools::RMSE(loop.lm.ann)
      
      k.fold.mRF.df <- rbind(k.fold.mRF.df,loop.k.fold.mRF.df  )
      k.fold.ann.df <- rbind(k.fold.ann.df,loop.k.fold.ann.df  )
      
      
      
      
      
    } #end if
  } #end for params
} #end for k folds

write.csv(k.fold.mRF.df,paste(dump.ml.fld,
                              "S2_mRF_Kfolds_1500_samples.csv",sep = "/"))
write.csv(k.fold.ann.df,paste(dump.ml.fld,
                              "S2_ann_Kfolds_1500_samples.csv",sep = "/"))


train.df <- df.3000
fold.selection <- fold.3000
for (i in 1:nrfolds){
  
  temp.train.df <- train.df[which(fold.selection != i),]
  temp.valid.df <- train.df[which(fold.selection == i),]
  
  
  
  #creating the mRF model
  temp.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                    data=temp.train.df,block.size = 25)
  temp.pred.mRF <- predict(temp.mRF,
                           newdata=temp.valid.df)
  temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))
  
  #creating the ANN model
  temp.X.mat <- as.matrix(temp.train.df[,c(6:ncol(temp.train.df))])
  temp.Y.mat <- as.matrix(temp.train.df[,c(1:5)])
  temp.X.mat.valid <- as.matrix(temp.valid.df[,c(6:ncol(temp.valid.df))])
  
  
  temp.ann <- neuralnetwork(X=temp.X.mat,y=temp.Y.mat,regression=T,
                            hidden.layers =net.st,
                            loss.type = "squared",
                            activ.functions = act.fn,
                            n.epochs = epochs,
                            standardize = T,
                            learn.rates = lrates,
                            optim.type = optype)
  
  #now i could iterate for each trait but im to lazy for that
  temp.ann.pred <- predict(temp.ann,temp.X.mat.valid)
  
  #head(temp.mv.pred.mRF)
  #head(temp.ann.pred)
  
  
  for (j in 1:length(param.names)){
    
    print(paste("Storing statistics for ",
                i,"ith fold of variable",param.names[j]))
    
    
    if ((i == 1) & (j == 1)){
      
      #the creates a data frame to store the outputs
      #j<- 1
      k.fold.mRF.df <- data.frame(kf_th=i,Model="mRF",variable=param.names[j])
      k.fold.ann.df <- data.frame(kf_th=i,Model="ann",variable=param.names[j])
      
      #this fetches the variables
      mRF.pred <- temp.mv.pred.mRF[,j]
      ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      lm.mRF <- lm(mRF.pred ~ temp.valid.df[,j])
      lm.ann <- lm(ann.pred ~ temp.valid.df[,j])
      
      
      temp.mv.pred.mRF[,j]
      
      #stores the data
      k.fold.mRF.df$Slpe <- lm.mRF$coefficients[[2]]
      k.fold.mRF.df$Inte <- lm.mRF$coefficients[[1]]
      k.fold.mRF.df$Rsqr <- summary(lm.mRF)$r.square
      k.fold.mRF.df$MAE  <- DescTools::MAE(lm.mRF)
      k.fold.mRF.df$MAPE <- DescTools::MAPE(lm.mRF)
      k.fold.mRF.df$SMAPE <- DescTools::SMAPE(lm.mRF)
      k.fold.mRF.df$MSE <- DescTools::MSE(lm.mRF)
      k.fold.mRF.df$RMSE <- DescTools::RMSE(lm.mRF)
      
      k.fold.ann.df$Slpe <- lm.ann$coefficients[[2]]
      k.fold.ann.df$Inte <- lm.ann$coefficients[[1]]
      k.fold.ann.df$Rsqr <- summary(lm.ann)$r.square
      k.fold.ann.df$MAE  <- DescTools::MAE(lm.ann)
      k.fold.ann.df$MAPE <- DescTools::MAPE(lm.ann)
      k.fold.ann.df$SMAPE <- DescTools::SMAPE(lm.ann)
      k.fold.ann.df$MSE <- DescTools::MSE(lm.ann)
      k.fold.ann.df$RMSE <- DescTools::RMSE(lm.ann)
      
      
      #plot checking
      #plot(temp.valid.df[,i],temp.valid.df[,i])
      #points(temp.valid.df[,i],mRF.pred,pch=19,col="blue")
      #points(temp.valid.df[,i],ann.pred,pch=19,col="red")
      
      
      
    } else {
      
      #print(j)
      #we have to repeat above and append the dataframe
      
      loop.k.fold.mRF.df <- data.frame(kf_th=i, Model="mRF", variable=param.names[j])
      loop.k.fold.ann.df <- data.frame(kf_th=i, Model="ann", variable=param.names[j])
      
      #this fetches the variables
      loop.mRF.pred <- temp.mv.pred.mRF[,j]
      loop.ann.pred <- temp.ann.pred$predictions[,j]
      
      #fits a models
      loop.lm.mRF <- lm(loop.mRF.pred ~ temp.valid.df[,j])
      loop.lm.ann <- lm(loop.ann.pred ~ temp.valid.df[,j])
      
      #stores the data
      loop.k.fold.mRF.df$Slpe <- loop.lm.mRF$coefficients[[2]]
      loop.k.fold.mRF.df$Inte <- loop.lm.mRF$coefficients[[1]]
      loop.k.fold.mRF.df$Rsqr <- summary(loop.lm.mRF)$r.square
      loop.k.fold.mRF.df$MAE  <- DescTools::MAE(loop.lm.mRF)
      loop.k.fold.mRF.df$MAPE <- DescTools::MAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$SMAPE <- DescTools::SMAPE(loop.lm.mRF)
      loop.k.fold.mRF.df$MSE <- DescTools::MSE(loop.lm.mRF)
      loop.k.fold.mRF.df$RMSE <- DescTools::RMSE(loop.lm.mRF)
      
      loop.k.fold.ann.df$Slpe <- loop.lm.ann$coefficients[[2]]
      loop.k.fold.ann.df$Inte <- loop.lm.ann$coefficients[[1]]
      loop.k.fold.ann.df$Rsqr <- summary(loop.lm.ann)$r.square
      loop.k.fold.ann.df$MAE  <- DescTools::MAE(loop.lm.ann)
      loop.k.fold.ann.df$MAPE <- DescTools::MAPE(loop.lm.ann)
      loop.k.fold.ann.df$SMAPE <- DescTools::SMAPE(loop.lm.ann)
      loop.k.fold.ann.df$MSE <- DescTools::MSE(loop.lm.ann)
      loop.k.fold.ann.df$RMSE <- DescTools::RMSE(loop.lm.ann)
      
      k.fold.mRF.df <- rbind(k.fold.mRF.df,loop.k.fold.mRF.df  )
      k.fold.ann.df <- rbind(k.fold.ann.df,loop.k.fold.ann.df  )
      
      
      
      
      
    } #end if
  } #end for params
} #end for k folds

write.csv(k.fold.mRF.df,paste(dump.ml.fld,
                              "S2_mRF_Kfolds_3000_samples.csv",sep = "/"))
write.csv(k.fold.ann.df,paste(dump.ml.fld,
                              "S2_ann_Kfolds_3000_samples.csv",sep = "/"))


######################################
######## from here we want to summarize the results
#######################################

list.files(dump.ml.fld)
mRF.kfolds.0500 <- read.csv(paste(dump.ml.fld,"S2_mRF_Kfolds_0500_samples.csv",sep="/"))[,-c(1)]
mRF.kfolds.1500 <- read.csv(paste(dump.ml.fld,"S2_mRF_Kfolds_1500_samples.csv",sep="/"))[,-c(1)]
mRF.kfolds.3000 <- read.csv(paste(dump.ml.fld,"S2_mRF_Kfolds_3000_samples.csv",sep="/"))[,-c(1)]
ann.kfolds.0500 <- read.csv(paste(dump.ml.fld,"S2_ann_Kfolds_0500_samples.csv",sep="/"))[,-c(1)]
ann.kfolds.1500 <- read.csv(paste(dump.ml.fld,"S2_ann_Kfolds_1500_samples.csv",sep="/"))[,-c(1)]
ann.kfolds.3000 <- read.csv(paste(dump.ml.fld,"S2_ann_Kfolds_3000_samples.csv",sep="/"))[,-c(1)]


mRF.kfolds.0500$nsamples <- 1
mRF.kfolds.1500$nsamples <- 2
mRF.kfolds.3000$nsamples <- 3
ann.kfolds.0500$nsamples <- 1
ann.kfolds.1500$nsamples <- 2
ann.kfolds.3000$nsamples <- 3

names(mRF.kfolds.0500)
library(reshape2)
#the next function is quite stupid...
#agg.mRF.kfolds.0500 <- aggregate(mRF.kfolds.0500[,c()], by=list(mRF.kfolds.0500$Model,
#                                                          mRF.kfolds.0500$variable,
#                                                          mRF.kfolds.0500$nsamples),FUN=mean)

#lets melt instead, easier for ggplot
mRF.kfolds <- rbind(mRF.kfolds.0500,
                    mRF.kfolds.1500,
                    mRF.kfolds.3000) 
ann.kfolds <- rbind(ann.kfolds.0500,
                    ann.kfolds.1500,
                    ann.kfolds.3000) 

#in the end i should have 8*75 (75rows for each 8 features)
melt.mRF.kfolds <- melt(data = mRF.kfolds,id.vars = c("kf_th","Model","variable","nsamples"))
melt.ann.kfolds <- melt(data = ann.kfolds,id.vars = c("kf_th","Model","variable","nsamples"))

names(melt.mRF.kfolds)[3]<-"Trait"
names(melt.ann.kfolds)[3]<-"Trait"

melt.kfolds <- rbind(melt.mRF.kfolds,melt.ann.kfolds)[,-1]
head(melt.kfolds)

library(dplyr)
8*2*5*3 # 8 features * 2 models * 5 traits  * 3 sample groups
summ.kfolds <- summarize(group_by(melt.kfolds, Model, Trait,nsamples,variable), mean=mean(value), sd=sd(value))
summ.kfolds.df <- as.data.frame(summ.kfolds)

#lets divide everything
unique(melt.kfolds$variable)
#seperation sets
melt.kfolds$ModelSamples <- paste(melt.kfolds$Model,melt.kfolds$nsamples,sep="_")
melt.kfolds$SamplesColors <- NA
melt.kfolds$SamplesColors[melt.kfolds$nsamples == 1] <- "A"
melt.kfolds$SamplesColors[melt.kfolds$nsamples == 2] <- "B"
melt.kfolds$SamplesColors[melt.kfolds$nsamples == 3] <- "C"
melt.kfolds$ModelName <- NA 
melt.kfolds$ModelName[melt.kfolds$Model == "mRF"] <- "Random Forest"
melt.kfolds$ModelName[melt.kfolds$Model == "ann"] <- "Artificial Neural Network"

rs.Slpe.df <- melt.kfolds[melt.kfolds$variable=="Slpe",]
rs.Inte.df <- melt.kfolds[melt.kfolds$variable=="Inte",]
rs.Rsqr.df <- melt.kfolds[melt.kfolds$variable=="Rsqr",]

rs.MAE.df  <- melt.kfolds[melt.kfolds$variable=="MAE",]
rs.MAPE.df <- melt.kfolds[melt.kfolds$variable=="MAPE",]
rs.SMAPE.df <- melt.kfolds[melt.kfolds$variable=="SMAPE",]
rs.MSE.df  <- melt.kfolds[melt.kfolds$variable=="MSE",]
rs.RMSE.df <- melt.kfolds[melt.kfolds$variable=="RMSE",]

library(ggplot2)
library(ggthemes)
model.labels <- c("Random Forest","Artificial Neural Network")

#plot times
#Slope
tiff(paste(dump.ml.fld,"S2_Kfolds_Slope.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))

p<-ggplot(rs.Slpe.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Slope ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)
 
dev.off()

#intercep factor
tiff(paste(dump.ml.fld,"S2_Kfolds_Inter.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.Inte.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Intercept factor ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()



#R squared
tiff(paste(dump.ml.fld,"S2_Kfolds_rsquare.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.Rsqr.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "R-squared ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()


#Mean absolute error - notice that the values actually depend on the dimension of the variable 
tiff(paste(dump.ml.fld,"S2_Kfolds_MAE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.MAE.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Mean absolute error ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()


summ.kfolds.df[1:10,]
#Mean absolute percentage error
tiff(paste(dump.ml.fld,"S2_Kfolds_MAPE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.MAPE.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Mean absolute percentage error ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()

#SMAPE
summ.kfolds.df[1:10,]
tiff(paste(dump.ml.fld,"S2_Kfolds_SMAPE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.SMAPE.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Symmetric mean absolute percentage error ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()

#MSE
summ.kfolds.df[1:10,]
tiff(paste(dump.ml.fld,"S2_Kfolds_MSE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.MSE.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Mean Squared error ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()

#RMSE
summ.kfolds.df[1:10,]
tiff(paste(dump.ml.fld,"S2_Kfolds_RMSE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
p<-ggplot(rs.RMSE.df, aes(x=Trait, y=value,fill=SamplesColors))+ 
  geom_boxplot()+
  scale_fill_discrete(name = "nr samples", labels = c("500", "1500", "3000"))+
  scale_y_continuous(name = "Root Mean Squared error ")+
  facet_wrap(~ModelName)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=16))
plot(p)

dev.off()


##MAKING A TABLE TO SEE THE SUMMARY

#check the summary
dim(summ.kfolds.df)
head(summ.kfolds.df)

final.sum.df  = dcast(summ.kfolds.df, Model+Trait + nsamples ~  variable, value.var = "mean")
final.sum.df.sd  = dcast(summ.kfolds.df, Model+Trait + nsamples ~  variable, value.var = "sd") 
write.csv(final.sum.df ,paste(dump.ml.fld,
                              "S2_Kfolds_Summary_mean.csv",sep = "/"))
write.csv(final.sum.df.sd ,paste(dump.ml.fld,
                              "S2_Kfolds_Summary_std.csv",sep = "/"))


final.sum.df <- data.frame(Trait=NA)
final.sum.df$Trait 

#############################################
#############################################



