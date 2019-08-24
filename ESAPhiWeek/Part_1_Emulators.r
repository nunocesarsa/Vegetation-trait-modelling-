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

#for kfolding
library(dismo)

#for the rmse
library(Metrics)

#for handling data frame more efficiently
library(reshape2)

#group summary
library(nlme)

#for std y column
library(matrixStats)

#for raster stuff
library(raster)
library(maptools)

par(mfrow=c(1,1))
gc()

setwd("D:/ESAPhiWeek")

#first we set a bunch of initial parameters we want to predict
#these are loosely based on:
#https://www.sciencedirect.com/science/article/pii/S0303243415000100
#https://www.mdpi.com/2072-4292/11/13/1597/htm
#https://www.mdpi.com/2072-4292/8/2/119

#its a 5 parameter estimate
param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
                         5,80, #Cab
                         0,25, #Car
                         0.005,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.02, #Cm
                         0.1,8), #LAI
                         #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)


#lets generate a space of biophysical parameters and their spectral response with prosail
prosail.runs <- 1000*nrow(param.maxmin) #1000 per each trait

#using this, we generate a LHS of the trait space
LHS <- Latinhyper(param.maxmin,prosail.runs)

param.list <- data.frame(Cab=LHS[,1],
                         Car=LHS[,2],
                         Cw=LHS[,3],
                         Cm=LHS[,4],
                         LAI=LHS[,5])


#we accept all other prosail parameters as default and we convert out spectral response to S resolution
mRF.spclib <- PROSAIL(parameterList = param.list)
mRF.s2.spclib <- spectralResampling(mRF.spclib,
                                    "Sentinel2",response_function = TRUE)
#we should remove the bands not usually used in S2 (only if offered up to 20m res) 
#b 2,3,4,5,6,7,8a,11,12
mRF.s2.spc.sel <- mRF.s2.spclib[,c(2,3,4,5,6,7,9,12,13)]

#we create a DF with everything
train.df <- cbind(param.list,as.data.frame(mRF.s2.spc.sel))

#we should also rename the bands to something proper
names(train.df) <- c(names(train.df)[1:5],
                     "B02","B03","B04",
                     "B05","B06","B07",
                     "B8A","B11","B12")

#self testing the model
#First we K fold the model
set.seed(1000) #this forces the same starting point

nrfolds <- 10 #a number of folds to use
fold.selection <- kfold(train.df,nrfolds)

table(fold.selection) # just to check how much data in each fold
#which(fold.selection != 1) #we use this to find each folder

#also a quicker function to get the r2
rsq <- function (x, y) cor(x, y) ^ 2

for (i in 1:nrfolds){
  
  temp.train.df <- train.df[which(fold.selection != i),]
  temp.valid.df <- train.df[which(fold.selection == i),]
  
  #creating the mRF model
  temp.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                    data=temp.train.df,block.size = 50)
  
  temp.pred.mRF <- predict(temp.mRF,
                           newdata=temp.valid.df)
  
  temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))
  
  #now i could iterate for each trait but im to lazy for that
  
  
  for (j in 1:length(names(param.list))){
    print(paste("Storing r2 and rmse for the",
                i,"ith fold of variable",names(param.list)[j]))
    
    temp.r2 <- rsq(x = temp.valid.df[,j],
                   y = temp.mv.pred.mRF[,j])
    
    temp.RMSE <- rmse(actual = temp.valid.df[,j],
                      predicted = temp.mv.pred.mRF[,j])
    
    if ((i == 1) & (j == 1)){
      #the first iteration creates a df to hold the output of each k fold test
      
      
      k.fold.test.df <- data.frame(kf_th=i,variable=names(param.list)[j],
                                   rsquare=temp.r2,RMSE=temp.RMSE)
      
    } else {
      #while all other iterations just append the result to the end
      temp.df<- data.frame(kf_th=i, variable=names(param.list)[j],
                           rsquare=temp.r2,RMSE=temp.RMSE)
      
      k.fold.test.df <- rbind(k.fold.test.df,temp.df)
      
    }
    

    
  }
}

#once this is done, we need summarize by each k-fold
names(k.fold.test.df) <- c("kf_th","biotrait","rsquare" ,"RMSE" )

k.fold.test.df.redux <- k.fold.test.df[,c(2,3,4)]

aggregate(k.fold.test.df.redux, by=list(k.fold.test.df.redux$biotrait), 
          FUN=mean, na.rm=TRUE)
library(nlme)
k.fold.summary <- gsummary(k.fold.test.df.redux,FUN=mean,groups=k.fold.test.df.redux$biotrait)

#are these good RMSE? - depends on the scale of my data actually

bb <- colMeans(param.list)
avg.trait.val <- colMeans(param.list)

k.fold.summary$TraitAvg <- avg.trait.val
k.fold.summary$TraitStd <- colSds(x = as.matrix(param.list))
k.fold.summary$TraitMin <- param.maxmin[,2]
k.fold.summary$TraitMax <- param.maxmin[,1]
#k.fold.summary$RMSE_to_max <- k.fold.summary$RMSE/param.maxmin[,2]
#k.fold.summary$RMSE_to_min <- k.fold.summary$RMSE/param.maxmin[,1]
#k.fold.summary$RMSE_to_midpoint <- k.fold.summary$RMSE/(param.maxmin[,2]-param.maxmin[,1])
#k.fold.summary$RMSE2AvgTrait <- k.fold.summary$RMSE/avg.trait.val

k.fold.summary

#We can now try generate uncertainity maps bsed on the model. 

#step one, we retrain a model with all the data
#step two, we plot the the model against another dataset
#we use the error by val to generate uncertainity

#training the full model
full.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                  data=train.df,block.size = 50)


set.seed(10000)

#generating a test dataset
param.maxmin.unc <- param.maxmin

#lets generate a space of biophysical parameters and their spectral response with prosail
prosail.runs.unc <- 1000*nrow(param.maxmin.unc) #1000 per each trait

#using this, we generate a LHS of the trait space
LHS.unc <- Latinhyper(param.maxmin.unc,prosail.runs.unc)

param.list.unc <- data.frame(Cab=LHS[,1],
                             Car=LHS[,2],
                             Cw=LHS[,3],
                             Cm=LHS[,4],
                             LAI=LHS[,5])


#generating a test dataset prosail output
mRF.spclib.unc <- PROSAIL(parameterList = param.list.unc)
mRF.s2.spclib.unc <- spectralResampling(mRF.spclib.unc,
                                    "Sentinel2",response_function = TRUE)

#removing the bands as before
#b 2,3,4,5,6,7,8a,11,12
mRF.s2.spc.sel.unc <- mRF.s2.spclib.unc[,c(2,3,4,5,6,7,9,12,13)]

unc.df <- cbind(param.list.unc,as.data.frame(mRF.s2.spc.sel.unc))
#we should also rename the bands to something proper
names(unc.df) <- c(names(unc.df)[1:5],
                     "B02","B03","B04",
                     "B05","B06","B07",
                     "B8A","B11","B12")

pred.mRF.unc <- predict(full.mRF,
                        newdata=unc.df)

mv.pred.mRF.unc <- as.data.frame(get.mv.predicted(pred.mRF.unc))

#lets bring things together
names(mv.pred.mRF.unc)<- paste("pred",names(mv.pred.mRF.unc),sep="_")
unc.df <- cbind(param.list.unc,mv.pred.mRF.unc)

par(mfrow=c(3,2))
for(i in 1:ncol(param.list.unc)){
  
  print(i)
  plot(param.list.unc[,i],
       mv.pred.mRF.unc[,i],
       xlab=names(unc.df)[i],
       ylab=names(mv.pred.mRF.unc)[i])
}

#talk with mitra - can we use the CI of the linear model as an estimate of uncertainity


#using THE CONFIDENCE INTERVAL to create a raster prediction and its uncertainities

#first we generate 100x100x5 object where each ith z dimension is one of the parameters
set.seed(5000)

nl <- 25
nc <- 25
m.Cab <- matrix(rtnorm(nl*nc,mean=55,sd=20, lower=param.maxmin[1,1],upper=param.maxmin[1,2]),nl,nc)
m.Car <- matrix(rtnorm(nl*nc,mean=12.5,sd=7,lower=param.maxmin[2,1],upper=param.maxmin[2,2]),nl,nc)
m.Cw  <- matrix(rtnorm(nl*nc,mean=0.013,sd=0.005,lower=param.maxmin[3,1],upper=param.maxmin[3,2]),nl,nc)
m.Cm  <- matrix(rtnorm(nl*nc,mean=0.013,sd=0.005,lower=param.maxmin[4,1],upper=param.maxmin[4,2]),nl,nc)
m.LAI <- matrix(rtnorm(nl*nc,mean=4,sd=2.5,lower=param.maxmin[5,1],upper=param.maxmin[5,2]),nl,nc)

#and a place to receive the predictions
m.pred.Cab <- m.Cab*0
m.pred.Car <- m.Cab*0
m.pred.Cw  <- m.Cab*0
m.pred.Cm  <- m.Cab*0
m.pred.LAI <- m.Cab*0
#and a place to receive uncertanities
m.unc.Cab <- m.Cab*0
m.unc.Car <- m.Cab*0
m.unc.Cw  <- m.Cab*0
m.unc.Cm  <- m.Cab*0
m.unc.LAI <- m.Cab*0

#then we create a linear model for each simulated prosail vs mRF (simulated)
lm.Cab <- lm(pred_Cab~Cab,data=unc.df)
lm.Car <- lm(pred_Car~Car,data=unc.df)
lm.Cw  <- lm(pred_Cw~ Cw, data=unc.df)
lm.Cm  <- lm(pred_Cm~ Cm, data=unc.df)
lm.LAI <- lm(pred_LAI~LAI,data=unc.df)

#once we have these matrices, we can run prosail on each col line.
for (i in 1:ncol(m.Cab)){
  
  for (j in 1:nrow(m.Cab)){
    print(paste("Predicting: col",i,"line",j ))
    
    m.param.list <- data.frame(Cab=m.Cab[j,i],Car=m.Car[j,i],
                               Cw=m.Cw[j,i],Cm=m.Cm[j,i],
                               LAI=m.LAI[j,i])
    
    #generating a test dataset prosail output
    m.spclib.unc <- PROSAIL(parameterList = m.param.list)
    m.s2.spclib.unc <- spectralResampling(m.spclib.unc,
                                            "Sentinel2",response_function = TRUE)
    #removing the bands as before
    #b 2,3,4,5,6,7,8a,11,12
    m.s2.spc.sel.unc <- m.s2.spclib.unc[,c(2,3,4,5,6,7,9,12,13)]
    
    m.df <- cbind(m.param.list,as.data.frame(m.s2.spc.sel.unc))
    
    names(m.df) <- c(names(m.df)[1:5],
                       "B02","B03","B04",
                       "B05","B06","B07",
                       "B8A","B11","B12")
    
    m.pred.mRF.unc <- predict(full.mRF,
                              newdata=m.df)
    
    m.mv.pred.mRF.unc <- as.data.frame(get.mv.predicted(m.pred.mRF.unc))
    
    m.pred.Cab[j,i] <- m.mv.pred.mRF.unc[1,1]
    m.pred.Car[j,i] <- m.mv.pred.mRF.unc[1,2]
    m.pred.Cw[j,i]  <- m.mv.pred.mRF.unc[1,3]
    m.pred.Cm[j,i]  <- m.mv.pred.mRF.unc[1,4]
    m.pred.LAI[j,i] <- m.mv.pred.mRF.unc[1,5]
    
    #now we hve the predictions and the original values, lets see
    
    temp.test.Cab <- predict(lm.Cab,newdata=data.frame(Cab=m.pred.Cab[j,i]),interval="prediction")
    temp.test.Car <- predict(lm.Car,newdata=data.frame(Car=m.pred.Car[j,i]),interval="prediction")
    temp.test.Cw  <- predict(lm.Cw, newdata=data.frame(Cw=m.pred.Cw[j,i  ]),interval="prediction")
    temp.test.Cm  <- predict(lm.Cm, newdata=data.frame(Cm=m.pred.Cm[j,i  ]),interval="prediction")
    temp.test.LAI <- predict(lm.LAI,newdata=data.frame(LAI=m.pred.LAI[j,i]),interval="prediction")
    
    m.unc.Cab[j,i] <- temp.test.Cab[,3]-temp.test.Cab[,2]
    m.unc.Car[j,i] <- temp.test.Car[,3]-temp.test.Car[,2]
    m.unc.Cw[j,i]  <- temp.test.Cw[,3] -temp.test.Cw[,2 ]
    m.unc.Cm[j,i]  <- temp.test.Cm[,3] -temp.test.Cm[,2 ]
    m.unc.LAI[j,i] <- temp.test.LAI[,3]-temp.test.LAI[,2]
    

     
  }
}

#once all the predictions are done, lets create an uncertainity rst using our original full model vs its training



temp.cab.rst <- raster(m.Cab)
temp.cab.pred.rst <- raster(m.pred.Cab)
temp.cab.unc.rst <- raster(m.unc.Cab)
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.cab.rst,main="Cab - original")
plot(temp.cab.pred.rst, main="Cab - predicted")
plot(abs(temp.cab.rst-temp.cab.pred.rst),main="Absolute difference")
plot(temp.cab.unc.rst,main="CI of the prediction")
#plot((temp.cab.unc.rst-cellStats(temp.cab.unc.rst,min))/
       #(cellStats(temp.cab.unc.rst,max)-cellStats(temp.cab.unc.rst,min)))

temp.Car.rst <- raster(m.Car)
temp.Car.pred.rst <- raster(m.pred.Car)
temp.Car.unc.rst <- raster(m.unc.Car)
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Car.rst,main="Car - original")
plot(temp.Car.pred.rst, main="Car - predicted")
plot(abs(temp.Car.rst-temp.Car.pred.rst),main="Absolute difference")
plot(temp.Car.unc.rst,main="CI of the prediction")

temp.Cw.rst <- raster(m.Cw)
temp.Cw.pred.rst <- raster(m.pred.Cw)
temp.Cw.unc.rst <- raster(m.unc.Cw)
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Cw.rst,main="Cw - original")
plot(temp.Cw.pred.rst, main="Cw - predicted")
plot(abs(temp.Cw.rst-temp.Cw.pred.rst),main="Absolute difference")
plot(temp.Cw.unc.rst,main="CI of the prediction")

temp.Cm.rst <- raster(m.Cm)
temp.Cm.pred.rst <- raster(m.pred.Cm)
temp.Cm.unc.rst <- raster(m.unc.Cm)
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Cm.rst,main="Cm - original")
plot(temp.Cm.pred.rst, main="Cm - predicted")
plot(abs(temp.Cm.rst-temp.Cm.pred.rst),main="Absolute difference")
plot(temp.Cm.unc.rst,main="CI of the prediction")


temp.LAI.rst <- raster(m.LAI)
temp.LAI.pred.rst <- raster(m.pred.LAI)
temp.LAI.unc.rst <- raster(m.unc.LAI)
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.LAI.rst,main="LAI - original")
plot(temp.LAI.pred.rst, main="LAI - predicted")
plot(abs(temp.LAI.rst-temp.LAI.pred.rst),main="Absolute difference")
plot(temp.LAI.unc.rst,main="CI of the prediction")

#STOP FOR NOW
