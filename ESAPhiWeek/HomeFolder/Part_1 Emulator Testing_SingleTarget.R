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
#

#change by joris
param.maxmin.t <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  .999, exp(-120/100), #Cab
  .999, exp(-20/100), #Car
  .999,0.01, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  .999,0.01, #Cm
  .999,0.01), #LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

#lets generate a space of biophysical parameters and their spectral response with prosail
prosail.runs <- 100*nrow(param.maxmin) #1000 per each trait

LHS <- Latinhyper(param.maxmin.t,prosail.runs)
head(LHS)

param.maxmin[1,] <- -100*log(param.maxmin.t[1,])
param.maxmin[2,] <- -100*log(param.maxmin.t[2,])
param.maxmin[3,] <- -(1/50)*log(param.maxmin.t[3,])
param.maxmin[4,] <- -(1/100)*log(param.maxmin.t[4,])
param.maxmin[5,] <- -(2)*log(param.maxmin.t[5,])
param.maxmin


LHS.test <- Latinhyper(param.maxmin.t,prosail.runs*5)




#using this, we generate a LHS of the trait space
#LHS <- Latinhyper(param.maxmin,prosail.runs)


param.list <- data.frame(Cab=-100*log(LHS[,1]),
                         Car=-100*log(LHS[,2]),
                         Cw=-(1/50)*log(LHS[,3]),
                         Cm=-(1/100)*log(LHS[,4]),
                         LAI=-(2)*log(LHS[,5]))
param.list.test <- data.frame(Cab=-100*log(LHS.test[,1]),
                         Car=-100*log(LHS.test[,2]),
                         Cw=-(1/50)*log(LHS.test[,3]),
                         Cm=-(1/100)*log(LHS.test[,4]),
                         LAI=-(2)*log(LHS.test[,5]))


head(param.list)

#we accept all other prosail parameters as default and we convert out spectral response to S resolution
mRF.spclib <- PROSAIL(parameterList = param.list)
mRF.s2.spclib <- spectralResampling(mRF.spclib,
                                    "Sentinel2",response_function = TRUE)

mRF.spclib.test <- PROSAIL(parameterList = param.list.test)
mRF.s2.spclib.test <- spectralResampling(mRF.spclib.test,
                                    "Sentinel2",response_function = TRUE)


plot(mRF.spclib)
plot(mRF.s2.spclib,ylim=c(0,0.05))
#we should remove the bands not usually used in S2 (only if offered up to 20m res) 
#b 2,3,4,5,6,7,8a,11,12
mRF.s2.spc.sel <- mRF.s2.spclib[,c(2,3,4,5,6,7,9,12,13)]
mRF.s2.spc.sel.test <- mRF.s2.spclib.test[,c(2,3,4,5,6,7,9,12,13)]

#we create a DF with everything
train.df <- cbind(param.list,as.data.frame(mRF.s2.spc.sel))
train.df.test <- cbind(param.list.test,as.data.frame(mRF.s2.spc.sel.test))

#we should also rename the bands to something proper
names(train.df) <- c(names(train.df)[1:5],
                     "B02","B03","B04",
                     "B05","B06","B07",
                     "B8A","B11","B12")

#names(train.df.test) <- c(names(train.df.test)[1:5],
#                     "B02","B03","B04",
#                     "B05","B06","B07",
#                     "B8A","B11","B12")

#self testing the model
#First we K fold the model
set.seed(1000) #this forces the same starting point

nrfolds <- 10 #a number of folds to use
fold.selection <- kfold(train.df,nrfolds)

table(fold.selection) # just to check how much data in each fold
#which(fold.selection != 1) #we use this to find each folder

#also a quicker function to get the r2
rsq <- function (x, y) cor(x, y) ^ 2

train.df.original <- train.df


#do this section if you ant to transform your variables
head(train.df)
train.df$Cab <- exp(-train.df$Cab/100)
boxplot(train.df$Cab)
train.df$Car <- exp(-train.df$Car/100)
boxplot(train.df$Car)
train.df$LAI <- exp(-train.df$LAI/2)
boxplot(train.df$LAI)
train.df$Cm <- exp(-100*train.df$Cm)
boxplot(train.df$Cm)
train.df$Cw <- exp(-50*train.df$Cw)
boxplot(train.df$Cw)

head(train.df)

for (i in 1:nrfolds){
  
  blk.size = 1
  
  Cab.temp.train.df <- train.df[which(fold.selection != i),-c(2,3,4,5)]
  Cab.temp.valid.df <- train.df[which(fold.selection == i),-c(2,3,4,5)]
  
  Car.temp.train.df <- train.df[which(fold.selection != i),-c(1,3,4,5)]
  Car.temp.valid.df <- train.df[which(fold.selection == i),-c(1,3,4,5)]
  
  Cw.temp.train.df <- train.df[which(fold.selection != i),-c(1,2,4,5)]
  Cw.temp.valid.df <- train.df[which(fold.selection == i),-c(1,2,4,5)]
  
  Cm.temp.train.df <- train.df[which(fold.selection != i),-c(1,2,3,5)]
  Cm.temp.valid.df <- train.df[which(fold.selection == i),-c(1,2,3,5)]
  
  LAI.temp.train.df <- train.df[which(fold.selection != i),-c(1,2,3,4)]
  LAI.temp.valid.df <- train.df[which(fold.selection == i),-c(1,2,3,4)] 
  
  temp.train.df <- train.df[which(fold.selection != 1),]
  temp.valid.df <- train.df[which(fold.selection == i),] 
  
  #creating the mRF model - MULTI TASK MODEL 
  temp.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                    data=temp.train.df,block.size = blk.size)
  
  temp.pred.mRF <- predict(temp.mRF,
                           newdata=temp.valid.df)
 
  temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))
  
  
  #creating the mRF model - SINGLE TASK
  Cab.temp.mRF <- rfsrc(Multivar(Cab)~.,
                        data=Cab.temp.train.df,block.size = blk.size)
  
  Cab.temp.pred.mRF <- predict(Cab.temp.mRF,
                               newdata=Cab.temp.valid.df)
  
  Cab.temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(Cab.temp.pred.mRF))
  
  
  Car.temp.mRF <- rfsrc(Multivar(Car)~.,
                        data=Car.temp.train.df,block.size = blk.size)
  
  Car.temp.pred.mRF <- predict(Car.temp.mRF,
                               newdata=Car.temp.valid.df)
  
  Car.temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(Car.temp.pred.mRF))  
  
  
  Cw.temp.mRF <- rfsrc(Multivar(Cw)~.,
                        data=Cw.temp.train.df,block.size =blk.size)
  
  Cw.temp.pred.mRF <- predict(Cw.temp.mRF,
                               newdata=Cw.temp.valid.df)
  
  Cw.temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(Cw.temp.pred.mRF))  
  
  
  Cm.temp.mRF <- rfsrc(Multivar(Cm)~.,
                       data=Cm.temp.train.df,block.size = blk.size)
  
  Cm.temp.pred.mRF <- predict(Cm.temp.mRF,
                              newdata=Cm.temp.valid.df)
  
  Cm.temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(Cm.temp.pred.mRF))  
  
  
  LAI.temp.mRF <- rfsrc(Multivar(LAI)~.,
                       data=LAI.temp.train.df,block.size = blk.size)
  
  LAI.temp.pred.mRF <- predict(LAI.temp.mRF,
                              newdata=LAI.temp.valid.df)
  
  LAI.temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(LAI.temp.pred.mRF))  
  

  #now i could iterate for each trait but im to lazy for that - I ACTUALLY WASN'T IT SEEMS LOL
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
  
  
  
  
  
  #lets try to adapt to store the single target predictions:
  #well, if we adapt our output structure to be the same as the multivariate then we dont need to change almost anything
  
  
  st.temp.valid.df <- temp.valid.df 
  st.temp.mv.pred.mRF <- cbind(Cab.temp.mv.pred.mRF,
                               Car.temp.mv.pred.mRF,
                               Cw.temp.mv.pred.mRF,
                               Cm.temp.mv.pred.mRF,
                               LAI.temp.mv.pred.mRF)

  for (k in 1:length(names(param.list))){
    print(paste("Single Target - Storing r2 and rmse for the",
                i,"ith fold of variable",names(param.list)[k]))
    
    st.temp.r2 <- rsq(x = st.temp.valid.df[,k],
                   y = st.temp.mv.pred.mRF[,k])
    
    st.temp.RMSE <- rmse(actual = st.temp.valid.df[,k],
                      predicted = st.temp.mv.pred.mRF[,k])
    
    if ((i == 1) & (k == 1)){
      #the first iteration creates a df to hold the output of each k fold test
      
      
      k.fold.test.df.st <- data.frame(kf_th=i,variable=names(param.list)[k],
                                   rsquare=st.temp.r2,RMSE=st.temp.RMSE)
      
    } else {
      #while all other iterations just append the result to the end
      temp.df.st <- data.frame(kf_th=i, variable=names(param.list)[k],
                           rsquare=st.temp.r2,RMSE=st.temp.RMSE)
      
      k.fold.test.df.st <- rbind(k.fold.test.df.st,temp.df.st)
      
    }
    
    
    
  }
  
}

j
i

#once this is done, we need summarize by each k-fold
names(k.fold.test.df) <- c("kf_th","biotrait","rsquare" ,"RMSE" )
names(k.fold.test.df.st) <- c("kf_th","biotrait","rsquare" ,"RMSE" )
k.fold.test.df.redux <- k.fold.test.df[,c(2,3,4)]
k.fold.test.df.redux.st <- k.fold.test.df.st[,c(2,3,4)]

library(nlme)
k.fold.summary <- gsummary(k.fold.test.df.redux,FUN=mean,groups=k.fold.test.df.redux$biotrait)
avg.trait.val <- colMeans(train.df[,c(1:5)])
k.fold.summary$TraitAvg <- avg.trait.val
k.fold.summary$TraitStd <- colSds(x = as.matrix(param.list))
k.fold.summary$TraitMin <- param.maxmin[,2]
k.fold.summary$TraitMax <- param.maxmin[,1]

k.fold.summary.st <- gsummary(k.fold.test.df.redux.st,FUN=mean,groups=k.fold.test.df.redux.st$biotrait)
k.fold.summary.st$TraitAvg <- avg.trait.val
k.fold.summary.st$TraitStd <- colSds(x = as.matrix(param.list))
k.fold.summary.st$TraitMin <- param.maxmin[,2]
k.fold.summary.st$TraitMax <- param.maxmin[,1]

k.fold.summary
k.fold.summary.st


train.df.test$Cab <- exp(-train.df.test$Cab/100)
boxplot(train.df.test$Cab)
train.df.test$Car <- exp(-train.df.test$Car/100)
boxplot(train.df.test$Car)
train.df.test$LAI <- exp(-train.df.test$LAI/2)
boxplot(train.df.test$LAI)
train.df.test$Cm <- exp(-100*train.df.test$Cm)
boxplot(train.df.test$Cm)
train.df.test$Cw <- exp(-50*train.df.test$Cw)
boxplot(train.df.test$Cw)


final.test <- train.df.test[,c(1:5)]

head(train.df.test)
temp.mRF$family
bbb.temp.pred.mRF <- predict(temp.mRF,
                             newdata=train.df.test)

bbb.pred <- get.mv.predicted(bbb.temp.pred.mRF )

head(final.test)
head(bbb.pred)
final.test <- cbind(final.test,bbb.pred)
head(final.test)

rsq(final.test[,3],final.test[,8])
mse(final.test[,3],final.test[,8])
plot(final.test[,3],final.test[,8])





k.fold.summary.st > k.fold.summary


bb <- tune.rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
           data=temp.train.df)

cc <- tune.rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                data=temp.train.df,doBest = T)

#implement a mean percentage error

#are the


aggregate(k.fold.test.df.redux, by=list(k.fold.test.df.redux$biotrait), 
          FUN=mean, na.rm=TRUE)




names(param.list)[k]
rsq(x = st.temp.valid.df[,k],
    y = st.temp.mv.pred.mRF[,k])

rsq(x = temp.valid.df[,j],
    y = temp.mv.pred.mRF[,j])

rmse(actual = st.temp.valid.df[,k],
     predicted = st.temp.mv.pred.mRF[,k])

rmse(actual = temp.valid.df[,j],
     predicted = temp.mv.pred.mRF[,j])

#tsting
temp.train.df <- train.df[which(fold.selection != 1),-c(2,3,4,5)]
temp.valid.df <- train.df[which(fold.selection == 1),-c(2,3,4,5)] 
temp.mRF <- rfsrc(Multivar(Cab)~.,
                  data=temp.train.df,block.size = 50)

temp.pred.mRF <- predict(temp.mRF,
                         newdata=temp.valid.df)

temp.mv.pred.mRF <- as.data.frame(get.mv.predicted(temp.pred.mRF))

names(train.df)