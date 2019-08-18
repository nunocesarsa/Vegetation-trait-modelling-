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


names(unc.df)
test.lm <- lm(pred_Cab~Cab,data=unc.df)
par(mfrow=c(1,1))
plot(param.list.unc[,1],
     mv.pred.mRF.unc[,1],
     xlab=names(unc.df)[1],
     ylab=names(mv.pred.mRF.unc)[1])

plot(test.lm)

#this works to produce a confidence interval but its not that trustworthy.. for me..
test.lm <- lm(pred_Cab~Cab,data=unc.df)
newdata=data.frame(Cab=c(-10,20,30,40,50,70,100,200,300,1000))
predict(test.lm, 
        newdata, 
        interval="prediction") 

cc <- predict(test.lm, 
              newdata, 
              interval="prediction") 
cc[,3]-cc[,2]
