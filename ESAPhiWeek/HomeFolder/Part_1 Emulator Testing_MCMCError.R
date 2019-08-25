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
write.csv(k.fold.summary,"10Fold_crossval.csv")
write.csv2(k.fold.summary,"10Fold_crossval_csv2.csv")

#We can now try generate uncertainity maps bsed on the model. 

#step one, we retrain a model with all the data
#step two, we plot the the model against another dataset
#we use the error by val to generate uncertainity

#training the full model
full.mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                  data=train.df,block.size = 50)

plot(full.mRF,m.target="Cab")
plot(full.mRF,m.target="Car")
plot(full.mRF,m.target="Cw")
plot(full.mRF,m.target="Cm")
plot(full.mRF,m.target="LAI")

#causes exception error and full crash of R
#plot(quantreg(Multivar(Cab,Car,Cw,Cm,LAI)~.,
#              data=train.df,))

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


mean_stdE <- function(x) sd(x)/sqrt(length(x))


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
    
    #then we create a linear model for each simulated prosail vs mRF (simulated)
    #nr of runs of k - like a monte carlo
    kruns <- 50
    submodel.size <- 100
    temp.test.Cab <- c(1:kruns)
    temp.test.Car <- c(1:kruns)
    temp.test.Cw  <- c(1:kruns)
    temp.test.Cm  <- c(1:kruns)
    temp.test.LAI <- c(1:kruns)
    
    for (k in 1:kruns){
      
      unc.df.redux <- unc.df[sample(1:nrow(unc.df),submodel.size),]
      
      
      lm.Cab <- lm(pred_Cab~Cab,data=unc.df.redux )
      lm.Car <- lm(pred_Car~Car,data=unc.df.redux )
      lm.Cw  <- lm(pred_Cw~ Cw, data=unc.df.redux )
      lm.Cm  <- lm(pred_Cm~ Cm, data=unc.df.redux )
      lm.LAI <- lm(pred_LAI~LAI,data=unc.df.redux )
      
      temp.test.Cab[k] <- predict(lm.Cab,newdata=data.frame(Cab=m.pred.Cab[j,i]))#,interval="prediction")
      temp.test.Car[k] <- predict(lm.Car,newdata=data.frame(Car=m.pred.Car[j,i]))#,interval="prediction")
      temp.test.Cw[k]  <- predict(lm.Cw, newdata=data.frame(Cw=m.pred.Cw[j,i  ]))#,interval="prediction")
      temp.test.Cm[k]  <- predict(lm.Cm, newdata=data.frame(Cm=m.pred.Cm[j,i  ]))#,interval="prediction")
      temp.test.LAI[k] <- predict(lm.LAI,newdata=data.frame(LAI=m.pred.LAI[j,i]))#,interval="prediction")
      

    }
  
    
    m.unc.Cab[j,i] <- mean_stdE(temp.test.Cab)
    m.unc.Car[j,i] <- mean_stdE(temp.test.Car)
    m.unc.Cw[j,i]  <- mean_stdE(temp.test.Cw)
    m.unc.Cm[j,i]  <- mean_stdE(temp.test.Cm)
    m.unc.LAI[j,i] <- mean_stdE(temp.test.LAI)
    

     
  }
}

warnings()
#once all the predictions are done, lets create an uncertainity rst using our original full model vs its training



### Cab
temp.cab.rst <- raster(m.Cab)
temp.cab.pred.rst <- raster(m.pred.Cab)
temp.cab.unc.rst <- raster(m.unc.Cab)

names(cab.stack) <- c("Cab","Prediction","Abs_dif","MC StdE")
maxv <- ceiling(maxValue(temp.cab.rst)+1)
minv <- 0
brks <- seq(minv,maxv,by=(maxv-minv)/10)
nbrks <- length(brks)-1
r.range <- c(minv, maxv)
#colfunc<-colorRampPalette(c("red","yellow","royalblue","darkgreen"))
#colfunc<-colorRampPalette(c("springgreen", "royalblue", "yellow", "red"))
colfunc<-colorRampPalette(c("burlywood1","orange","green"))
                          
par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.cab.rst,breaks=brks,col=colfunc(nbrks), legend = F, zlim=c(minv,maxv),main="Cab - original")
plot(temp.cab.pred.rst,breaks=brks,col=colfunc(nbrks), legend = T, zlim=c(minv,maxv), main="Cab - predicted")
plot(abs(temp.cab.rst-temp.cab.pred.rst),main="Absolute difference")
plot(temp.cab.unc.rst,main="stdE of model mean")

# generate palette
# generate palette
#colfunc<-colorRampPalette(c("springgreen", "royalblue", "yellow", "red"))

plot(cab.stack)

### Car
temp.Car.rst <- raster(m.Car)
temp.Car.pred.rst <- raster(m.pred.Car)
temp.Car.unc.rst <- raster(m.unc.Car)

names(Car.stack) <- c("Car","Prediction","Abs_dif","MC StdE")
maxv <- ceiling(maxValue(temp.Car.rst)+1)
minv <- 0
brks <- seq(minv,maxv,by=(maxv-minv)/10)
nbrks <- length(brks)-1
r.range <- c(minv, maxv)
colfunc<-colorRampPalette(c("burlywood1","orange","green"))
#colfunc<-colorRampPalette(c("red","yellow","royalblue","darkgreen"))
#colfunc<-colorRampPalette(c("springgreen", "royalblue", "yellow", "red"))


par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Car.rst,breaks=brks,col=colfunc(nbrks), legend =F, zlim=c(minv,maxv),main="Car - original")
plot(temp.Car.pred.rst,breaks=brks,col=colfunc(nbrks), legend = T, zlim=c(minv,maxv), main="Car - predicted")
plot(abs(temp.Car.rst-temp.Car.pred.rst),main="Absolute difference")
plot(temp.Car.unc.rst,main="stdE of model mean")



####Cw
temp.Cw.rst <- raster(m.Cw)
temp.Cw.pred.rst <- raster(m.pred.Cw)
temp.Cw.unc.rst <- raster(m.unc.Cw)

maxv <- maxValue(temp.Cw.rst)
minv <- 0
brks <- seq(minv,maxv,by=(maxv-minv)/10)
nbrks <- length(brks)-1
r.range <- c(minv, maxv)
colfunc<-colorRampPalette(c("burlywood1","orange","green"))

par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Cw.rst,breaks=brks,col=colfunc(nbrks), legend = F, zlim=c(minv,maxv),main="Cw - original")
plot(temp.Cw.pred.rst,breaks=brks,col=colfunc(nbrks), legend = T, zlim=c(minv,maxv), main="Cw - predicted")
plot(abs(temp.Cw.rst-temp.Cw.pred.rst),main="Absolute difference")
plot(temp.Cw.unc.rst,main="stdE of model mean")


####Cm
temp.Cm.rst <- raster(m.Cm)
temp.Cm.pred.rst <- raster(m.pred.Cm)
temp.Cm.unc.rst <- raster(m.unc.Cm)

maxv <- maxValue(temp.Cw.rst)
minv <- 0
brks <- seq(minv,maxv,by=(maxv-minv)/10)
nbrks <- length(brks)-1
r.range <- c(minv, maxv)
colfunc<-colorRampPalette(c("burlywood1","orange","green"))

par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.Cm.rst,breaks=brks,col=colfunc(nbrks), legend = F, zlim=c(minv,maxv),main="Cm - original")
plot(temp.Cm.pred.rst,breaks=brks,col=colfunc(nbrks), legend = T, zlim=c(minv,maxv), main="Cm - predicted")
plot(abs(temp.Cm.rst-temp.Cm.pred.rst),main="Absolute difference")
plot(temp.Cm.unc.rst,main="stdE of model mean")



temp.LAI.rst <- raster(m.LAI)
temp.LAI.pred.rst <- raster(m.pred.LAI)
temp.LAI.unc.rst <- raster(m.unc.LAI)

maxv <- ceiling(maxValue(temp.LAI.rst)+1)
minv <- 0
brks <- seq(minv,maxv,by=(maxv-minv)/10)
nbrks <- length(brks)-1
r.range <- c(minv, maxv)
colfunc<-colorRampPalette(c("burlywood1","orange","green"))



par(mfrow=c(2,2)) #THIS IS A GOOD IDEA TO VISUALIZE THE DIFFERENCES AND RELATE THAT TO THE PARAMETERS
plot(temp.LAI.rst,breaks=brks,col=colfunc(nbrks), legend = F, zlim=c(minv,maxv),main="LAI - original")
plot(temp.LAI.pred.rst,breaks=brks,col=colfunc(nbrks), legend = T, zlim=c(minv,maxv), main="LAI - predicted")
plot(abs(temp.LAI.rst-temp.LAI.pred.rst),main="Absolute difference")
plot(temp.LAI.unc.rst,main="stdE of model mean")


#plot((temp.cab.unc.rst-cellStats(temp.cab.unc.rst,min))/
       #(cellStats(temp.cab.unc.rst,max)-cellStats(temp.cab.unc.rst,min)))








#STOP FOR NOW








m.pred.mRF.unc

#testing another way to estimate the uncertainity

#the idea is to produce RMSE of a local prediction. We use the RMSE around this area to make the interval
library(lmtest)
library(sandwich)

coeftest(food.ols, vcov = vcovHC(food.ols, "HC1")) 


head(unc.df)
unc.df.redux <- unc.df[sample(1:nrow(unc.df),150),]

library(nlme)
#m01 <- gls(wow~poly(wav,3), data=mp)#, correlation = corARMA(p=1))

temp.lai.fit <- gls(pred_LAI~LAI,data=unc.df.redux  )#, correlation = corARMA(p=1))
temp.lai.fit.within <- predict(temp.lai.fit)

V <- vcov(temp.lai.fit)
X <- model.matrix(pred_LAI~LAI,data=unc.df.redux  )
se.fit <- sqrt(diag(X %*% V %*% t(X)))

lwr=temp.lai.fit.within-1.96*se.fit
upr=temp.lai.fit.within+1.96*se.fit

lwr

ggplot(data = unc.df.redux  , aes(x = LAI, y = pred_LAI, ymin = lwr, ymax = upr)) + 
  geom_point() + 
  geom_ribbon(alpha = .3)    # set opacity so points are visible


#generate 30 samples of the CAB, test the validation agains the prediction and use those CI.




plot(temp.lai.fit.ci)
?plot.predict
temp.lai.fit.ci <- intervals(temp.lai.fit)

plot(temp.lai.fit)
summary(temp.lai.fit)
head(temp.lai.fit.ci$coef )
bb <- temp.lai.fit.ci$coef 

temp.lm <-lm(pred_Cab~Cab,data=unc.df)
temp.lm <-lm(pred_Cab~Cab,data=unc.df)

test.data <- c(1,5,4,10,9)
mean(test.data)
sd(test.data)
std.error(test.data)



newdata=data.frame(Cab=c(-10,20,30,40,50,70,100,200,300,1000))
predict(temp.lm, 
        newdata, 
        interval="confidence") 

sum.temp.lm <- summary(temp.lm)
sum.temp.lm$residuals

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
        interval="confidence") 

predict(test.lm, 
        newdata, 
        interval="prediction") 

predict(test.lm, 
        newdata)

predict(lm.Cab, 
        newdata, 
        interval="confidence")





head(mv.pred.mRF.unc)

unc.df <- cbind(param.list.unc,as.data.frame(mRF.s2.spc.sel.unc))
#we should also rename the bands to something proper
names(unc.df) <- c(names(unc.df)[1:5],
                   "B02","B03","B04",
                   "B05","B06","B07",
                   "B8A","B11","B12")

pred.mRF.unc <- predict(full.mRF,
                        newdata=unc.df)

mv.pred.mRF.unc <- as.data.frame(get.mv.predicted(pred.mRF.unc))


param.list <- data.frame(Cab=LHS[,1],
                         Car=LHS[,2],
                         Cw=LHS[,3],
                         Cm=LHS[,4],
                         LAI=LHS[,5])






###abandoned or testing code


data.frame(mv.pred.mRF.unc[1,1])
mv.pred.mRF.unc[1,1]

summary(test.lm)

plot(lm(mv.pred.mRF.unc[,1]~param.list.unc[,1]))

par(mfrow=c(1,1))
library(ggplot2)
x=param.list.unc[,5]
y=mv.pred.mRF.unc[,5]

ggplot(data.frame(x,y), aes(x=x,y=y)) + 
  geom_point() + 
  geom_smooth(method="lm",color='#2C3E50')



attach(stackloss)    # attach the data frame 
stackloss.lm = lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.)
stackloss.lm
newdata = data.frame(Air.Flow=72, Water.Temp=20, Acid.Conc.=85)
predict(stackloss.lm, newdata, interval="confidence") 
cc <- (stack.loss)

confint(test.lm)
