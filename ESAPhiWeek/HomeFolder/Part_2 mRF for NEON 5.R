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

param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
                         10,60, #Cab
                         5,25, #Car
                         0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.01,0.02, #Cm
                         0.5,6),#LAI
                         #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)

#creating a training space
train.n <- 500 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
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

#checking out the training and validation spectral spaces
par(mfrow=c(1,2))
plot(train.spclib)
plot(valid.spclib)

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
#if NA's appear they should be cleaned.. seems PROSAIL produces NA in extremes


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

#now we bring it into a DF
names(train.trait.df)
s2.mRF.train.df <- cbind(train.trait.df[,-c(6:8)],s2.train.spectr.20m.df) #we remove the fixed parameters
s2.mRF.valid.df <- cbind(valid.trait.df[,-c(6:8)],s2.valid.spectr.20m.df) #we remove the fixed parameters

#case of 7 params
names(s2.mRF.train.df )
#N.train.df   <- s2.mRF.train.df[,-c(2:7)]
Cab.train.df <- s2.mRF.train.df[,-c(2:5)]
Car.train.df <- s2.mRF.train.df[,-c(1,3:5)]
Cw.train.df  <- s2.mRF.train.df[,-c(1:2,4:5)]
Cm.train.df  <- s2.mRF.train.df[,-c(1:3,5)]
LAI.train.df <- s2.mRF.train.df[,-c(1:4)]
#hspot.train.df <- s2.mRF.train.df[,-c(1:6)]

#verify this
#names(N.train.df)
names(Cab.train.df)
names(Car.train.df)
names(Cw.train.df)
names(Cm.train.df)
names(LAI.train.df)
#names(hspot.train.df)

#training the single target regression model
#st.mRF.N    <- rfsrc(N~.  ,data=N.train.df,block.size = 50)
st.mRF.Cab  <- rfsrc(Cab~.,data=Cab.train.df,block.size = 50)
st.mRF.Car  <- rfsrc(Car~.,data=Car.train.df,block.size = 50)
st.mRF.Cw   <- rfsrc(Cw~. ,data=Cw.train.df,block.size = 50)
st.mRF.Cm   <- rfsrc(Cm~. ,data=Cm.train.df,block.size = 50)
st.mRF.LAI  <- rfsrc(LAI~.,data=LAI.train.df,block.size = 50)
#st.mRF.hsp  <- rfsrc(hspot~.,data=hspot.train.df,block.size = 50)

#creating predictions
st.test.df <- valid.trait.df[,-c(6:8)]
names(st.test.df)


#st.test.df$stPred_N   <- get.mv.predicted(predict(st.mRF.N,newdata=  s2.valid.spectr.20m.df))
st.test.df$stPred_Cab <- get.mv.predicted(predict(st.mRF.Cab,newdata=s2.valid.spectr.20m.df))
st.test.df$stPred_Car <- get.mv.predicted(predict(st.mRF.Car,newdata=s2.valid.spectr.20m.df))
st.test.df$stPred_Cw  <- get.mv.predicted(predict(st.mRF.Cw,newdata= s2.valid.spectr.20m.df))
st.test.df$stPred_Cm  <- get.mv.predicted(predict(st.mRF.Cm,newdata= s2.valid.spectr.20m.df))
st.test.df$stPred_LAI <- get.mv.predicted(predict(st.mRF.LAI,newdata=s2.valid.spectr.20m.df))
#st.test.df$stPred_hsp <- get.mv.predicted(predict(st.mRF.hsp,newdata=s2.valid.spectr.20m.df))

#lets plot these
names(valid.trait.df.clean)
names(valid.trait.df)
temp.df <- data.frame(trait=names(valid.trait.df[,-c(6:8)]),
                      slope=NA,inter=NA,Rsqur=NA,
                      RMSE=NA,MAE=NA,MAPE=NA)
temp.df
#separting the prediction from the table is actually useful now..
#names(st.test.df)
st.test.pred.df <- st.test.df[,6:ncol(st.test.df)]

rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }

head(valid.trait.df)
head(st.test.pred.df)

par(mfrow=c(2,3))
for (i in 1:nrow(temp.df)){
  #linear model stats
  x <- s2.mRF.valid.df[,i]
  y <- st.test.pred.df[,i]
  
  plot(x, y,
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  
  
  trait.fit <- lm(y~x)
  abline(trait.fit,lty=2,lwd=2)
  
  temp.df$slope[i] <- trait.fit$coefficients[[2]]
  temp.df$inter[i] <- trait.fit$coefficients[[1]]
  temp.df$Rsqur[i] <- summary(trait.fit)$r.square
  temp.df$RMSE[i]  <- RMSE.custom(trait.fit$residuals)
  temp.df$MAE[i]   <- mae.custom(trait.fit$residuals)
  temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

temp.df$Mean <- colMeans(s2.mRF.valid.df[,1:5])
temp.df$std  <- colSds(as.matrix(s2.mRF.valid.df[,1:5]))
temp.df

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
mtext("Single target models",cex=1.5)

########################################################
############# Multi output Random Forest################
########################################################
#MULTIOUTPUT Random Forest
names(s2.mRF.train.df)
names(s2.mRF.train.df)
dim(s2.mRF.train.df)

subsample <- s2.mRF.train.df[sample(1:nrow(s2.mRF.train.df),1000),]

head(s2.mRF.train.df)

#train time
mRF.Block50 <- rfsrc(Multivar(Cab,Car,
                              Cw,Cm,LAI)~.,data=s2.mRF.train.df,block.size = 50)
                     #nsplit=0,
                     #samptype = "swor",
                     #nodesize=1,
                     #mtry=8)
                     #ensemble='oob',
                     #ntree=5000,
                     #importance="none")


mRF.test.pred.df <- as.data.frame(get.mv.predicted(predict(mRF.Block50,
                                                           newdata=s2.valid.spectr.20m.df)))


#and a df for the fit summary
mRF.temp.df <- data.frame(trait=names(valid.trait.df[,-c(6:8)]),
                          slope=NA,inter=NA,Rsqur=NA,RMSE=NA,MAE=NA,MAPE=NA)

par(mfrow=c(2,3))
for (i in 1:nrow(temp.df)){
  
  #linear model stats
  x <- s2.mRF.valid.df[,i]
  y <- mRF.test.pred.df[,i]
  
  plot(y,x,
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=temp.df[i,1],pch=19,cex=.5,col="salmon")
  
  mRF.trait.fit <- lm(y~x)
  abline(mRF.trait.fit,lty=2,lwd=2)
  
  mRF.temp.df$slope[i] <- mRF.trait.fit$coefficients[[2]]
  mRF.temp.df$inter[i] <- mRF.trait.fit$coefficients[[1]]
  mRF.temp.df$Rsqur[i] <- summary(mRF.trait.fit)$r.square
  mRF.temp.df$RMSE[i]  <- RMSE.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAE[i]   <- mae.custom(mRF.trait.fit$residuals)
  mRF.temp.df$MAPE[i]  <- mean(abs((y-x)/x)*100)
}

#these are the same for all datasets
mRF.temp.df$Mean <- colMeans(s2.mRF.valid.df[,1:5])
mRF.temp.df$std  <- colSds(as.matrix(s2.mRF.valid.df[,1:5]))
mRF.temp.df

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
mtext("Multi target models",cex=1.5)


plot(nn)
head(nn.pred)

########################################################
############# Multi Task learning ######################
########################################################
#attempting a transformation
s2.mRF.train.df.t <- s2.mRF.train.df
s2.mRF.train.df.t$Cab <- exp(-s2.mRF.train.df.t$Cab/100)
s2.mRF.train.df.t$Car <- exp(-s2.mRF.train.df.t$Car/100)
s2.mRF.train.df.t$Cw <- exp(-s2.mRF.train.df.t$Cw/50)
s2.mRF.train.df.t$Cm <- exp(-s2.mRF.train.df.t$Cm/100)
s2.mRF.train.df.t$LAI <- exp(-s2.mRF.train.df.t$LAI/2)
head(s2.mRF.train.df.t)
s2.mRF.valid.df.t <- s2.mRF.valid.df
s2.mRF.valid.df.t$Cab <- exp(-s2.mRF.valid.df.t$Cab/100)
s2.mRF.valid.df.t$Car <- exp(-s2.mRF.valid.df.t$Car/100)
s2.mRF.valid.df.t$Cw <- exp(-s2.mRF.valid.df.t$Cw/50)
s2.mRF.valid.df.t$Cm <- exp(-s2.mRF.valid.df.t$Cm/100)
s2.mRF.valid.df.t$LAI <- exp(-s2.mRF.valid.df.t$LAI/2)
head(s2.mRF.valid.df.t)

#normalizing the variables
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

s2.mRF.train.df.n <- s2.mRF.train.df
s2.mRF.train.df.n$Cab <- normalize(s2.mRF.train.df.n$Cab)
s2.mRF.train.df.n$Car <- normalize(s2.mRF.train.df.n$Car)
s2.mRF.train.df.n$Cw <- normalize(s2.mRF.train.df.n$Cw)
s2.mRF.train.df.n$Cm <- normalize(s2.mRF.train.df.n$Cm)
s2.mRF.train.df.n$LAI <- normalize(s2.mRF.train.df.n$LAI)
head(s2.mRF.train.df.t)
s2.mRF.valid.df.n <- s2.mRF.valid.df
s2.mRF.valid.df.n$Cab <- normalize(s2.mRF.valid.df.n$Cab)
s2.mRF.valid.df.n$Car <- normalize(s2.mRF.valid.df.n$Car)
s2.mRF.valid.df.n$Cw <-  normalize(s2.mRF.valid.df.n$Cw)
s2.mRF.valid.df.n$Cm <-  normalize(s2.mRF.valid.df.n$Cm)
s2.mRF.valid.df.n$LAI <- normalize(s2.mRF.valid.df.n$LAI)
head(s2.mRF.valid.df.t)









##########################################
#checking neuralnet ######################
###########################################

#very intuitive and interesting

library(neuralnet)

sub.train.df <- s2.mRF.train.df.n[sample(1:nrow(s2.mRF.train.df.n),500),]
sub.valid.df <- s2.mRF.valid.df.n[sample(1:nrow(s2.mRF.valid.df.n),1000),]

names(sub.train.df)
head(sub.train.df)

softplus.fun <- function(x) log(1+exp(x))
lin.fun <- function(x) (x)
reLU.fun <- function(x){
  if (x < 0) {0
    }else (x)
  
}
reLU.fun(-1)
reLU.fun(2)
library(sigmoid)

nn=neuralnet(Cab+Car+Cw+Cm+LAI~.,data=sub.train.df ,
             algorithm='rprop+',
             hidden=c(9,5),act.fct = "tanh",
             err.fct = 'sse',
             linear.output = F)
nn.pred <- predict(nn, newdata=sub.valid.df, rep = 1, all.units = FALSE)


par(mfrow=c(2,3))
plot(sub.valid.df[,1],nn.pred[,1])
plot(sub.valid.df[,2],nn.pred[,2])
plot(sub.valid.df[,3],nn.pred[,3])
plot(sub.valid.df[,4],nn.pred[,4])
plot(sub.valid.df[,5],nn.pred[,5])

nn$err.fct()
plot(nn)
head(nn.pred)

library(ANN2)
neuralnetwork(X, y, hidden.layers, regression = FALSE,
              standardize = TRUE, loss.type = "log", huber.delta = 1,
              activ.functions = "tanh", step.H = 5, step.k = 100,
              optim.type = "sgd", learn.rates = 1e-04, L1 = 0, L2 = 0,
              sgd.momentum = 0.9, rmsprop.decay = 0.9, adam.beta1 = 0.9,
              adam.beta2 = 0.999, n.epochs = 100, batch.size = 32,
              drop.last = TRUE, val.prop = 0.1, verbose = TRUE,
              random.seed = NULL)

sub.train.df <- s2.mRF.train.df[sample(1:nrow(s2.mRF.train.df),500),]
sub.valid.df <- s2.mRF.valid.df[sample(1:nrow(s2.mRF.valid.df),1000),]

sub.train.df <- s2.mRF.train.df.t[sample(1:nrow(s2.mRF.train.df.t),500),]
sub.valid.df <- s2.mRF.valid.df.t[sample(1:nrow(s2.mRF.valid.df.t),1000),]
head(sub.train.df)

names(sub.train.df)
Y.mat <- as.matrix(sub.train.df[,c(1:5)])
X.mat <- as.matrix(sub.train.df[,c(6:14)])
X.mat.test <- as.matrix(sub.valid.df[,c(6:14)])
library(ANN2)
nn2 <- neuralnetwork(X=X.mat,y=Y.mat,regression=T,hidden.layers =c(10,6),loss.type = "squared",
                     activ.functions = c("relu",'tanh'),n.epochs = 10000,standardize = T)

nn2.pred <- predict(nn2,X.mat.test)
#head(nn2.pred)

par(mfrow=c(2,3))
plot(sub.valid.df[,1],nn2.pred$predictions[,1])
plot(sub.valid.df[,2],nn2.pred$predictions[,2])
plot(sub.valid.df[,3],nn2.pred$predictions[,3])
plot(sub.valid.df[,4],nn2.pred$predictions[,4])
plot(sub.valid.df[,5],nn2.pred$predictions[,5])

print(nn2)


plot(nn2)
o <- tune(Multivar(N,Cab,Car,
                   Cw,Cm,LAI,hspot) ~ ., subsample )
par(mfrow=c(1,1))
library(akima)
if (library("akima", logical.return = TRUE)) {
  ## nice little wrapper for plotting results
  plot.tune <- function(o, linear = TRUE) {
    x <- o$results[,1]
    y <- o$results[,2]
    z <- o$results[,3]
    so <- interp(x=x, y=y, z=z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]
    y0 <- y[idx]
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette =
                     colorRampPalette(c("yellow", "red")),
                   xlab = "nodesize",
                   ylab = "mtry",
                   main = "OOB error for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                     points(x,y,pch=16,cex=.25)})
  }
  ## plot the surface
  
  plot.tune(o)
}



library(RMTL)
names(s2.mRF.train.df)
X <- as.matrix(as.matrix(s2.mRF.train.df[,8:16]))
Y <- as.matrix(s2.mRF.train.df[,c(1:7)])

MTL.L21 <- MTL(X, Y, type = "Classification", Regularization = "L21",
               Lam1 = 0.1, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, 
                                                                  tol = 10^-3,
                                                                  maxIter = 1000), 
               G = NULL, k = 2)

data<-Create_simulated_data(Regularization="L21", type="Regression")

library(MultivariateRandomForest)
names(s2.mRF.train.df)
help(build_forest_predict)

row.ind <- sample(1:nrow(s2.mRF.valid.df),1000)
test.table <- s2.mRF.valid.df[row.ind ,]
head(test.table)
dim(as.matrix(s2.mRF.train.df[,c(1:7)]))
dim(as.matrix(s2.mRF.train.df[,8:16]))
dim(as.matrix(test.table[,8:16]))

mRF.new <- build_forest_predict(trainX=as.matrix(s2.mRF.train.df[,8:16]),
                                trainY=as.matrix(s2.mRF.train.df[,c(1:7)]),
                                n_tree=1000,
                                m_feature =5,
                                min_leaf=50,
                                testX = as.matrix(test.table))


model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 100, activation = 'tanh') %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units = 1, activation = 'sigmoid')





