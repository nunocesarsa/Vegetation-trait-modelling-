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
#param.maxmin <- matrix(c(10,80, #Cab
#                         5,40, #Car
#                         0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
#                         0.01,0.02, #Cm
#                         0.5,7),#LAI
#                       #0.05,0.1), #hotstop
#                       nrow=5,ncol = 2,byrow = T)

param.maxmin <- matrix(c(5,100, #Cab
                         5,50, #Car
                         0.005,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.03, #Cm
                         0.5,9.5),#LAI
                       #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)



#now we create a set of initial parameters to be able to start prosail
LHS.3000 <- Latinhyper(param.maxmin,3000)
LHS.valid <- Latinhyper(param.maxmin,6000)


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

spclib.plist.3000 <-PROSAIL(parameterList = train.par.df.3000 )
spclib.plist.valid <- PROSAIL(parameterList = valid.par.df)

s2.spclib.3000 <- spectralResampling(spclib.plist.3000,"Sentinel2",response_function = T)
s2.spclib.valid <- spectralResampling(spclib.plist.valid,"Sentinel2",response_function = T)

s2.spclib.3000.20m <- s2.spclib.3000[,c(2,3,4,5,6,7,9,12,13)]
s2.spclib.valid.20m <- s2.spclib.valid[,c(2,3,4,5,6,7,9,12,13)]

s2.spectra.3000 <- as.data.frame(spectra(s2.spclib.3000.20m))
s2.spectra.valid <- as.data.frame(spectra(s2.spclib.valid.20m))

s2.band.names <- c("B02","B03","B04",
                   "B05","B06","B07",
                   "B8A","B11","B12")

names(s2.spectra.3000) <- s2.band.names
names(s2.spectra.valid) <- s2.band.names

df.3000 <-  cbind(train.par.df.3000,s2.spectra.3000)
df.valid <- cbind(valid.par.df,s2.spectra.valid)

#setting up for ANN
Y.mat.3000 <- as.matrix(df.3000[,c(1:5)])
X.mat.3000 <- as.matrix(df.3000[,c(6:ncol(df.3000))])

X.mat.valid <- as.matrix(df.valid[,c(6:ncol(df.valid))])

#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c("tanh","relu")
lrates <- 0.005
epochs <- 3000
optype <- "adam"

#training the models
#creating and predicting the random forest
mRF <- rfsrc(Multivar(Cab,Car,Cw,Cm,LAI)~.,
                  data=df.3000 ,block.size = 25)
pred.mRF <- predict(mRF,
                    newdata=df.valid)
mv.pred.mRF <- as.data.frame(get.mv.predicted(pred.mRF))

#creating and predicting the ANN
ann <- neuralnetwork(X=X.mat.3000,y=Y.mat.3000,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype
                     )

ann.pred <- predict(ann,X.mat.valid)

#let's plot
head(df.valid)

mRF.df <- as.data.frame(mv.pred.mRF)
ann.df <- as.data.frame(ann.pred)

names(mRF.df) <- paste("mRF_",names(df.valid[,1:5]),sep="")
names(ann.df) <- paste("ANN_",names(df.valid[,1:5]),sep="")

out.df <- cbind(df.valid[,1:5],mRF.df,ann.df)

head(out.df)


tiff(paste(dump.ml.fld,"S2_MacLearn_Newdata_improvedplot.tif",sep = "/"),
     units="px", width = 2048, height = 1024, res=124,
     compression = c("lzw"))

par(mfrow=c(2,6))

line.size <- 3

for (i in 1:length(names(df.valid[,1:5]))){
  print(i)
  
  #i<- 1
  #auxiliary
  xlab <- paste("Reference",names(df.valid[,1:5]))
  ylab <- paste("Predicted",names(df.valid[,1:5]))
  
  #we also add 2 lines for each regression
  mrf.lm <- lm(out.df[,i+5]~out.df[,i])
  ann.lm <- lm(out.df[,i+10]~out.df[,i])
  
  #plots an empty background to force boundaries
  plot(out.df[,i], out.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab=xlab[i],ylab=ylab[i],
       main="Random Forest",pch=19,cex=.8,col="white")
  
  points(out.df[,i], out.df[,i+5],
       #xlim=c(0,1),ylim=c(0,1),
       pch=19,cex=.2,col="blue3")
  
  abline(mrf.lm,lty=2,lwd=line.size )
  #abline(ann.lm,lty=2,lwd=2,col="red4")

  #plots an empty background to force boundaries
  plot(out.df[,i], out.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab=xlab[i],ylab=ylab[i],
       main="Artificial neural network",pch=19,cex=.8,col="white")
  
  points(out.df[,i], out.df[,i+10],
         #xlim=c(0,1),ylim=c(0,1),
         pch=19,cex=.2,col="red4")
  
  abline(ann.lm,lty=2,lwd=line.size )
}

dev.off()


##################horizontal plot

tiff(paste(dump.ml.fld,"S2_MacLearn_Newdata_improvedplot_horizontal.tif",sep = "/"),
     units="px", width = 2048, height = 1024, res=124,
     compression = c("lzw"))

par(mfrow=c(2,5))

line.size <- 3

#first plot cycle
for (i in 1:length(names(df.valid[,1:5]))){
  print(i)
  
  #i<- 1
  #auxiliary
  xlab <- paste("Reference",names(df.valid[,1:5]))
  ylab <- paste("Predicted",names(df.valid[,1:5]))
  
  #we also add 2 lines for each regression
  mrf.lm <- lm(out.df[,i+5]~out.df[,i])
  #ann.lm <- lm(out.df[,i+10]~out.df[,i])
  
  #plots an empty background to force boundaries
  plot(out.df[,i], out.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab=xlab[i],ylab=ylab[i],
       #main="Random Forest",
       pch=19,cex=.8,col="white")
  
  points(out.df[,i], out.df[,i+5],
         #xlim=c(0,1),ylim=c(0,1),
         pch=19,cex=.2,col="blue3")
  
  abline(mrf.lm,lty=2,lwd=line.size )
  #abline(ann.lm,lty=2,lwd=2,col="red4")
  
  
}

for (i in 1:length(names(df.valid[,1:5]))){
  print(i)
  
  #i<- 1
  #auxiliary
  xlab <- paste("Reference",names(df.valid[,1:5]))
  ylab <- paste("Predicted",names(df.valid[,1:5]))
  
  #we also add 2 lines for each regression
  mrf.lm <- lm(out.df[,i+5]~out.df[,i])
  ann.lm <- lm(out.df[,i+10]~out.df[,i])
  
  
  #plots an empty background to force boundaries
  plot(out.df[,i], out.df[,i],
       #xlim=c(0,1),ylim=c(0,1),
       xlab=xlab[i],ylab=ylab[i],
       #main="Artificial neural network",
       pch=19,cex=.8,col="white")
  
  points(out.df[,i], out.df[,i+10],
         #xlim=c(0,1),ylim=c(0,1),
         pch=19,cex=.2,col="red4")
  
  abline(ann.lm,lty=2,lwd=line.size )
}



dev.off()