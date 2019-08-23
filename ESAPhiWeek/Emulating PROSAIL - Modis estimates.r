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

#for some metrics on acc
library(Metrics)

#setwd
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"
dir.create(dump.fld) #it gives out a warning if the folder exits

#first we create 2 sets of variables, one for training and one for validation

#the limits of prosail params come from here
param.maxmin <- matrix(c(5,80, #Cab
                         1,25, #Car
                         0.005,0.02, #Cw
                         0.005,0.02, #Cm
                         0.5,8), #LAI
                       nrow=5,ncol = 2,byrow = T)

param.maxmin <- matrix(c(0.8,2.5, #N, leaf layers
                         5,80, #Cab
                         1,25, #Car
                         0.005,0.02, #Cw
                         0.005,0.02, #Cm
                         0.1,8,#LAI
                         0,90), #mean leaf angle distribution
                       nrow=7,ncol = 2,byrow = T)





#creating a training space
train.n <- 500 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
train.LHS <- Latinhyper(param.maxmin,train.n)

valid.n <- 1000 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin,valid.n)

#thre is no need to linearize or normalize before since the LHS is generated independenly for each SET of parameters
#an alternative would be a non random generalization of the parameters
#e.g. on a grid -> this would naturally reduce any randomness thrown into the training and instead ensure that the 
# model is capturing all the different variations 
#but would likely strongly 
#stable.grid <- Grid(param.maxmin,train.n)

#Now we can go into generating the data
#stable.param.table <- data.frame(Cab=stable.grid [,1],
#                                Car=stable.grid [,2],
#                                Cw=stable.grid [,3],
#                                Cm=stable.grid [,4],
#                                LAI=stable.grid [,5])


train.param.table <- data.frame(N=train.LHS[,1],
                                Cab=train.LHS[,2],
                                Car=train.LHS[,3],
                                Cw=train.LHS[,4],
                                Cm=train.LHS[,5],
                                LAI=train.LHS[,6],
                                lidfa=train.LHS[,7],
                                TypeLidf = 0 ,
                                tts = 0,
                                tto = 30)

valid.param.table <- data.frame(N=train.LHS[,1],
                                Cab=valid.LHS[,2],
                                Car=valid.LHS[,3],
                                Cw=valid.LHS[,4],
                                Cm=valid.LHS[,5],
                                LAI=valid.LHS[,6],
                                lidfa=valid.LHS[,7],
                                TypeLidf = 0 ,
                                tts = 0,
                                tto = 30)

#lets generate prosail outputs
#stable.spclib <- PROSAIL(parameterList = stable.param.table)
#stable.spectr <- spectra(stable.spclib)

train.spclib <- PROSAIL(parameterList = train.param.table)
train.spectr <- spectra(train.spclib)

valid.spclib <- PROSAIL(parameterList = valid.param.table)
valid.spectr <- spectra(valid.spclib)

#lets look at our inputs
par(mfrow=c(1,2))
#plot(stable.spclib,main="Grid generated")
plot(train.spclib,main="Training LHS generated")
plot(valid.spclib,main="Validation LHS generated")
par(mfrow=c(1,1))
#some of these go under the limit of 0 on the reflectance


#first, lets generate a RF algorithm to predict spectra from traits (this as substitution from the
#prosail generator)

#first we decompose the spectras into components
#m.stable.spectr <- as.matrix(stable.spectr) 
m.train.spectr <- as.matrix(train.spectr) 
m.valid.spectr <- as.matrix(valid.spectr) 

#notice that this is an SVD decomposition
#pca.stable.spectr <- prcomp(m.stable.spectr)
pca.train.spectr <- prcomp(m.train.spectr)
#pca.valid.spectr <- prcomp(m.valid.spectr)

#pca.stable.cmeans <- colMeans(m.stable.spectr)
pca.train.cmeans <- colMeans(m.train.spectr)
#pca.valid.cmeans <- colMeans(m.valid.spectr)


#lets check the comulative variance
#pca.stable.vars <- apply(pca.stable.spectr$x,2,var)
#pca.stable.prop <- pca.stable.vars/sum(pca.stable.vars)

pca.train.vars <- apply(pca.train.spectr$x,2,var)
pca.train.prop <- pca.train.vars/sum(pca.train.vars)

#pca.valid.vars <- apply(pca.valid.spectr$x,2,var)
#pca.valid.prop <- pca.valid.vars/sum(pca.valid.vars)

#cumsum(pca.stable.prop)[1:10]
cumsum(pca.train.prop)[1:10]
#cumsum(pca.valid.prop)[1:10]

#lets stay at .999 expained variance
nComp = 7

#pca.stable.df <- cbind(pca.stable.spectr$x[,1:nComp],stable.param.table)
pca.train.df <- cbind(pca.train.spectr$x[,1:nComp],train.param.table[,1:7])
pca.valid.df <- cbind(pca.valid.spectr$x[,1:nComp],valid.param.table[,1:7])

head(pca.train.df)

#lets now build the multi output RF based on these PC
#mRF_pca.stable <- rfsrc(Multivar(PC1,PC2,PC3,PC4,PC5)~.,data=pca.stable.df,block.size = 30)
mRF_pca.train <- rfsrc(Multivar(PC1,PC2,PC3,PC4,PC5,PC6,PC7)~.,data=pca.train.df,block.size = 30)

#and predict
#pred.mRF_pca.stable <- predict(mRF_pca.stable,newdata=valid.param.table)
head(valid.param.table)
pred.mRF_pca.train <- predict(mRF_pca.train,newdata=valid.param.table[,1:7])

#out.mRF_pca.stable <- get.mv.predicted(pred.mRF_pca.stable)
out.mRF_pca.train <- get.mv.predicted(pred.mRF_pca.train)

#each column is a predicted part of the output, we can now reconstruct the output
#spec.pca.stable.pred <-  scale(out.mRF_pca.stable[,1:nComp] %*% t(pca.stable.spectr$rotation[,1:nComp]), 
#                               center = -pca.stable.cmeans ,
#                               scale = FALSE)

spec.pca.train.pred <-  scale(out.mRF_pca.train[,1:nComp] %*% t(pca.train.spectr$rotation[,1:nComp]), 
                               center = -pca.train.cmeans ,
                               scale = FALSE)

dim(df.param.list.test.spectra)
par(mfrow=c(1,1))
plot(m.valid.spectr [1,],type="l",ylab="Reflectance",xlab="Band index",main="Visual on simulation 1")
#lines(spec.pca.stable.pred[1,],col=2,lty=2)
lines(spec.pca.train.pred[1,],col=2,lty=3)
legend("topright",c("Original","mRF - LHS training"),col=c(1,2),lty=1:2, cex=0.8)
#legend("topright",c("Original","mRF - Grid training","mRF - LHS training"),col=c(1,2,3),lty=1:3, cex=0.8)

#lets compare the total sum and mean Spectral angle difference
#first we must convert them to a speclib
full.wavelength <- wavelength(valid.spclib)
#pred.stable.speclib <- speclib(spec.pca.stable.pred,full.wavelength )
pred.train.speclib <- speclib(spec.pca.train.pred,full.wavelength )

par(mfrow=c(1,2))
plot(valid.spclib,main="Original")
plot(pred.train.speclib,main="Predicted space")

#Spectral differences

#we will calculate the SAM differeces between each run of the model 

sam.pred2valid.df <- data.frame(RunNr=seq(1:valid.n),mRF_Grid=NA,mRF_LHS=NA)
for(i in 1:nrow(sam.pred2valid.df)){
  
  #sam.pred2valid.df$mRF_Grid[i] <- sam(pred.stable.speclib[i,],valid.spclib[i,])
  sam.pred2valid.df$mRF_LHS[i] <- sam(pred.train.speclib[i,],valid.spclib[i,])
}


summary(sam.pred2valid.df )
write.csv(sam.pred2valid.df,
          paste(dump.fld,
                "FullTable_mRF_Direct_error_Hyperspectral_7params.csv",sep="/"))
write.csv(summary(sam.pred2valid.df ),
          paste(dump.fld,
                "Summary_mRF_Direct_error_Hyperspectral_7params.csv",sep="/"))

#we can now convert everything to a sentinel 2 data and see the band by band correlation

#the validatoin
s2.valid.speclib <- spectralResampling(valid.spclib,
                                             "MODIS",response_function = F)

#the predictions
#s2.pred.stable.speclib <- spectralResampling(pred.stable.speclib,
#                                             "MODIS",response_function = F)
s2.pred.train.speclib <- spectralResampling(pred.train.speclib,
                                            "MODIS",response_function = F)
s2.pred.train.speclib@wavelength

get.sensor.characteristics("MODIS", response_function = FALSE)
#lets use only the 20mbands
sentinel.band.sel <- c(2,3,4,5,6,7,8,9,12,13)
MODIS.band.sel <- c(1,2,3,4,5,6,7)

#if sentinel
#s2.valid.speclib.sel <- s2.valid.speclib[,sentinel.band.sel]
#s2.pred.stable.speclib <- s2.pred.stable.speclib[,sentinel.band.sel]
#s2.pred.train.speclib <- s2.pred.train.speclib[,sentinel.band.sel]

#if MODIS
s2.valid.speclib.sel <- s2.valid.speclib[,MODIS.band.sel]
#s2.pred.stable.speclib <- s2.pred.stable.speclib[,MODIS.band.sel]
s2.pred.train.speclib <- s2.pred.train.speclib[,MODIS.band.sel]


band.names <-  c("B02","B03","B04",
                 "B05","B06","B07",
                 "B08","B8A","B11",
                 "B12")
band.names.modis <-  paste("Band_",seq(1,7),sep="")



s2.valid.spectra <- spectra(s2.valid.speclib.sel)
#s2.stable.spectra <- spectra(s2.pred.stable.speclib)
s2.train.spectra <- spectra(s2.pred.train.speclib)
#lets remove outliers before plotting
for (k in 1:ncol(s2.valid.speclib.sel)){
  #print(k)
  #k <- 1
  outlier.vals <- boxplot(s2.valid.spectra[,k],plot=F)$out
  
  lin.index <- which(s2.valid.spectra[,k] %in% outlier.vals )
  #lin.index
  if (length(lin.index)>0){
    s2.valid.spectra <- s2.valid.spectra[-lin.index,]
    s2.train.spectra <- s2.train.spectra[-lin.index,]
    
  }
  
 
  
  
  #plot(s2.valid.spectra[,k],s2.train.specta[,k],
  #     #xlim=c(0,1),ylim=c(0,1),
  #     xlab="Reference",ylab="Predicted",
  #     main=band.names.modis[k],pch=19,cex=.5)
  #main=band.names[k],pch=19)
  #points(s2.valid.spectra[,k],s2.train.specta[,k],
  #     main=band.names[k],pch=19,col="red")
}

temp.df <- data.frame(band=band.names.modis,slope=NA,intercept=NA,R2=NA,RMSE=NA,MAE=NA)
rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }
par(mfrow=c(2,4))

for (k in 1:ncol(s2.valid.speclib.sel)){
  print(k)
  plot(s2.valid.spectra[,k],s2.train.spectra[,k],
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=band.names.modis[k],pch=19,cex=.5,col="lightsalmon")
  
       #main=band.names[k],pch=19)
  #points(s2.valid.spectra[,k],s2.train.specta[,k],
  #     main=band.names[k],pch=19,col="red")
  
  bb <- (s2.valid.spectra)[,k]
  cc <- (s2.valid.spectra)[,k]
  best.fit <- lm(bb~cc)
  band.fit <- lm(s2.train.spectra[,k]~cc)
  
  abline(best.fit,lty=2,lwd=1.5)
  abline(band.fit,lty=1,col="blue",lwd=1.5)
  
  #print(band.names.modis[k])
  #print(band.fit$coefficients)
  #print(rsq(cc,s2.train.spectra[,k]))
  #print(RMSE(s2.train.spectra[,k],obs = s2.valid.spectra[,k]))
  temp.df$slope[k] <- band.fit$coefficients[[2]]
  temp.df$intercept[k] <- band.fit$coefficients[[1]]
  temp.df$R2[k]   <- summary(band.fit)$r.square
  temp.df$RMSE[k] <- RMSE.custom(band.fit$residuals)
  temp.df$MAE[k] <- mae.custom(band.fit$residuals)
  
  #temp.df$RMSE[k] <- rmse(s2.valid.spectra[,k],s2.train.spectra[,k])
  #temp.df$MSE[k]  <- mse(s2.valid.spectra[,k],s2.train.spectra[,k])
}

temp.df
write.csv(temp.df,
          paste(dump.fld,
                "Summary_mRF_ModisSimul.csv",sep="/"))





##########################################3
#single objective modelling
##################################################


#pca.stable.df <- cbind(pca.stable.spectr$x[,1:nComp],stable.param.table)
pca.train.df <- cbind(pca.train.spectr$x[,1:nComp],train.param.table[,1:7])
pca.valid.df <- cbind(pca.valid.spectr$x[,1:nComp],valid.param.table[,1:7])

#we predict each component of the PCA in function of the multivariate trait space
head(pca.train.df)
pca.PC1.df  <- pca.train.df[,c(1,8:14)]
pca.PC2.df  <- pca.train.df[,c(2,8:14)]
pca.PC3.df  <- pca.train.df[,c(3,8:14)]
pca.PC4.df  <- pca.train.df[,c(4,8:14)]
pca.PC5.df  <- pca.train.df[,c(5,8:14)]
pca.PC6.df  <- pca.train.df[,c(6,8:14)]
pca.PC7.df  <- pca.train.df[,c(7,8:14)]

#training the single target model - by multivariate traits im giving it a looooooooooot of leverage
mRF_pca1 <- rfsrc(PC1~.,data=pca.PC1.df,block.size = 30)
mRF_pca2 <- rfsrc(PC2~.,data=pca.PC2.df,block.size = 30)
mRF_pca3 <- rfsrc(PC3~.,data=pca.PC3.df,block.size = 30)
mRF_pca4 <- rfsrc(PC4~.,data=pca.PC4.df,block.size = 30)
mRF_pca5 <- rfsrc(PC5~.,data=pca.PC5.df,block.size = 30)
mRF_pca6 <- rfsrc(PC6~.,data=pca.PC6.df,block.size = 30)
mRF_pca7 <- rfsrc(PC7~.,data=pca.PC7.df,block.size = 30)


#now we slowly recreate a multi target output from single target models
#pred.mRF_pca.stable <- predict(mRF_pca.stable,newdata=valid.param.table)
head(valid.param.table[,1:7])
pred.mRF_pca1 <- predict(mRF_pca1,newdata=valid.param.table[,1:7])
pred.mRF_pca2 <- predict(mRF_pca2,newdata=valid.param.table[,1:7])
pred.mRF_pca3 <- predict(mRF_pca3,newdata=valid.param.table[,1:7])
pred.mRF_pca4 <- predict(mRF_pca4,newdata=valid.param.table[,1:7])
pred.mRF_pca5 <- predict(mRF_pca5,newdata=valid.param.table[,1:7])
pred.mRF_pca6 <- predict(mRF_pca6,newdata=valid.param.table[,1:7])
pred.mRF_pca7 <- predict(mRF_pca7,newdata=valid.param.table[,1:7])

#getting the predictions
out.pred.mRF_pca1 <- get.mv.predicted(pred.mRF_pca1)
out.pred.mRF_pca2 <- get.mv.predicted(pred.mRF_pca2)
out.pred.mRF_pca3 <- get.mv.predicted(pred.mRF_pca3)
out.pred.mRF_pca4 <- get.mv.predicted(pred.mRF_pca4)
out.pred.mRF_pca5 <- get.mv.predicted(pred.mRF_pca5)
out.pred.mRF_pca6 <- get.mv.predicted(pred.mRF_pca6)
out.pred.mRF_pca7 <- get.mv.predicted(pred.mRF_pca7)

#and lets make the house once again
single.pred.PCA <- cbind(out.pred.mRF_pca1,out.pred.mRF_pca2,out.pred.mRF_pca3,out.pred.mRF_pca4,
                         out.pred.mRF_pca5,out.pred.mRF_pca6,out.pred.mRF_pca7) 

head(single.pred.PCA)
single.pred.spectra <-  scale(single.pred.PCA[,1:nComp] %*% t(pca.train.spectr$rotation[,1:nComp]), 
                              center = -pca.train.cmeans ,
                              scale = FALSE)

par(mfrow=c(1,1))
plot(m.valid.spectr [1,],type="l",ylab="Reflectance",xlab="Band index",main="Visual on simulation 1")
lines(single.pred.spectra[1,],col=2,lty=3)
legend("topright",c("Original","stRF - LHS training -"),col=c(1,2),lty=1:2, cex=0.8)
#legend("topright",c("Original","mRF - Grid training","mRF - LHS training"),col=c(1,2,3),lty=1:3, cex=0.8)

#full.wavelength <- wavelength(valid.spclib)
#pred.stable.speclib <- speclib(spec.pca.stable.pred,full.wavelength )
single.pred.speclib <- speclib(single.pred.spectra,full.wavelength )

par(mfrow=c(1,2))
plot(valid.spclib,main="Original")
plot(pred.train.speclib,main="Predicted space")


sam.single.pred2valid.df <- data.frame(RunNr=seq(1:valid.n),stRF_Grid=NA,stRF_LHS=NA)
for(i in 1:nrow(sam.single.pred2valid.df)){
  
  #sam.pred2valid.df$mRF_Grid[i] <- sam(pred.stable.speclib[i,],valid.spclib[i,])
  sam.single.pred2valid.df$stRF_LHS[i] <- sam(single.pred.speclib[i,],valid.spclib[i,])
}


summary(sam.single.pred2valid.df )
write.csv(sam.single.pred2valid.df ,
          paste(dump.fld,
                "FullTable_stRF_Direct_error_Hyperspectral_7params.csv",sep="/"))
write.csv(summary(sam.single.pred2valid.df ),
          paste(dump.fld,
                "Summary_stRF_Direct_error_Hyperspectral_7params.csv",sep="/"))

###################################################
### now we go into modelling MODIS responses - single target training
###################################################


#Ressample the validation into MODIS bands - notice its not a gaussian response but just an averaged response
s2.valid.speclib <- spectralResampling(valid.spclib,
                                       "MODIS",response_function = F)

#the predictions
#s2.pred.stable.speclib <- spectralResampling(pred.stable.speclib,
#                                             "MODIS",response_function = F)
s2.st.pred.speclib <- spectralResampling(single.pred.speclib,
                                            "MODIS",response_function = F)
s2.st.pred.speclib@wavelength

get.sensor.characteristics("MODIS", response_function = FALSE)
#lets use only the 20mbands
sentinel.band.sel <- c(2,3,4,5,6,7,8,9,12,13)
MODIS.band.sel <- c(1,2,3,4,5,6,7)

#if sentinel
#s2.valid.speclib.sel <- s2.valid.speclib[,sentinel.band.sel]
#s2.pred.stable.speclib <- s2.pred.stable.speclib[,sentinel.band.sel]
#s2.pred.train.speclib <- s2.pred.train.speclib[,sentinel.band.sel]

#if MODIS
s2.valid.speclib.sel <- s2.valid.speclib[,MODIS.band.sel]
#s2.pred.stable.speclib <- s2.pred.stable.speclib[,MODIS.band.sel]
s2.st.pred.speclib <- s2.st.pred.speclib[,MODIS.band.sel]


band.names <-  c("B02","B03","B04",
                 "B05","B06","B07",
                 "B08","B8A","B11",
                 "B12")
band.names.modis <-  paste("Band_",seq(1,7),sep="")



s2.valid.spectra <- spectra(s2.valid.speclib.sel)
#s2.stable.spectra <- spectra(s2.pred.stable.speclib)
s2.st.pred.spectra <- spectra(s2.st.pred.speclib)
#lets remove outliers before plotting
for (k in 1:ncol(s2.valid.speclib.sel)){
  #print(k)
  #k <- 1
  outlier.vals <- boxplot(s2.valid.spectra[,k],plot=F)$out
  
  lin.index <- which(s2.valid.spectra[,k] %in% outlier.vals )
  #lin.index
  if (length(lin.index)>0){
    s2.valid.spectra <- s2.valid.spectra[-lin.index,]
    s2.st.pred.spectra <- s2.st.pred.spectra[-lin.index,]
    
  }
  
  
  
  
  #plot(s2.valid.spectra[,k],s2.train.specta[,k],
  #     #xlim=c(0,1),ylim=c(0,1),
  #     xlab="Reference",ylab="Predicted",
  #     main=band.names.modis[k],pch=19,cex=.5)
  #main=band.names[k],pch=19)
  #points(s2.valid.spectra[,k],s2.train.specta[,k],
  #     main=band.names[k],pch=19,col="red")
}

temp.df <- data.frame(band=band.names.modis,slope=NA,intercept=NA,R2=NA,RMSE=NA,MAE=NA)
rsq <- function (x, y) cor(x, y) ^ 2
RMSE.custom <- function(error) { sqrt(mean(error^2)) }
mae.custom  <- function(error) { mean(abs(error)) }
par(mfrow=c(2,4))

for (k in 1:ncol(s2.valid.speclib.sel)){
  print(k)
  plot(s2.valid.spectra[,k], s2.st.pred.spectra[,k],
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=band.names.modis[k],pch=19,cex=.5,col="lightsalmon")
  
  #main=band.names[k],pch=19)
  #points(s2.valid.spectra[,k],s2.train.specta[,k],
  #     main=band.names[k],pch=19,col="red")
  
  bb <- (s2.valid.spectra)[,k]
  cc <- (s2.valid.spectra)[,k]
  best.fit <- lm(bb~cc)
  band.fit <- lm(s2.st.pred.spectra[,k]~cc)
  
  abline(best.fit,lty=2,lwd=1.5)
  abline(band.fit,lty=1,col="blue",lwd=1.5)
  
  #print(band.names.modis[k])
  #print(band.fit$coefficients)
  #print(rsq(cc,s2.train.spectra[,k]))
  #print(RMSE(s2.train.spectra[,k],obs = s2.valid.spectra[,k]))
  temp.df$slope[k] <- band.fit$coefficients[[2]]
  temp.df$intercept[k] <- band.fit$coefficients[[1]]
  temp.df$R2[k]   <- summary(band.fit)$r.square
  temp.df$RMSE[k] <- RMSE.custom(band.fit$residuals)
  temp.df$MAE[k] <- mae.custom(band.fit$residuals)
  
  #temp.df$RMSE[k] <- rmse(s2.valid.spectra[,k],s2.train.spectra[,k])
  #temp.df$MSE[k]  <- mse(s2.valid.spectra[,k],s2.train.spectra[,k])
}

temp.df
write.csv(temp.df,
          paste(dump.fld,
                "Summary_stRF_ModisSimul.csv",sep="/"))
