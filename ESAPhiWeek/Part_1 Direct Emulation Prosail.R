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

#setwd
setwd("D:/ESAPhiWeek/")

dump.fld <- "./dumps"
dir.create(dump.fld) #it gives out a warning if the folder exits

#first we create 2 sets of variables, one for training and one for validation

#the limits of prosail params come from here
param.maxmin <- matrix(c(5,80, #Cab
                         1,25, #Car
                         0.005,0.02, #Cw
                         0.005,0.02, #Cm
                         0.1,8), #LAI
                       nrow=5,ncol = 2,byrow = T)






#creating a training space
train.n <- 500 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
train.LHS <- Latinhyper(param.maxmin,train.n)

valid.n <- 1500 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin,valid.n)

#thre is no need to linearize or normalize before since the LHS is generated independenly for each SET of parameters
#an alternative would be a non random generalization of the parameters
#e.g. on a grid -> this would naturally reduce any randomness thrown into the training and instead ensure that the 
# model is capturing all the different variations 
#but would likely strongly 
stable.grid <- Grid(param.maxmin,train.n)

#Now we can go into generating the data
stable.param.table <- data.frame(Cab=stable.grid [,1],
                                Car=stable.grid [,2],
                                Cw=stable.grid [,3],
                                Cm=stable.grid [,4],
                                LAI=stable.grid [,5])

train.param.table <- data.frame(Cab=train.LHS[,1],
                                Car=train.LHS[,2],
                                Cw=train.LHS[,3],
                                Cm=train.LHS[,4],
                                LAI=train.LHS[,5])

valid.param.table <- data.frame(Cab=valid.LHS[,1],
                                Car=valid.LHS[,2],
                                Cw=valid.LHS[,3],
                                Cm=valid.LHS[,4],
                                LAI=valid.LHS[,5])

#lets generate prosail outputs
stable.spclib <- PROSAIL(parameterList = stable.param.table)
stable.spectr <- spectra(stable.spclib)

train.spclib <- PROSAIL(parameterList = train.param.table)
train.spectr <- spectra(train.spclib)

valid.spclib <- PROSAIL(parameterList = valid.param.table)
valid.spectr <- spectra(valid.spclib)

#lets look at our inputs
par(mfrow=c(1,3))
plot(stable.spclib,main="Grid generated")
plot(train.spclib,main="Training LHS generated")
plot(valid.spclib,main="Validation LHS generated")
par(mfrow=c(1,1))
#some of these go under the limit of 0 on the reflectance


#first, lets generate a RF algorithm to predict spectra from traits (this as substitution from the
#prosail generator)

#first we decompose the spectras into components
m.stable.spectr <- as.matrix(stable.spectr) 
m.train.spectr <- as.matrix(train.spectr) 
m.valid.spectr <- as.matrix(valid.spectr) 

#notice that this is an SVD decomposition
pca.stable.spectr <- prcomp(m.stable.spectr)
pca.train.spectr <- prcomp(m.train.spectr)
pca.valid.spectr <- prcomp(m.valid.spectr)

pca.stable.cmeans <- colMeans(m.stable.spectr)
pca.train.cmeans <- colMeans(m.train.spectr)
pca.valid.cmeans <- colMeans(m.valid.spectr)


#lets check the comulative variance
pca.stable.vars <- apply(pca.stable.spectr$x,2,var)
pca.stable.prop <- pca.stable.vars/sum(pca.stable.vars)

pca.train.vars <- apply(pca.train.spectr$x,2,var)
pca.train.prop <- pca.train.vars/sum(pca.train.vars)

pca.valid.vars <- apply(pca.valid.spectr$x,2,var)
pca.valid.prop <- pca.valid.vars/sum(pca.valid.vars)

cumsum(pca.stable.prop)[1:10]
cumsum(pca.train.prop)[1:10]
cumsum(pca.valid.prop)[1:10]

#lets stay at .995 expained variance
nComp = 5

pca.stable.df <- cbind(pca.stable.spectr$x[,1:nComp],stable.param.table)
pca.train.df <- cbind(pca.train.spectr$x[,1:nComp],train.param.table)
pca.valid.df <- cbind(pca.valid.spectr$x[,1:nComp],valid.param.table)

head(pca.stable.df)

#lets now build the multi output RF based on these PC
mRF_pca.stable <- rfsrc(Multivar(PC1,PC2,PC3,PC4,PC5)~.,data=pca.stable.df,block.size = 1)
mRF_pca.train <- rfsrc(Multivar(PC1,PC2,PC3,PC4,PC5)~.,data=pca.train.df,block.size = 1)

#and predict
pred.mRF_pca.stable <- predict(mRF_pca.stable,newdata=valid.param.table)
pred.mRF_pca.train <- predict(mRF_pca.train,newdata=valid.param.table)

out.mRF_pca.stable <- get.mv.predicted(pred.mRF_pca.stable)
out.mRF_pca.train <- get.mv.predicted(pred.mRF_pca.train)

#each column is a predicted part of the output, we can now reconstruct the output
spec.pca.stable.pred <-  scale(out.mRF_pca.stable[,1:nComp] %*% t(pca.stable.spectr$rotation[,1:nComp]), 
                               center = -pca.stable.cmeans ,
                               scale = FALSE)

spec.pca.train.pred <-  scale(out.mRF_pca.train[,1:nComp] %*% t(pca.train.spectr$rotation[,1:nComp]), 
                               center = -pca.train.cmeans ,
                               scale = FALSE)

dim(df.param.list.test.spectra)
par(mfrow=c(1,1))
plot(m.valid.spectr [1,],type="l",ylab="Reflectance",xlab="Band index",main="Visual on simulation 1")
lines(spec.pca.stable.pred[1,],col=2,lty=2)
lines(spec.pca.train.pred[1,],col=3,lty=3)
legend("topright",c("Original","mRF - Grid training","mRF - LHS training"),col=c(1,2,3),lty=1:3, cex=0.8)

#lets compare the total sum and mean Spectral angle difference
#first we must convert them to a speclib
full.wavelength <- wavelength(stable.spclib)
pred.stable.speclib <- speclib(spec.pca.stable.pred,full.wavelength )
pred.train.speclib <- speclib(spec.pca.train.pred,full.wavelength )

#Spectral differences

#we will calculate the SAM differeces between each run of the model 

sam.pred2valid.df <- data.frame(RunNr=seq(1:valid.n),mRF_Grid=NA,mRF_LHS=NA)
for(i in 1:nrow(sam.pred2valid.df)){
  
  sam.pred2valid.df$mRF_Grid[i] <- sam(pred.stable.speclib[i,],valid.spclib[i,])
  sam.pred2valid.df$mRF_LHS[i] <- sam(pred.train.speclib[i,],valid.spclib[i,])
}


summary(sam.pred2valid.df )
write.csv(sam.pred2valid.df,
          paste(dump.fld,
                "FullTable_mRF_Direct_error_Hyperspectral.csv",sep="/"))
write.csv(summary(sam.pred2valid.df ),
          paste(dump.fld,
                "Summary_mRF_Direct_error_Hyperspectral.csv",sep="/"))

#we can now convert everything to a sentinel 2 data and see the band by band correlation

#the validatoin
s2.valid.speclib <- spectralResampling(valid.spclib,
                                             "MODIS",response_function = F)

#the predictions
s2.pred.stable.speclib <- spectralResampling(pred.stable.speclib,
                                             "MODIS",response_function = F)
s2.pred.train.speclib <- spectralResampling(pred.train.speclib,
                                            "MODIS",response_function = F)

#lets use only the 20mbands

s2.valid.speclib.sel <- s2.valid.speclib[,c(2,3,4,5,6,7,8,9,12,13)]

s2.pred.stable.speclib <- s2.pred.stable.speclib[,c(2,3,4,5,6,7,8,9,12,13)]
s2.pred.train.speclib <- s2.pred.train.speclib[,c(2,3,4,5,6,7,8,9,12,13)]

band.names <-  c("B02","B03","B04",
                 "B05","B06","B07",
                 "B08","B8A","B11",
                 "B12")
band.names.modis <-  paste("Band_",seq(1,19),sep="")

s2.valid.spectra <- spectra(s2.valid.speclib)
s2.stable.spectra <- spectra(s2.pred.stable.speclib)
s2.train.specta <- spectra(s2.pred.train.speclib)
s2.valid.spectra <- spectra(s2.valid.speclib)
s2.stable.spectra <- spectra(s2.pred.stable.speclib)
s2.train.specta <- spectra(s2.pred.train.speclib)



dim(s2.valid.spectra)
par(mfrow=c(2,3))
for (k in 1:ncol(s2.valid.spectra)){
  print(k)
  plot(s2.valid.spectra[,k],s2.train.specta[,k],
       #xlim=c(0,1),ylim=c(0,1),
       xlab="Reference",ylab="Predicted",
       main=band.names.modis [k],pch=19)
       #main=band.names[k],pch=19)
  #points(s2.valid.spectra[,k],s2.train.specta[,k],
  #     main=band.names[k],pch=19,col="red")
}

plot(s2.valid.speclib)
s2.valid.speclib
dim(s2.valid.spectra)


s2.valid.speclib@wavelength
s2.pred.train.speclib@wavelength

plot(s2.valid.speclib[,1],s2.pred.stable.speclib[,1])



names(s2.valid.speclib)<- band.names
names(s2.valid.speclib)<- band.names
names(s2.valid.speclib)<- band.names

pred.stable.speclib

bb <- spectra(pred.train.speclib)
cc <- spectra(valid.spclib)

tt <- bb[1,]
tt.1 <- bb[2,]
tt.1 <- bb[,1]

sam(spectra(pred.train.speclib)[1,],spectra(valid.spclib)[1,])

sam(pred.train.speclib[1,],valid.spclib[1,])
sam(pred.train.speclib[7501,],valid.spclib[7500,])
diff.stable.speclib

dim(spec.pca.stable.pred)

#testing
data(spectral_data)
spectra <- spectra(spectral_data)
wavelength <- spectral_data$wavelength

sam(spec.pca.stable.pred[1,],m.valid.spectr [1,])

#running the pca
spec.pca <- prcomp(my.spectr.plist.matrix)
spec.means = colMeans(my.spectr.plist.matrix)

nComp = 5
spec.pred = spec.pca$x[,1:nComp] %*% t(spec.pca$rotation[,1:nComp])
spec.pred = scale(spec.pred, center = -spec.means , scale = FALSE)



