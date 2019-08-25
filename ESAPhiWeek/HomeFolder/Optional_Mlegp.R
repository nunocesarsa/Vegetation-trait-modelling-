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


#
library(mlegp)
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
head(param.list.test)

#we accept all other prosail parameters as default and we convert out spectral response to S resolution
mRF.spclib <- PROSAIL(parameterList = param.list)
mRF.s2.spclib <- spectralResampling(mRF.spclib,
                                    "Sentinel2",response_function = TRUE)

mRF.spclib.test <- PROSAIL(parameterList = param.list.test)
mRF.s2.spclib.test <- spectralResampling(mRF.spclib.test,
                                         "Sentinel2",response_function = TRUE)
#mlep tutorial
x = -5:5
z1 = 10 - 5*x + rnorm(length(x))
z2 = 7 * sin(x) + rnorm(length(x))
fitMulti = mlegp(x, cbind(z1,z2))
plot(fitMulti)
x
z2
#my case
mRF.s2.spc.sel <- mRF.s2.spclib[,c(2,3,4,5,6,7,9,12,13)]
bands.s2 <- as.data.frame(mRF.s2.spc.sel)
head(mRF.s2.spc.sel)

fitMulti.gp <- mlegp(param.list[,1],bands.s2)

par(mfrow=c(1,1))
x <- seq(0,1,l=10)
y <- abs(sin(2*pi*x))^.8
plot(x, y)

lm_mod <- lm(y ~ x)
plot(x, y)
abline(a=lm_mod$coef[1], b=lm_mod$coef[2], col='red')

library(GauPro)
gp <- GauPro(x, y, parallel=FALSE)
gp

plot(x, y)
curve(gp$predict(x), add=T, col=2)


plot(x, y)
y
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x)+2*gp$predict(x, se=T)$se, add=T, col=4)
curve(gp$predict(x)-2*gp$predict(x, se=T)$se, add=T, col=4)

bands.s2.m <- as.matrix(bands.s2)
head(bands.s2)
gp <- GauPro(X=bands.s2.m, Z=LHS, N=500,5,parallel=FALSE)
gp <- GauPro(X=bands2.s2.m, Z=LHS, N=500,5,parallel=FALSE,verbose=0)


gp <- GauPro(X=bands2.s2.m, Z=LHS, N=500,5,parallel=FALSE,verbose=0)


x <- matrix(seq(0,1,length.out = n), ncol=1)
y <- sin(2*pi*x) + rnorm(n,0,1e-1)


t.LHS <- t(LHS)
t.bands2.s2.m <- t(bands.s2.m)

head(bands.s2.m )


library(kernlab)
cab.df <- cbind(Cab=LHS[,1],bands.s2)
head(cab.df)
gp.cab <- gausspr(Cab~., data=cab.df)#,kernel="anovadot")

test.spectra <- as.data.frame(mRF.s2.spclib.test)

test.cab <- predict(gp.cab,test.spectra)

plot(LHS.test[,1],test.cab)

cab.car.df <- cbind(Cab=LHS[,1],Car=LHS[,2],bands.s2)
head(cab.car.df)
gp.cab <- gausspr(Cab+Car~., data=cab.car.df)

test.cab <- predict(gp.cab,test.spectra)

#no java no packag
library(extraTrees)
library(rJava)
sessionInfo()

library(GPfit)

#this library expects that the parameters are between 0 and 1
scale_norm(x = c(-1, 4, 10, 182))
# lower bound extended beyond -1
# upper bound still range of data
scale_norm(x = c(-1, 4, 10, 182), range = c(-100, 100))

#its a 5 parameter estimate
param.maxmin.wide <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  0.1,120, #Cab
  0.01,25, #Car
  0.00005,0.09, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.00005,0.05, #Cm
  0.1,9.5), #LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

#its a 5 parameter estimate
param.maxmin.wide.norm <- t(apply(param.maxmin,1,scale_norm))
param.maxmin.wide.norm 

#can we invert them back?
un_norm <- function (xmin,xmax,p) p*(xmax-xmin)+xmin
un_norm(0.1,120,0.5)


#we run a latin hypercube on our normalized space
LHS.norm <- Latinhyper(param.maxmin.wide.norm,prosail.runs*5)

#we have to reconvert each column by the conversion matric above
head(LHS.norm)

LHS.unnorm <- LHS.norm
LHS.unnorm[,1] <- LHS.norm[,1]*(param.maxmin.wide[1,2]-param.maxmin.wide[1,1])+param.maxmin.wide[1,1]
LHS.unnorm[,2] <- LHS.norm[,2]*(param.maxmin.wide[2,2]-param.maxmin.wide[2,1])+param.maxmin.wide[2,1]
LHS.unnorm[,3] <- LHS.norm[,3]*(param.maxmin.wide[3,2]-param.maxmin.wide[3,1])+param.maxmin.wide[3,1]
LHS.unnorm[,4] <- LHS.norm[,4]*(param.maxmin.wide[4,2]-param.maxmin.wide[4,1])+param.maxmin.wide[4,1]
LHS.unnorm[,5] <- LHS.norm[,5]*(param.maxmin.wide[5,2]-param.maxmin.wide[5,1])+param.maxmin.wide[5,1]

head(LHS.unnorm)
summary(LHS.unnorm)

parameterList.gp <- data.frame(Cab=LHS.unnorm[,1],
                               Car=LHS.unnorm[,2],
                               Cw=LHS.unnorm[,3],
                               Cm=LHS.unnorm[,4],
                               LAI=LHS.unnorm[,5])

#given these parameters, lets run prosail
mgp.spclib <- PROSAIL(parameterList = parameterList.gp)
mgp.s2.spclib <- spectralResampling(mgp.spclib,
                                    "Sentinel2",response_function = TRUE)

mgp.s2.spc.sel <- mgp.s2.spclib[,c(2,3,4,5,6,7,9,12,13)]

gp.train.df <- cbind(LHS.norm,as.data.frame(mgp.s2.spc.sel ))

m.LHS.norm <- as.matrix(LHS.norm)
m.spectra  <- as.matrix(as.data.frame(mgp.s2.spc.sel))

head(m.spectra)

dim(gp.train.df)
dim(LHS.norm)
dim(mgp.s2.spc.sel)

mlegp(t(LHS), t(m.spectra ), constantMean = 1, nugget = NULL, nugget.known = 0,
      min.nugget = 0, param.names = NULL, gp.names = NULL,
      PC.UD = NULL, PC.num = NULL, PC.percent = NULL,
      simplex.ntries = 5, simplex.maxiter = 500, simplex.reltol = 1e-8,
      BFGS.maxiter = 500, BFGS.tol = 0.01, BFGS.h = 1e-10, seed = 0,
      verbose = 1, parallel = FALSE)

library(RMTL)
data<-Create_simulated_data(Regularization="Graph", type="Regression",t=5,9,500)
#train a MTL model
#cold-start
model<-MTL(data$X, data$Y, type="Regression", Regularization="L21",
           Lam1=0.1, Lam2=0, opts=list(init=0, tol=10^-6, maxIter=1500))
plot(model)
xx <- data$X
yy <- data$Y
head(m.spectra )
dim(m.spectra)
yy[[1]] <- 
#lets create a feature set

#on my case i have 5 variables and 9 features, multiple observations

#lets create an object with the MTL structure
rtm.data <- list()
rtm.data$X <- m.spectra
dim(m.LHS.norm)
dim(m.spectra)

rtm.model <- MTL(m.spectra,m.LHS.norm)

#neuralnet
head(gp.train.df)
names(gp.train.df)[1:5] <- names(param.list)
names(gp.train.df)
library(neuralnet)
neuralnet(formula, data, hidden = 1, threshold = 0.01,
         stepmax = 1e+05, rep = 1, startweights = NULL,
         learningrate.limit = NULL, learningrate.factor = list(minus = 0.5,
                                                               plus = 1.2), learningrate = NULL, lifesign = "none",
         lifesign.step = 1000, algorithm = "rprop+", err.fct = "sse",
         act.fct = "logistic", linear.output = TRUE, exclude = NULL,
         constant.weights = NULL, likelihood = FALSE)

nn_mdl <- neuralnet(as.formula(Cab+Car+Cw+Cm+LAI~.),hidden=c(5,3),data=gp.train.df)
nn_mdl
#we run a latin hypercube on our normalized space
LHS.norm.test <- Latinhyper(param.maxmin.wide.norm,prosail.runs*10)
LHS.unnorm.test <- LHS.norm.test
LHS.unnorm.test[,1] <- LHS.norm.test[,1]*(param.maxmin.wide[1,2]-param.maxmin.wide[1,1])+param.maxmin.wide[1,1]
LHS.unnorm.test[,2] <- LHS.norm.test[,2]*(param.maxmin.wide[2,2]-param.maxmin.wide[2,1])+param.maxmin.wide[2,1]
LHS.unnorm.test[,3] <- LHS.norm.test[,3]*(param.maxmin.wide[3,2]-param.maxmin.wide[3,1])+param.maxmin.wide[3,1]
LHS.unnorm.test[,4] <- LHS.norm.test[,4]*(param.maxmin.wide[4,2]-param.maxmin.wide[4,1])+param.maxmin.wide[4,1]
LHS.unnorm.test[,5] <- LHS.norm.test[,5]*(param.maxmin.wide[5,2]-param.maxmin.wide[5,1])+param.maxmin.wide[5,1]

parameterList.test <- data.frame(Cab=LHS.unnorm.test[,1],
                               Car=LHS.unnorm.test[,2],
                               Cw=LHS.unnorm.test[,3],
                               Cm=LHS.unnorm.test[,4],
                               LAI=LHS.unnorm.test[,5])

mgp.spclib.test <- PROSAIL(parameterList = parameterList.test)
mgp.s2.spclib.test <- spectralResampling(mgp.spclib.test,
                                    "Sentinel2",response_function = TRUE)

mgp.s2.spc.sel.test <- mgp.s2.spclib.test[,c(2,3,4,5,6,7,9,12,13)]
df.s2.spc.sel.test <- as.data.frame(mgp.s2.spc.sel.test)

nn.predict <- predict(nn_mdl,newdata=df.s2.spc.sel.test)
