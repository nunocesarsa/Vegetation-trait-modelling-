#library for PROSAIL and also Spectral Angle Mapper

library(hsdar)

#library for optimization procedure - you still should read what this actually does lol
library(SoilHyP)

#library that has some extra random generation tools
library(MCMCglmm)

#This should help increase the CPU cores used
library(snowfall)
#library(parallel)
#sfParallel() #should respond 1 cpu
#sfCpus()
#detectCores()
#lets use always half the CPU core available
#half_cpu_cores <- round(detectCores()/2)
#sfInit(parallel = TRUE, cpus = half_cpu_cores, slaveOutfile = "slave.out",nostart = F)
#stopifnot( sfCpus() == 4 )
#stopifnot( sfParallel() == TRUE )
#sfStop()


#first step, creating a function to minimize
#we minimize the Spectral angle mapper - meaning, the difference between two generated spectra from PROSAIL model
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"
set.seed(1) #ITS IMPORTANT TO DO THIS - ensures you always have the same outputs.. but also implies you also have to run the entire scrippt unless you keep returing the seed to one later



#this function of SAM shouls compare agasint a given SPECTRAL LIBRARY
sam.fun.spclib <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  init.param <- data.frame(Cab=Cab_init,
                           Cw=Cw_init,
                           LAI=LAI_init)
  
  #this fetches new prosail parameters.. if you want to minimize towards an observation, then this, should be the observation
  
  if (is.null(target.spectra)==T){
    #print("Warning: its minimizing to the default PROSAIL parameters") #uncomment if you want to get a spammy console
    
    my.spectra <- PROSAIL(parameterList = init.param)
    #prosail default parameters 
    def.prosail <- PROSAIL() 
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,def.prosail)
    fin.val <- spec.dist[1]
  }
  
  if (is.null(target.spectra)==F) {
    
    #print("Attempting to find the vegetation traits") #uncomment if you want to get a spammy console
    #print(prosail.param)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,target.spectra)
    fin.val <- spec.dist[1]
    
  }
  
  
  return(fin.val)
}

#now we set up the upper and lower limits (important! to ensure convergence)
lower_lim=c(30,0.005,0.5)
upper_lim=c(60,0.03,4)

#first we create our initial set of parameters, this could be anything random, hopefully within the limits
init.par <- c(Cab=rnorm(1,mean=40,sd=15),
              Cw=rnorm(1,mean=0.02,sd=0.005),
              LAI=rtnorm(1,mean=1,sd=3,lower=0,upper=6)) #these are generated from random samples

init.par <- c(40,0.02,2) #same as before just simpler.. and not random

prosail.runs <- 20
df.param.list <- data.frame(Cab=rtnorm(prosail.runs,mean=40,sd=15,lower=0,upper=100),
                            Cw=rtnorm(prosail.runs,mean=0.02,sd=0.005,lower=0,upper=0.05),
                            LAI=rtnorm(prosail.runs,mean=3,sd=2,lower=0,upper=6))
head(df.param.list)
#these are the parameters we want to achieve in the end, but lets say we did not know them
#all we had was some prior expectation of what the values would be

#imagine from literature we knew this
init.param <- data.frame(Cab=40,
                         Cw=0.02,
                         LAI=2)
init.param <- c(40, #i haven't made a function to convert df to numeric list
                0.02,
                2)


#and we inferred a CI around the mean as this one:
lower_lim_list=c(0,0,0)
upper_lim_list=c(100,0.05,6)

#lets create a DF to receive the results
out.df <- df.param.list
out.df$SCE_Cab <- NA
out.df$SCE_Cw <- NA
out.df$SCE_LAI <- NA

out.df$LBFGSB_Cab <- NA
out.df$LBFGSB_Cw  <- NA
out.df$LBFGSB_LAI <- NA
out.df

for( i in 1:nrow(out.df)){
  print(paste("Estimating sample:", i))
  
  #print(df.param.list[i,])
  tgt.spec <- PROSAIL(parameterList = df.param.list[i,]) #notice, i have to generate the spectra but these could be your own spectra samples...
  
  
  
  SCEoptim.pred <- SCEoptim(FUN = sam.fun.spclib,
                            par = init.param,
                            #method="L-BFGS-B",
                            lower=lower_lim_list,
                            upper=upper_lim_list,
                            target.spectra=tgt.spec)
  
  LBFGSB.pred <- optim(init.param, sam.fun.spclib, gr = NULL,
                       target.spectra=tgt.spec,
                       method = c( "L-BFGS-B"),
                       lower = lower_lim_list, upper = upper_lim_list,
                       control = list(), hessian = FALSE)
  
  #storing the outputs
  out.df$SCE_Cab[i] <- SCEoptim.pred$par[1]
  out.df$SCE_Cw[i]  <- SCEoptim.pred$par[2]
  out.df$SCE_LAI[i] <- SCEoptim.pred$par[3]
  
  out.df$LBFGSB_Cab[i] <- LBFGSB.pred$par[1]
  out.df$LBFGSB_Cw[i]  <- LBFGSB.pred$par[2]
  out.df$LBFGSB_LAI[i] <- LBFGSB.pred$par[3]
  
}



#out.df.omit <- na.omit(out.df)

par(mfrow=c(1,4))
plot(out.df$Cab,out.df$Cab,pch=19,xlab="Reference",ylab="Prediction",main="Cab",cex=1.5)
points(out.df$Cab,out.df$SCE_Cab,pch=2,col="blue",cex=.75)
points(out.df$Cab,out.df$LBFGSB_Cab,pch=3,col="red",cex=.75)

plot(out.df$Cw,out.df$Cw,pch=19,xlab="Reference",ylab="Prediction",main="Cw",cex=1.5)
points(out.df$Cw,out.df$SCE_Cw,pch=2,col="blue",cex=.75)
points(out.df$Cw,out.df$LBFGSB_Cw,pch=3,col="red",cex=.75)

plot(out.df$LAI,out.df$SCE_LAI,pch=19,xlab="Reference",ylab="Prediction",main="LAI",cex=1.5)
points(out.df$LAI,out.df$SCE_LAI,pch=2,col="blue",cex=.75)
points(out.df$LAI,out.df$LBFGSB_LAI,pch=3,col="red",cex=.75)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','SCE','L-BFGS-B'), 
       pch=c(19,2,3), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red'))

write.csv(out.df,paste(dump.fld,
                       "Inverse_optimization_3Params_SCEvsLBFGSB.csv",sep="/"))



#we have trained by optimization, lets try the same but using now mRF approach
#first we generate a training space

#first we create our initial set of parameters, this could be anything random, hopefully within the limit
n.params <- 3

train.par.df.300 <- data.frame(Cab=rnorm(100*n.params,mean=40,sd=15),
                           Cw=rnorm(100*n.params,mean=0.02,sd=0.005),
                           LAI=rtnorm(100*n.params,mean=1,sd=3,lower=0,upper=6))
train.par.df.900 <- data.frame(Cab=rnorm(300*n.params,mean=40,sd=15),
                               Cw=rnorm(300*n.params,mean=0.02,sd=0.005),
                               LAI=rtnorm(300*n.params,mean=1,sd=3,lower=0,upper=6))
train.par.df.1500 <- data.frame(Cab=rnorm(500*n.params,mean=40,sd=15),
                               Cw=rnorm(500*n.params,mean=0.02,sd=0.005),
                               LAI=rtnorm(500*n.params,mean=1,sd=3,lower=0,upper=6))

#the previous uses a truncated normal distribution which might not necessarely include samples in the whole
#dimension of the solution space

#same as before but using a LHS approach (tried grid elsewhere and results were awkward)
df.param.list
summary(df.param.list)
#lets train given any param list
summary(df.param.list[,1])[[6]] 
#what i am doing here is broadening the space of the training so that i ensure some likelihood that the LHS will capture points that include all the extremes
min.Cab <- summary(df.param.list[,1])[[1]] - summary(df.param.list[,1])[[1]]*.3
min.Cw <- summary(df.param.list[,2])[[1]] - summary(df.param.list[,2])[[1]]*.3
min.LAI <- summary(df.param.list[,3])[[1]] - summary(df.param.list[,3])[[1]]*.3
max.Cab <- summary(df.param.list[,1])[[6]] + summary(df.param.list[,1])[[6]]*.3
max.Cw <- summary(df.param.list[,2])[[6]] + summary(df.param.list[,2])[[6]]*.3
max.LAI <- summary(df.param.list[,3])[[6]] + summary(df.param.list[,3])[[6]]*.3

param.maxmin <- matrix(c(min.Cab,max.Cab,
                         min.Cw ,max.Cw,
                         min.LAI,max.LAI),
                       nrow=3,ncol = 2,byrow = T)

LHS.300 <- Latinhyper(param.maxmin,100*n.params)
LHS.900 <- Latinhyper(param.maxmin,300*n.params)
LHS.1500 <- Latinhyper(param.maxmin,500*n.params)

#using the latin hypercube approach for sample generation
train.par.df.300 <- data.frame(Cab=LHS.300[,1],
                               Cw= LHS.300[,2],
                               LAI=LHS.300[,3])
train.par.df.900 <- data.frame(Cab=LHS.900[,1],
                               Cw= LHS.900[,2],
                               LAI=LHS.900[,3])
train.par.df.1500 <- data.frame(Cab=LHS.1500[,1],
                                Cw= LHS.1500[,2],
                                LAI=LHS.1500[,3])


#we convert it to a spectral object
spclib.plist.300  <-PROSAIL(parameterList = train.par.df.300 )
spclib.plist.900  <-PROSAIL(parameterList = train.par.df.900 )
spclib.plist.1500 <-PROSAIL(parameterList = train.par.df.1500 )

spectra.300  <- as.data.frame(spectra(spclib.plist.300))
spectra.900  <- as.data.frame(spectra(spclib.plist.900))
spectra.1500 <- as.data.frame(spectra(spclib.plist.1500))


#lets combine all the data
df.300  <-  cbind(train.par.df.300,spectra.300)
df.900  <-  cbind(train.par.df.900,spectra.900)
df.1500 <-  cbind(train.par.df.1500,spectra.1500)

#we are training with hyperspectral data
mRF.300 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.300,block.size = 10)
mRF.900 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.900,block.size = 10)
mRF.1500 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.1500,block.size = 10)

#we will train the best with sentinel data.. or we can even try to predict the same traits when trained with
#sentinel data

#now lets generate a spectra of our original traind ata
df.param.list.spclib <- PROSAIL(parameterList=df.param.list)
df.param.list.spectr <- as.data.frame(spectra(df.param.list.spclib))

#we can now try to predict based on each model
mRF.pred.300  <- predict(mRF.300,newdata=df.param.list.spectr )
mRF.pred.900  <- predict(mRF.900,newdata=df.param.list.spectr )
mRF.pred.1500 <- predict(mRF.1500,newdata=df.param.list.spectr )

#now we have to extrac the tables from the predicted model object
mv.mRF.pred.300 <- get.mv.predicted(mRF.pred.300)
mv.mRF.pred.900 <- get.mv.predicted(mRF.pred.900)
mv.mRF.pred.1500 <- get.mv.predicted(mRF.pred.1500)


#and finally we can bring it to the original df and add it to the plot
mv.mRF.pred.300
#just to keep a table of the previous results
out.df.stored <- out.df

out.df$mRF_Cab300 <- mv.mRF.pred.300[,1]
out.df$mRF_Cw300  <- mv.mRF.pred.300[,2]
out.df$mRF_LAI300 <- mv.mRF.pred.300[,3]
out.df$mRF_Cab900 <- mv.mRF.pred.900[,1]
out.df$mRF_Cw900  <- mv.mRF.pred.900[,2]
out.df$mRF_LAI900 <- mv.mRF.pred.900[,3]
out.df$mRF_Cab1500 <- mv.mRF.pred.1500[,1]
out.df$mRF_Cw1500  <- mv.mRF.pred.1500[,2]
out.df$mRF_LAI1500 <- mv.mRF.pred.1500[,3]


par(mfrow=c(2,4))
plot(out.df$Cab,out.df$Cab,pch=19,xlab="Reference",ylab="Prediction",main="Cab",cex=1.5)
points(out.df$Cab,out.df$SCE_Cab,pch=2,col="blue",cex=.75)
points(out.df$Cab,out.df$LBFGSB_Cab,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_Cab~out.df$Cab),col="blue")
abline(lm(out.df$LBFGSB_Cab~out.df$Cab),col="red")

plot(out.df$Cw,out.df$Cw,pch=19,xlab="Reference",ylab="Prediction",main="Cw",cex=1.5)
points(out.df$Cw,out.df$SCE_Cw,pch=2,col="blue",cex=.75)
points(out.df$Cw,out.df$LBFGSB_Cw,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_Cw~out.df$Cw),col="blue")
abline(lm(out.df$LBFGSB_Cw~out.df$Cw),col="red")

plot(out.df$LAI,out.df$SCE_LAI,pch=19,xlab="Reference",ylab="Prediction",main="LAI",cex=1.5)
points(out.df$LAI,out.df$SCE_LAI,pch=2,col="blue",cex=.75)
points(out.df$LAI,out.df$LBFGSB_LAI,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_LAI~out.df$LAI),col="blue")
abline(lm(out.df$LBFGSB_LAI~out.df$LAI),col="red")

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','SCE','L-BFGS-B'), 
       pch=c(19,2,3), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red'))


plot(out.df$Cab,out.df$Cab,pch=19,xlab="Reference",ylab="Prediction",main="Cab",cex=1.5)
points(out.df$Cab,out.df$mRF_Cab300,pch=13,col="blue",cex=.75)
points(out.df$Cab,out.df$mRF_Cab900,pch=15,col="red",cex=.75)
points(out.df$Cab,out.df$mRF_Cab1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_Cab300~out.df$Cab),col="blue")
abline(lm(out.df$mRF_Cab900~out.df$Cab),col="red")
abline(lm(out.df$mRF_Cab1500~out.df$Cab),col="green")

plot(out.df$Cw,out.df$Cw,pch=19,xlab="Reference",ylab="Prediction",main="Cw",cex=1.5)
points(out.df$Cw,out.df$mRF_Cw300,pch=13,col="blue",cex=.75)
points(out.df$Cw,out.df$mRF_Cw900,pch=15,col="red",cex=.75)
points(out.df$Cw,out.df$mRF_Cw1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_Cw300~out.df$Cw),col="blue")
abline(lm(out.df$mRF_Cw900~out.df$Cw),col="red")
abline(lm(out.df$mRF_Cw1500~out.df$Cw),col="green")

plot(out.df$LAI,out.df$SCE_LAI,pch=19,xlab="Reference",ylab="Prediction",main="LAI",cex=1.5)
points(out.df$LAI,out.df$mRF_LAI300,pch=13,col="blue",cex=.75)
points(out.df$LAI,out.df$mRF_LAI900,pch=15,col="red",cex=.75)
points(out.df$LAI,out.df$mRF_LAI1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_LAI300~out.df$LAI),col="blue")
abline(lm(out.df$mRF_LAI900~out.df$LAI),col="red")
abline(lm(out.df$mRF_LAI1500~out.df$LAI),col="green")

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','mRF (n=300)','mRF (n=900)','mRF (n=1500)'), 
       pch=c(19,13,15,23), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red','green'))

