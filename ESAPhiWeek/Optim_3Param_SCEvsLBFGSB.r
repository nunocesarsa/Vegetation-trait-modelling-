#library for PROSAIL and also Spectral Angle Mapper

library(hsdar)

#library for optimization procedure - you still should read what this actually does lol
library(SoilHyP)

#library that has some extra random generation tools
library(MCMCglmm)

#This should help increase the CPU cores used
library(snowfall)
library(parallel)
sfParallel() #should respond 1 cpu
sfCpus()
detectCores()
#lets use always half the CPU core available
half_cpu_cores <- round(detectCores()/2)
sfInit(parallel = TRUE, cpus = half_cpu_cores, slaveOutfile = "slave.out",nostart = F)
#stopifnot( sfCpus() == 4 )
#stopifnot( sfParallel() == TRUE )
#sfStop()


#first step, creating a function to minimize
#we minimize the Spectral angle mapper - meaning, the difference between two generated spectra from PROSAIL model
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"



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



out.df

warnings()
out.df.omit <- na.omit(out.df)

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
