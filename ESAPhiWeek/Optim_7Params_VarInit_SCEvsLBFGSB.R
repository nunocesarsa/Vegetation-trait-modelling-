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
gc()
#lets use always half the CPU core available
half_cpu_cores <- round(detectCores()/2)
sfInit(parallel = TRUE, cpus = half_cpu_cores, slaveOutfile = "slave.out",nostart = F)
#stopifnot( sfCpus() == 4 )
#stopifnot( sfParallel() == TRUE )
#sfStop()

#we will first initiate a minimization procedure to estimate the 7 parameters on 50 random points and count the time to finish that process
#then we will do the same with mRF approach
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"



#first step, creating a function to minimize
#we minimize the Spectral angle mapper - meaning, the difference between two generated spectra from PROSAIL model


#this function of SAM shouls compare agasint a given SPECTRAL LIBRARY
sam.fun.spclib <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  #i could probably just change this to make the input more reasonable but cba
  init.param <- data.frame(N=list.params[1],
                           Cab=list.params[2],
                           Car=list.params[3],
                           Cw=list.params[4],
                           Cm=list.params[5],
                           LAI=list.params[6],
                           lidfa=list.params[7],
                           TypeLidf = 0 ,
                           tts = 0,
                           tto = 30)
  
  #this fetches new prosail parameters.. if you want to minimize towards an observation, then this, should be the observation
  
  if (is.null(target.spectra)==T){
    print("Warning: its minimizing to the default PROSAIL parameters") #uncomment if you want to get a spammy console
    
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
param.maxmin <- matrix(c(0.8,2.5, #N, leaf layers
                         5,80, #Cab
                         1,25, #Car
                         0.005,0.02, #Cw
                         0.005,0.02, #Cm
                         0.1,8,#LAI
                         0,90), #mean leaf angle distribution
                       nrow=7,ncol = 2,byrow = T)

param.maxmin <- matrix(c(1,1.5, #N, leaf layers
                         40,60, #Cab
                         10,20, #Car
                         0.01,0.02, #Cw
                         0.01,0.02, #Cm
                         2,5,#LAI
                         30,60), #mean leaf angle distribution
                       nrow=7,ncol = 2,byrow = T)



#this creates a set of target spectra we would like to predict
valid.n <- 10 #* nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin,valid.n)

valid.param.table <- data.frame(N=valid.LHS[,1],
                                Cab=valid.LHS[,2],
                                Car=valid.LHS[,3],
                                Cw=valid.LHS[,4],
                                Cm=valid.LHS[,5],
                                LAI=valid.LHS[,6],
                                lidfa=valid.LHS[,7],
                                TypeLidf = 0 ,
                                tts = 0,
                                tto = 30)
                             

#first we create our initial set of parameters, this could be anything random, hopefully within the limits
init.LHS <- Latinhyper(param.maxmin,10) #this generates 10 different possible starting points

init.par <- data.frame(N=init.LHS[,1],
                       Cab=init.LHS[,2],
                       Car=init.LHS[,3],
                       Cw=init.LHS[,4],
                       Cm=init.LHS[,5],
                       LAI=init.LHS[,6],
                       lidfa=init.LHS[,7])

#minimizes to default prosail (notice im varying only 3 params) (Cab = 40,Cw = 0.01,LAI = 1)

#now we set up the upper and lower limits (important! to ensure convergence)
lower_lim=param.maxmin[,1]
upper_lim=param.maxmin[,2]

#testing for the first validation
c1 <- SCEoptim(FUN = sam.fun.spclib,
               par = as.numeric(init.par[1,]), #starts in position one
               #method="L-BFGS-B",
               lower=lower_lim,
               upper=upper_lim,
               target.spectra=PROSAIL(valid.param.table[1,]),
               control = list(ncomplex = 10)) #this should be the first spectra created with our validation dataset

#lets see what happens given different starting points
df.initParamsTest <- init.par
#this stupid loop creates a test target data set (all rows are the same, we are just seeying the impact of initial positions)
objective.param <- valid.param.table[1,]
for (i in 1:nrow(init.par)){objective.param <- rbind(objective.param,valid.param.table[1,])}
objective.param <- objective.param[-1,]

#for checking time
ptm <- proc.time()

test.df <- data.frame(InitialPos = rep(NA,nrow(init.par)),
                      p_N = NA, d_N = NA,
                      p_Cab = NA,d_Cab = NA,
                      p_Car = NA,d_Car =NA,
                      p_Cw = NA, d_Cw =NA,
                      p_Cm = NA, d_Cm =NA,
                      p_LAI = NA,d_LAI =NA,
                      p_lidfa = NA,d_lidfa =NA,
                      niter = NA)

test.df.LBFGSB <- data.frame(InitialPos = rep(NA,nrow(init.par)),
                      p_N = NA, d_N = NA,
                      p_Cab = NA,d_Cab = NA,
                      p_Car = NA,d_Car =NA,
                      p_Cw = NA, d_Cw =NA,
                      p_Cm = NA, d_Cm =NA,
                      p_LAI = NA,d_LAI =NA,
                      p_lidfa = NA,d_lidfa =NA,
                      niter = NA)

for (i in 1:nrow(init.par)){
  
  
  print(paste("SCE optimiation for set: ",i))
  c1 <- SCEoptim(FUN = sam.fun.spclib,
                 par = as.numeric(df.initParamsTest[i,]), #starts in position one
                 #method="L-BFGS-B",
                 lower=lower_lim,
                 upper=upper_lim,
                 target.spectra=PROSAIL(objective.param[i,]),
                 control = list(ncomplex = 20))
  
  test.df$InitialPos <- 1
  test.df$p_N[i] <- c1$par[1]
  test.df$d_N[i] <- abs(c1$par[1] - objective.param[i,1])
  test.df$p_Cab[i]<- c1$par[2]
  test.df$d_Cab[i]<- abs(c1$par[2] - objective.param[i,2])
  test.df$p_Car[i]<- c1$par[3]
  test.df$d_Car[i]<- abs(c1$par[3] - objective.param[i,3])
  test.df$p_Cw[i] <- c1$par[4]
  test.df$d_Cw[i] <- abs(c1$par[4] - objective.param[i,4])
  test.df$p_Cm[i] <- c1$par[5]
  test.df$d_Cm[i] <- abs(c1$par[5] - objective.param[i,5])
  test.df$p_LAI[i]<- c1$par[6]
  test.df$d_LAI[i]<- abs(c1$par[6] - objective.param[i,6])
  test.df$p_lidfa[i]<- c1$par[7]
  test.df$d_lidfa[i]<- abs(c1$par[7] - objective.param[i,7])
  test.df$niter <- c1$iterations
  
  
  
  print(paste("L-BFGS-B optimiation for set: ",i))
  c2 <-  optim(as.numeric(df.initParamsTest[i,]), sam.fun.spclib, gr = NULL,
               target.spectra=PROSAIL(objective.param[i,]),
               method = c( "L-BFGS-B"),
               lower = lower_lim, upper = upper_lim,
               control = list(), hessian = FALSE)
  
  
  test.df.LBFGSB$InitialPos <- 1
  test.df.LBFGSB$p_N[i] <- c2$par[1]
  test.df.LBFGSB$d_N[i] <- abs(c2$par[1] - objective.param[i,1])
  test.df.LBFGSB$p_Cab[i]<- c2$par[2]
  test.df.LBFGSB$d_Cab[i]<- abs(c2$par[2] - objective.param[i,2])
  test.df.LBFGSB$p_Car[i]<- c2$par[3]
  test.df.LBFGSB$d_Car[i]<- abs(c2$par[3] - objective.param[i,3])
  test.df.LBFGSB$p_Cw[i] <- c2$par[4]
  test.df.LBFGSB$d_Cw[i] <- abs(c2$par[4] - objective.param[i,4])
  test.df.LBFGSB$p_Cm[i] <- c2$par[5]
  test.df.LBFGSB$d_Cm[i] <- abs(c2$par[5] - objective.param[i,5])
  test.df.LBFGSB$p_LAI[i]<- c2$par[6]
  test.df.LBFGSB$d_LAI[i]<- abs(c2$par[6] - objective.param[i,6])
  test.df.LBFGSB$p_lidfa[i]<- c2$par[7]
  test.df.LBFGSB$d_lidfa[i]<- abs(c2$par[7] - objective.param[i,7])
  test.df.LBFGSB$niter <- c2$counts[[1]]

}

write.csv(test.df,
          paste(dump.fld,
                "Inverse_optimization_SCEoptim_TighterL_Output.csv",sep="/"))
write.csv(test.df.LBFGSB,
          paste(dump.fld,
                "Inverse_optimization_LBFGSB_TighterL_Output.csv",sep="/"))
write.csv(objective.param,
          paste(dump.fld,
                "Inverse_optimization_Targets_TighterL_.csv",sep="/"))

#ploting the results of different starting positions vs 7 param traits optimization with broad limits
test.df.SCEOpt <- read.csv(paste(dump.fld,
                            "Inverse_optimization_SCEoptim_BroadL_Output.csv",sep="/"))
test.df.LBFGSB <- read.csv(paste(dump.fld,
                                  "Inverse_optimization_LBFGSB_BroadL_Output.csv",sep="/"))
objective.param <- read.csv(paste(dump.fld,
                                  "Inverse_optimization_Targets_BroadL.csv",sep="/"))

names(objective.param)
names(test.df.LBFGSB)
names(test.df.SCEOpt)
head(objective.param)
dim(test.df.SCEOpt)
seq(3,14,by=2)

#selecting some variables ot make it easier to plot
objective.param.sel <- objective.param[,2:8]
objective.param.sel
test.df.SCEOpt.sel <- test.df.SCEOpt[,seq(3,15,by=2)]
test.df.SCEOpt.sel
test.df.LBFGSB.sel <- test.df.LBFGSB[,seq(3,15,by=2)]
test.df.LBFGSB.sel


par(mfrow=c(2,4))
for (i in 1:ncol(objective.param.sel)){
  plot(objective.param.sel[1,i],pch=19,cex=1.5,
       main=names(objective.param.sel)[i], #gets the layer names
       xaxt='n',xlab="",ylab="", #removes unwanted stuff
       ylim = c(objective.param.sel[1,i]-objective.param.sel[1,i]*.75,
                objective.param.sel[1,i]+objective.param.sel[1,i]*.75)) #getting a broad interval around each value, e.g. an error of 50% around the value
  
  for (k in 1:nrow(test.df.SCEOpt.sel)){
    points(test.df.SCEOpt.sel[k,i],pch=2,col="blue")
    points(test.df.LBFGSB.sel[k,i],pch=3,col="red")
  }
  
}
#plotting the labels
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','SCE','L-BFGS-B'), 
       pch=c(19,2,3), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red'))
