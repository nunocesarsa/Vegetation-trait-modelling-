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


test.df <- data.frame(InitialPos = rep(NA,nrow(init.par)),
                      p_N = NA, d_N = NA,
                      p_Cab = NA,d_Cab = NA,
                      p_Car = NA,d_Car =NA,
                      p_Cw = NA, d_Cw =NA,
                      p_Cm = NA, d_Cm =NA,
                      p_LAI = NA,d_LAI =NA,
                      p_lidfa = NA,d_Car =NA,
                      niter = NA)
test.df
for (i in 1:nrow(init.par)){
  print(paste("testing for initial position nr: ",i))
  
  
  c1 <- SCEoptim(FUN = sam.fun.spclib,
                 par = as.numeric(df.initParamsTest[i,]), #starts in position one
                 #method="L-BFGS-B",
                 lower=lower_lim,
                 upper=upper_lim,
                 target.spectra=PROSAIL(objective.param[i,]),
                 control = list(ncomplex = 5))
  
  test.df$InitialPos <- 1
  test.df$p_N[i] <- c1$par[1]
  test.df$d_N[i] <- abs(c1$par[1] - df.initParamsTest[i,])
  test.df$p_Cab[i] <- c1$par[2]
  test.df$d_Cab[i] <- abs(c1$par[2] - df.initParamsTest[i,])
  test.df$p_Car[i] <- c1$par[3]
  test.df$d_Car[i] <- abs(c1$par[3] - df.initParamsTest[i,])
  test.df$p_Cw[i] <- c1$par[4]
  test.df$d_Cw[i] <- abs(c1$par[4] - df.initParamsTest[i,])
  test.df$p_Cm[i] <- c1$par[5]
  test.df$d_Cm[i] <- abs(c1$par[5] - df.initParamsTest[i,])
  test.df$p_LAI[i] <- c1$par[6]
  test.df$d_LAI[i] <- abs(c1$par[6] - df.initParamsTest[i,])
  test.df$p_lidfa[i] <- c1$par[7]
  test.df$d_lidfa[i] <- abs(c1$par[7] - df.initParamsTest[i,])
  test.df$niter <- c1$iterations

}








df.initParamsTest[,]
objective.param
test.df

c1$iterations
c1$par[1]

sam.fun.spclib ()
init.par

#minimizes to a target 1 (Cab = 50 - all else default)
init.prosail <- PROSAIL(Cab = 50)
c2 <- SCEoptim(FUN = sam.fun.spclib,
               par = init.par,
               #method="L-BFGS-B",
               lower=lower_lim,
               upper=upper_lim,
               target.spectra=init.prosail)
c2$par

#minimizes to target 2 (Cab = 50, Cw=0.02, LAI =3)
init.prosail <- PROSAIL(Cab=50,Cw=0.02,LAI=3)
c3 <- SCEoptim(FUN = sam.fun.spclib,
               par = init.par,
               #method="L-BFGS-B",
               lower=lower_lim,
               upper=upper_lim,
               target.spectra=init.prosail)
c3$par



#lets make an iterative example of the minimization - might take some time
#first creating a set of input parameters
prosail.runs <- 10
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
out.df$Cab_pred <- NA
out.df$Cw_pred <- NA
out.df$LAI_pred <- NA