
l#ibrary for PROSAIL and also Spectral Angle Mapper

library(hsdar)

#library for optimization procedure - you still should read what this actually does lol
library(SoilHyP)

#library that has some extra random generation tools
library(MCMCglmm)

#This should help increase the CPU cores used
library(snowfall)
library(parallel)
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
dir.create("./Out_Optim")
#dump.fld <- "./dumps"
dump.fld <- "./Out_Optim"

#setting seed
set.seed(2000)

#This functions calculates the difference between two different spectra
#the input should be 
sam.fun.spclib.3 <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  init.param <- data.frame(Cab=Cab_init,
                           Cw=Cw_init,
                           LAI=LAI_init)
  
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
    
    target.spectra <- spectralResampling(target.spectra,
                                         "Sentinel2",response_function = TRUE)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    
    my.spectra <- spectralResampling(PROSAIL(parameterList = init.param),
                                     "Sentinel2",response_function = TRUE)
    
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,target.spectra)
    fin.val <- spec.dist[1]
    
  }
  
  
  return(fin.val)
}
sam.fun.spclib.4 <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  Cm_init <- list.params[4]
  init.param <- data.frame(Cab=Cab_init,
                           Cw=Cw_init,
                           LAI=LAI_init,
                           Cm=Cm_init)
  
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
    
    target.spectra <- spectralResampling(target.spectra,
                                         "Sentinel2",response_function = TRUE)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    
    my.spectra <- spectralResampling(PROSAIL(parameterList = init.param),
                                     "Sentinel2",response_function = TRUE)
    
    
    #spectral angle mapper is used here
    spec.dist <- sam(my.spectra,target.spectra)
    fin.val <- spec.dist[1]
    
  }
  
  
  return(fin.val)
}
sam.fun.spclib.5 <- function(list.params,target.spectra=NULL){
  
  #this creates the set of data that we want to minimize to - prosail should receive it as a data.frame
  Cab_init <- list.params[1]
  Cw_init  <- list.params[2]
  LAI_init <- list.params[3]
  Cm_init <- list.params[4]
  Car_init <- list.params[5]
  init.param <- data.frame(Cab=Cab_init,
                           Cw=Cw_init,
                           LAI=LAI_init,
                           Cm=Cm_init,
                           Car=Car_init)
  
  #this fetches new prosail parameters.. if you want to minimize towards an observation, then this, should be the observation
  
  if (is.null(target.spectra)==T){
    #print("Warning: its minimizing to the default PROSAIL parameters") #uncomment if you want to get a spammy console
    
    target.spectra <- spectralResampling(target.spectra,
                                         "Sentinel2",response_function = TRUE)
    
    my.spectra <- PROSAIL(parameterList = init.param)
    
    my.spectra <- spectralResampling(PROSAIL(parameterList = init.param),
                                     "Sentinel2",response_function = TRUE)
    
    
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
param.maxmin <- matrix(c(5,100, #Cab
                         5,50, #Car
                         0.005,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.03, #Cm
                         0.5,9.5),#LAI
                       #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)

min.max.Cab <- param.maxmin[1,]
min.max.Car <- param.maxmin[2,]
min.max.Cw <- param.maxmin[3,]
min.max.Cm <- param.maxmin[4,]
min.max.LAI <- param.maxmin[5,]

#llim - lower limit
llim.3 <- c(param.maxmin[1,1],
            param.maxmin[3,1],
            param.maxmin[5,1])

llim.4 <- c(param.maxmin[1,1],
            param.maxmin[3,1],
            param.maxmin[5,1],
            param.maxmin[4,1])

llim.5 <- c(param.maxmin[1,1],
            param.maxmin[3,1],
            param.maxmin[5,1],
            param.maxmin[4,1],
            param.maxmin[2,1])
#ulim - upper limit
ulim.3 <- c(param.maxmin[1,2],
            param.maxmin[3,2],
            param.maxmin[5,2])

ulim.4 <- c(param.maxmin[1,2],
            param.maxmin[3,2],
            param.maxmin[5,2],
            param.maxmin[4,2])

ulim.5 <- c(param.maxmin[1,2],
            param.maxmin[3,2],
            param.maxmin[5,2],
            param.maxmin[4,2],
            param.maxmin[2,2])

min.max.Cab <- param.maxmin[1,]
min.max.Car <- param.maxmin[2,]
min.max.Cw <- param.maxmin[3,]
min.max.Cm <- param.maxmin[4,]
min.max.LAI <- param.maxmin[5,]
min.max.Cab[1]
#now we create a set of initial parameters to be able to start prosail
init.par.df <- data.frame(Cab=rtnorm(1,mean=mean(min.max.Cab),sd=mean(min.max.Cab)/4,lower=min.max.Cab[1],upper=min.max.Cab[2]), 
                          Cw= rtnorm(1,mean=mean(min.max.Cw),sd=mean(min.max.Cw)/4,lower=min.max.Cw[1],upper=min.max.Cw[2]),
                          LAI=rtnorm(1,mean=mean(min.max.LAI),sd=mean(min.max.LAI )/4,lower=min.max.LAI[1],upper=min.max.LAI[2]),
                          Cm= rtnorm(1,mean=mean(min.max.Cm),sd=mean(min.max.Cm)/4,lower=min.max.Cm[1],upper=min.max.Cm[2]),
                          Car=rtnorm(1,mean=mean(min.max.Car),sd=mean(min.max.Car)/4,lower=min.max.Car[1],upper=min.max.Car[2])) 
              
init.par.3 <- as.numeric(init.par.df[1,1:3])
init.par.4 <- as.numeric(init.par.df[1,1:4])
init.par.5 <- as.numeric(init.par.df[1,])

llim.5
ulim.5


#let's now generate a dataset of the 5 traits
prosail.runs <- 30
df.param.list <- data.frame(Cab=rtnorm(prosail.runs,mean=mean(min.max.Cab),sd=mean(min.max.Cab)/4,lower=min.max.Cab[1],upper=min.max.Cab[2]), 
                            Cw= rtnorm(prosail.runs,mean=mean(min.max.Cw),sd=mean(min.max.Cw)/4,lower=min.max.Cw[1],upper=min.max.Cw[2]),
                            LAI=rtnorm(prosail.runs,mean=mean(min.max.LAI),sd=mean(min.max.LAI )/4,lower=min.max.LAI[1],upper=min.max.LAI[2]),
                            Cm= rtnorm(prosail.runs,mean=mean(min.max.Cm),sd=mean(min.max.Cm)/4,lower=min.max.Cm[1],upper=min.max.Cm[2]),
                            Car=rtnorm(prosail.runs,mean=mean(min.max.Car),sd=mean(min.max.Car)/4,lower=min.max.Car[1],upper=min.max.Car[2])) 


#testing with 3 parameters
out.df.3 <- df.param.list[,c(1:3)]
out.df.3$SCE_Cab <- NA
out.df.3$SCE_Cw  <- NA
out.df.3$SCE_LAI <- NA

out.df.3$LBFGSB_Cab <- NA
out.df.3$LBFGSB_Cw  <- NA
out.df.3$LBFGSB_LAI <- NA

out.df.3$SCE_iter <- NA
out.df.3$LBFGSB_fncall <- NA

for( i in 1:nrow(out.df.3)){
  print(paste("Estimating sample:", i))
  
  #print(df.param.list[i,])
  tgt.spec <- PROSAIL(parameterList = df.param.list[i,])
  
  SCEoptim.pred <- SCEoptim(FUN = sam.fun.spclib.3,
                            par = init.par.3,
                            #method="L-BFGS-B",
                            lower=llim.3,
                            upper=ulim.3,
                            target.spectra=tgt.spec)
  
  LBFGSB.pred <- optim(init.par.3, sam.fun.spclib.3, gr = NULL,
                       target.spectra=tgt.spec,
                       method = c( "L-BFGS-B"),
                       lower = llim.3,
                       upper = ulim.3,
                       control = list(), hessian = FALSE)
  
  #storing the outputs
  out.df.3$SCE_Cab[i] <- SCEoptim.pred$par[1]
  out.df.3$SCE_Cw[i]  <- SCEoptim.pred$par[2]
  out.df.3$SCE_LAI[i] <- SCEoptim.pred$par[3]
  
  out.df.3$LBFGSB_Cab[i] <- LBFGSB.pred$par[1]
  out.df.3$LBFGSB_Cw[i]  <- LBFGSB.pred$par[2]
  out.df.3$LBFGSB_LAI[i] <- LBFGSB.pred$par[3]
  
  out.df.3$SCE_iter <- SCEoptim.pred$iterations
  out.df.3$LBFGSB_fncall <- LBFGSB.pred$counts[[1]]
  
}

#lets test for 4
out.df.4 <- df.param.list[,c(1:4)]
out.df.4$SCE_Cab <- NA
out.df.4$SCE_Cw  <- NA
out.df.4$SCE_LAI <- NA
out.df.4$SCE_Cm  <- NA

out.df.4$LBFGSB_Cab <- NA
out.df.4$LBFGSB_Cw  <- NA
out.df.4$LBFGSB_LAI <- NA
out.df.4$LBFGSB_Cm  <- NA

out.df.4$SCE_iter <- NA
out.df.4$LBFGSB_fncall <- NA

for( i in 1:nrow(out.df.4)){
  print(paste("Estimating sample:", i))
  
  #print(df.param.list[i,])
  tgt.spec <- PROSAIL(parameterList = df.param.list[i,])
  
  SCEoptim.pred <- SCEoptim(FUN = sam.fun.spclib.4,
                            par = init.par.4,
                            #method="L-BFGS-B",
                            lower=llim.4,
                            upper=ulim.4,
                            target.spectra=tgt.spec)
  
  LBFGSB.pred <- optim(init.par.4, sam.fun.spclib.4, gr = NULL,
                       target.spectra=tgt.spec,
                       method = c( "L-BFGS-B"),
                       lower = llim.4,
                       upper = ulim.4,
                       control = list(), hessian = FALSE)
  
  #storing the outputs
  out.df.4$SCE_Cab[i] <- SCEoptim.pred$par[1]
  out.df.4$SCE_Cw[i]  <- SCEoptim.pred$par[2]
  out.df.4$SCE_LAI[i] <- SCEoptim.pred$par[3]
  out.df.4$SCE_Cm[i] <- SCEoptim.pred$par[4]
  
  out.df.4$LBFGSB_Cab[i] <- LBFGSB.pred$par[1]
  out.df.4$LBFGSB_Cw[i]  <- LBFGSB.pred$par[2]
  out.df.4$LBFGSB_LAI[i] <- LBFGSB.pred$par[3]
  out.df.4$LBFGSB_Cm[i] <- LBFGSB.pred$par[4]
  
  out.df.4$SCE_iter <- SCEoptim.pred$iterations
  out.df.4$LBFGSB_fncall <- LBFGSB.pred$counts[[1]]
  
}

#lets test for 5
out.df.5 <- df.param.list[,c(1:5)]
out.df.5$SCE_Cab <- NA
out.df.5$SCE_Cw  <- NA
out.df.5$SCE_LAI <- NA
out.df.5$SCE_Cm  <- NA
out.df.5$SCE_Car <- NA

out.df.5$LBFGSB_Cab <- NA
out.df.5$LBFGSB_Cw  <- NA
out.df.5$LBFGSB_LAI <- NA
out.df.5$LBFGSB_Cm  <- NA
out.df.5$LBFGSB_Car <- NA

out.df.5$SCE_iter <- NA
out.df.5$LBFGSB_fncall <- NA

#sometimes this loop breaks down, seems to me the optimization calls for parameters completely outside PROSAIL range, generates NaN and thus fails SAM function
#for( i in 14:17){
for( i in 1:nrow(out.df.5)){
  print(paste("Estimating sample:", i))
  
  #print(df.param.list[i,])
  tgt.spec <- PROSAIL(parameterList = df.param.list[i,])
  
  SCEoptim.pred <- SCEoptim(FUN = sam.fun.spclib.5,
                            par = init.par.5,
                            #method="L-BFGS-B",
                            lower=llim.5,
                            upper=ulim.5,
                            target.spectra=tgt.spec)
  
  LBFGSB.pred <- optim(init.par.5, sam.fun.spclib.5, gr = NULL,
                       target.spectra=tgt.spec,
                       method = c( "L-BFGS-B"),
                       lower = llim.5,
                       upper = ulim.5,
                       control = list(), hessian = FALSE)
  
  #storing the outputs
  out.df.5$SCE_Cab[i] <- SCEoptim.pred$par[1]
  out.df.5$SCE_Cw[i]  <- SCEoptim.pred$par[2]
  out.df.5$SCE_LAI[i] <- SCEoptim.pred$par[3]
  out.df.5$SCE_Cm[i]  <- SCEoptim.pred$par[4]
  out.df.5$SCE_Car[i] <- SCEoptim.pred$par[5]
  
  out.df.5$LBFGSB_Cab[i] <- LBFGSB.pred$par[1]
  out.df.5$LBFGSB_Cw[i]  <- LBFGSB.pred$par[2]
  out.df.5$LBFGSB_LAI[i] <- LBFGSB.pred$par[3]
  out.df.5$LBFGSB_Cm[i]  <- LBFGSB.pred$par[4]
  out.df.5$LBFGSB_Car[i] <- LBFGSB.pred$par[5]
  
  out.df.5$SCE_iter <- SCEoptim.pred$iterations
  out.df.5$LBFGSB_fncall <- LBFGSB.pred$counts[[1]]
  
}


write.csv(out.df.3,paste(dump.fld,"Optim_S2_3Param_outdata_BroadL.csv",sep="/"))
write.csv(out.df.4,paste(dump.fld,"Optim_S2_4Param_outdata_BroadL.csv",sep="/"))
write.csv(out.df.5,paste(dump.fld,"Optim_S2_5Param_outdata_BroadL.csv",sep="/"))


########################################
sum.out.df <- data.frame(Trait=c("Cab",
                                 "Cw",
                                 "LAI",
                                 "Cm",
                                 "Car"))
sum.out.df$SCE_optim_3 <- NA
sum.out.df$SCE_optim_4 <- NA
sum.out.df$SCE_optim_5 <- NA

sum.out.df$LBF_optim_3 <- NA
sum.out.df$LBF_optim_4 <- NA
sum.out.df$LBF_optim_5 <- NA

library(DescTools)
for (i in 1:3){
  sum.out.df$SCE_optim_3[i] <- MAPE(ref=out.df.3[,i],x=out.df.3[,i+3])*100
  sum.out.df$LBF_optim_3[i] <- MAPE(ref=out.df.3[,i],x=out.df.3[,i+6])*100
}

for (i in 1:4){
  sum.out.df$SCE_optim_4[i] <- MAPE(ref=out.df.4[,i],x=out.df.4[,i+4])*100
  sum.out.df$LBF_optim_4[i] <- MAPE(ref=out.df.4[,i],x=out.df.4[,i+8])*100
}

for (i in 1:5){
  sum.out.df$SCE_optim_5[i] <- MAPE(ref=out.df.5[,i],x=out.df.5[,i+5])*100
  sum.out.df$LBF_optim_5[i] <- MAPE(ref=out.df.5[,i],x=out.df.5[,i+10])*100
}

write.csv(sum.out.df,paste(dump.fld,"Optim_S2_3Param_MAPE_Summary_BroadL.csv",sep="/"))

######################## plotting stuff
library(ggplot2)
library(ggthemes)
library(gridExtra)


cab.3.SCE <- ggplot(out.df.3, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=2)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()


cab.3.LBF <- ggplot(out.df.3, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=2) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

Cw.3.SCE <- ggplot(out.df.3, aes(x=Cw)) + 
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cw), color = "darkred",size=2)+
  geom_smooth(aes(x=Cw,y = SCE_Cw),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Cw.3.LBF <- ggplot(out.df.3, aes(x=Cw))+
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cw), color="steelblue",size=2) +
  geom_smooth(aes(x=Cw,y = LBFGSB_Cw),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

#LAI

LAI.3.SCE <- ggplot(out.df.3, aes(x=LAI)) + 
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "darkred",size=2)+
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

LAI.3.LBF <- ggplot(out.df.3, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_LAI), color="steelblue",size=2) +
  geom_smooth(aes(x=LAI,y = LBFGSB_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()


grid.arrange(cab.3.SCE, cab.3.LBF,
             Cw.3.SCE, Cw.3.LBF,
             LAI.3.SCE,LAI.3.LBF,
             nrow = 3,ncol=2,
             top="Target: 3 Parameters")

############################### 4 parameters plot

cab.4.SCE <- ggplot(out.df.4, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=2)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()


cab.4.LBF <- ggplot(out.df.4, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=2) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

Cw.4.SCE <- ggplot(out.df.4, aes(x=Cw)) + 
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cw), color = "darkred",size=2)+
  geom_smooth(aes(x=Cw,y = SCE_Cw),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Cw.4.LBF <- ggplot(out.df.4, aes(x=Cw))+
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cw), color="steelblue",size=2) +
  geom_smooth(aes(x=Cw,y = LBFGSB_Cw),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

LAI.4.SCE <- ggplot(out.df.4, aes(x=LAI)) + 
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "darkred",size=2)+
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

LAI.4.LBF <- ggplot(out.df.4, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_LAI), color="steelblue",size=2) +
  geom_smooth(aes(x=LAI,y = LBFGSB_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

Cm.4.SCE <- ggplot(out.df.4, aes(x=Cm)) + 
  geom_point(aes(y = Cm), color = "black",size=1) + 
  geom_smooth(aes(x=Cm,y = Cm),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cm), color = "darkred",size=2)+
  geom_smooth(aes(x=Cm,y = SCE_Cm),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Cm.4.LBF <- ggplot(out.df.4, aes(x=Cm))+
  geom_point(aes(y = Cm), color = "black",size=1) + 
  geom_smooth(aes(x=Cm,y = Cm),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cm), color="steelblue",size=2) +
  geom_smooth(aes(x=Cm,y = LBFGSB_Cm),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()



grid.arrange(cab.4.SCE, cab.4.LBF,
             Cw.4.SCE,  Cw.4.LBF,
             LAI.4.SCE, LAI.4.LBF,
             Cm.4.SCE,  Cm.4.LBF,
             nrow = 4,ncol=2,
             top="Target: 4 Parameters")

#################### plotting 5 parameters

cab.5.SCE <- ggplot(out.df.5, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=2)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()


cab.5.LBF <- ggplot(out.df.5, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=1) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=2) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

Cw.5.SCE <- ggplot(out.df.5, aes(x=Cw)) + 
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cw), color = "darkred",size=2)+
  geom_smooth(aes(x=Cw,y = SCE_Cw),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Cw.5.LBF <- ggplot(out.df.5, aes(x=Cw))+
  geom_point(aes(y = Cw), color = "black",size=1) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cw), color="steelblue",size=2) +
  geom_smooth(aes(x=Cw,y = LBFGSB_Cw),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

LAI.5.SCE <- ggplot(out.df.5, aes(x=LAI)) + 
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "darkred",size=2)+
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

LAI.5.LBF <- ggplot(out.df.5, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_LAI), color="steelblue",size=2) +
  geom_smooth(aes(x=LAI,y = LBFGSB_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()

Cm.5.SCE <- ggplot(out.df.5, aes(x=Cm)) + 
  geom_point(aes(y = Cm), color = "black",size=1) + 
  geom_smooth(aes(x=Cm,y = Cm),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cm), color = "darkred",size=2)+
  geom_smooth(aes(x=Cm,y = SCE_Cm),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Cm.5.LBF <- ggplot(out.df.5, aes(x=Cm))+
  geom_point(aes(y = Cm), color = "black",size=1) + 
  geom_smooth(aes(x=Cm,y = Cm),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cm), color="steelblue",size=2) +
  geom_smooth(aes(x=Cm,y = LBFGSB_Cm),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()


Car.5.SCE <- ggplot(out.df.5, aes(x=Car)) + 
  geom_point(aes(y = Car), color = "black",size=1) + 
  geom_smooth(aes(x=Car,y = Car),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Car), color = "darkred",size=2)+
  geom_smooth(aes(x=Car,y = SCE_Car),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()

Car.5.LBF <- ggplot(out.df.5, aes(x=Car))+
  geom_point(aes(y = Car), color = "black",size=1) + 
  geom_smooth(aes(x=Car,y = Car),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Car), color="steelblue",size=2) +
  geom_smooth(aes(x=Car,y = LBFGSB_Car),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()










#saving the plots
tiff(paste(dump.fld,"Optim_S2_5Params_BroadL.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
grid.arrange(cab.5.SCE, cab.5.LBF,
             Cw.5.SCE,  Cw.5.LBF,
             LAI.5.SCE, LAI.5.LBF,
             Cm.5.SCE,  Cm.5.LBF,
             Car.5.SCE, Car.5.LBF,
             nrow = 5,ncol=2,
             top="Target: 5 Parameters")
dev.off()




tiff(paste(dump.fld,"Optim_S2_4Params_BroadL.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))
grid.arrange(cab.4.SCE, cab.4.LBF,
             Cw.4.SCE,  Cw.4.LBF,
             LAI.4.SCE, LAI.4.LBF,
             Cm.4.SCE,  Cm.4.LBF,
             nrow = 4,ncol=2,
             top="Target: 4 Parameters")
dev.off()


tiff(paste(dump.fld,"Optim_S2_3Params_BroadL.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(cab.3.SCE, cab.3.LBF,
             Cw.3.SCE, Cw.3.LBF,
             LAI.3.SCE,LAI.3.LBF,
             nrow = 3,ncol=2,
             top="Target: 3 Parameters")
dev.off()











