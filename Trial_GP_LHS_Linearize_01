library(DiceDesign)
library(DiceKriging)
library(emoa)
library(GPareto) #make sure you have RTools installed for the c++/fortran etc compilers 

#PS: i recommend installing any dependency that it gives alarm off. 


library(hsdar) #prosail is here
library(ggplot2) #just for plotting

#dont play too much with this..
#Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
#Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/;C:/MinGW/bin/")

#generating an initial set of observations

#lets try to create a parameter list - we want to find 3 things variables
#LAI, Cab and Cw with default everything else

#we generate given a normal distribution - this can be seen as a prior
#LAI cant be less than 0.. so it must be a truncated normal distribution
library(MCMCglmm) #rtnorm comes from here

#lib of gaussian processes
library(mlegp)

#Latin hypercube comes from here
library(FME)




gc()
prosail.runs <- 10
#prosail.runs <- 100

#lets generate 3 matrices of values based on a LHS 
#matrix of format: param1 min,param 1 max ; each line a new parameter
param.maxmin <- matrix(c(20,60,
                         0.01,0.03,
                         0,3),
                       nrow=3,ncol = 2,byrow = T)

LHS <- Latinhyper(param.maxmin,prosail.runs)
#now, each column of the data represents 100 samples of one the parameters

#using random sampling based on normal distributions
param.list <- data.frame(Cab=rnorm(prosail.runs,mean=40,sd=15),
                         Cw=rnorm(prosail.runs,mean=0.02,sd=0.005),
                         LAI=rtnorm(prosail.runs,mean=1,sd=3,lower=0,upper=6))
#using the latin hypercube approach for sample generation
param.list <- data.frame(Cab=LHS[,1],
                         Cw=LHS[,2],
                         LAI=LHS[,3])


#old param list
#param.list <- data.frame(N = c(rep.int(seq(0.5, 1.5, 0.5), 2)),
#                         LAI = c(rep.int(0.5, 3), rep.int(1, 3)))
#FYI
#my.spclib <- PROSAIL(N = 1.5, Cab = 40, Car = 8, Cbrown = 0.0,
#                     Cw = 0.01, Cm = 0.009, psoil = 0, LAI = 1, 
#                     TypeLidf = 1, lidfa = -0.35, lidfb = -0.15,
#                     hspot = 0.01, tts = 30, tto = 10, psi = 0,
#                     parameterList = NULL, rsoil = NULL)

my.spclib.plist <- PROSAIL(parameterList = param.list)
plot(my.spclib.plist)


#before the fitting lets transform all spectras (linearize) - remember to transform the final output before

#Here you can make all transformations but for now only for the 3 parameters of interest

#transformation - imported from python

#Cab T and iT
"# Cab, posn 1
x_out[1] = np.exp(1./100. * x[1])

x_out = x*1.
# Cab, posn 1
x_out[1] = -100.*np.log ( x[1] )"

Cab_T <- function(vals){
  t_vals <- exp(vals/100)
  return(t_vals)
}
Cab_iT <- function(t_vals){
  vals <- 100*log(t_vals) #came with a minus but that is reversing the results...
  return(vals)
}

#Cw T and iT
"# Cw, posn 4
x_out[4] = np.exp(50. *x[4])

# Cw, posn 4
x_out[4] = (-1./50.)*np.log ( x[4] )"

Cw_T <- function(vals){
  t_vals <- exp(50*vals)
  return(t_vals)
}
Cw_iT <- function(t_vals){
  vals <- (1/50)*log(t_vals) #came with a minus but that is reversing the results...
  return(vals)
}


#LAI T and iT
"# LAI, posn 6"
"x_out[6] =  np.exp(1./2. *x[6])

# LAI, posn 6
x_out[6] = -2.*np.log ( x[6] )"

LAI_T <- function(vals){
  t_vals <- exp(vals/2)
  return(t_vals)
}

LAI_iT <- function(t_vals){
  vals <- (2)*log(t_vals) #came with a minus but that is reversing the results...
  return(vals)
}

#now the only change to the previous model is that we fit it towards the linearized values

#fitting a model - OLD ONE - FITS ON Untransformed values
#gp.fit.2 <- mlegp(spectra(my.spclib.plist),
#                  as.matrix(param.list)) #this works already but the results are awful

#transforming the functions
param.list.t <- param.list
param.list.t[,1]<-Cab_T(param.list[,1])
param.list.t[,2]<-Cab_T(param.list[,2])
param.list.t[,3]<-Cab_T(param.list[,3])

gp.fit.2 <- mlegp(spectra(my.spclib.plist),
                  as.matrix(param.list.t))

#not working, tottaly confused how to deal with pcas here - skip this small section
#lets try doing a PC decomposition beforehand - requires decomposing the input data beforehand
#singularValueImportance(t(as.matrix(param.list)))
#gp.fit.2 <- mlegp(spectra(my.spclib.plist),
#                  t(as.matrix(param.list)),PC.num =3)

#this was altered to predict real values - you have to remove the functions if you did not transform them
par(mfrow=c(1,3))
eval.df <- param.list
names(eval.df)
eval.df$Cab_t_selfPred <- Cab_iT(predict(gp.fit.2[[1]]))
eval.df$Cw_selfPred <- Cw_iT(predict(gp.fit.2[[2]]))
eval.df$LAI_selfPred <- LAI_iT(predict(gp.fit.2[[3]]))

#par(mfrow=c(1,3))
plot(eval.df[,1],eval.df[,4],main="Cab (selfpred)",xlab="Prosail",ylab="GP Self")
plot(eval.df[,2],eval.df[,5],main="Cw (selfpred)",xlab="Prosail",ylab="GP Self")
plot(eval.df[,3],eval.df[,6],main="LAI (selfpred)",xlab="Prosail",ylab="GP Self")

plot(eval.df[,1]-eval.df[,4],main="delta_Cab (selfpred)")
plot(eval.df[,2]-eval.df[,5],main="delta_Cw (selfpred)")
plot(eval.df[,3]-eval.df[,6],main="LAI (selfpred)")


#bb <- getSingularValues(t(as.matrix(param.list)))
#bb <- pcweights(t(as.matrix(param.list)), weights.num = NULL, cutoff = 99)


#ok, it seems we can predict our own data but can we predict new data??
#lets generate a new sample of reflecances


param.list.test <- data.frame(Cab=rnorm(prosail.runs,mean=40,sd=15),
                              Cw=rnorm(prosail.runs,mean=0.02,sd=0.005),
                              LAI=rtnorm(prosail.runs,mean=1,sd=3,lower=0,upper=6))

spec.lib.trial <- PROSAIL(parameterList = param.list.test)

new.mat <- cbind(spectra(spec.lib.trial,1:100))
dim(new.mat)

#if we used transformation, we must transform the output

pred.Cab <- Cab_iT(predict(gp.fit.2[[1]],newdata=new.mat))
pred.Cw <- Cw_iT(predict(gp.fit.2[[2]],newdata=new.mat))
pred.LAI <- LAI_iT(predict(gp.fit.2[[3]],newdata=new.mat))

test.df <- param.list.test
test.df$Cab_pred <- pred.Cab
test.df$Cw_pred <- pred.Cw
test.df$LAI_pred <- pred.LAI

par(mfrow=c(1,3))
plot(test.df[,1],test.df[,4],main="Cab newdata",xlab="Prosail",ylab="GP pred")
plot(test.df[,2],test.df[,5],main="Cw newdata",xlab="Prosail",ylab="GP pred")
plot(test.df[,3],test.df[,6],main="LAI newdata",xlab="Prosail",ylab="GP pred")

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(test.df[,1],test.df[,4])
RMSE(test.df[,2],test.df[,5])
RMSE(test.df[,3],test.df[,6])

plot(lm(test.df[,1]~test.df[,4]))

fit.lm <- lm(lm(test.df[,1]~test.df[,4]))
summary(fit.lm)
