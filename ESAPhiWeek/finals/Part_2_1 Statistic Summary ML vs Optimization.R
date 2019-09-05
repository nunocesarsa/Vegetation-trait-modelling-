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

#library for optimization procedure - you still should read what this actually does lol
library(SoilHyP)

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
dump.ml.opti.fld <- "./Out_OptimVsML"
dump.ml.fld <- "./Out_MachL"
dump.fld <- "./Out_Optim"



## loading the data
out.df <- read.csv(paste(dump.ml.opti.fld,"S2_Opti_Vs_ML.csv",sep = "/"))[,-1]

head(out.df)

sum.out.df <- data.frame(Statistic=rep(NA,8))
shape.df <- sum.out.df

sum.out.df$Statistic[1]<-"Slope"
sum.out.df$Statistic[2]<-"Intercept"
sum.out.df$Statistic[3]<-"R-sqr"
sum.out.df$Statistic[4]<-"MAE"
sum.out.df$Statistic[5]<-"MAPE"
sum.out.df$Statistic[6]<-"SMAPE"
sum.out.df$Statistic[7]<-"MSE"
sum.out.df$Statistic[8]<-"RMSE"

#storing SCE results per traits
names(out.df)
for (i in 1:5){
  x = out.df[,i]
  y = out.df[,i+5]
  
  trait.pred <- names(out.df)[i+5]
  print(names(out.df)[i+5])
  
  temp.lm <- lm(y~x)
  temp.lm$coefficients[[1]]
  temp.lm$coefficients[[2]]
  
  #storing results
  temp.df <- shape.df
  names(temp.df)<-trait.pred
  
  temp.df[1,1] <- temp.lm$coefficients[[2]]
  temp.df[2,1] <- temp.lm$coefficients[[1]]
  temp.df[3,1] <- summary(temp.lm)$r.squared
  
  temp.df[4,1] <- MAE(temp.lm)
  temp.df[5,1] <- MAPE(temp.lm)
  temp.df[6,1] <- SMAPE(temp.lm)
  temp.df[7,1] <- MSE(temp.lm)
  temp.df[8,1] <- RMSE(temp.lm)
  
  sum.out.df <- cbind(sum.out.df,temp.df)
}

#storing Ann results per traits
#names(out.df)
#head(out.df)
for (i in 1:5){
  x = out.df[,i]
  y = out.df[,i+10]
  
  trait.pred <- names(out.df)[i+10]
  print(names(out.df)[i+10])
  
  temp.lm <- lm(y~x)
  temp.lm$coefficients[[1]]
  temp.lm$coefficients[[2]]
  
  #storing results
  temp.df <- shape.df
  names(temp.df)<-trait.pred
  
  temp.df[1,1] <- temp.lm$coefficients[[2]]
  temp.df[2,1] <- temp.lm$coefficients[[1]]
  temp.df[3,1] <- summary(temp.lm)$r.squared
  
  temp.df[4,1] <- MAE(temp.lm)
  temp.df[5,1] <- MAPE(temp.lm)
  temp.df[6,1] <- SMAPE(temp.lm)
  temp.df[7,1] <- MSE(temp.lm)
  temp.df[8,1] <- RMSE(temp.lm)
  
  sum.out.df <- cbind(sum.out.df,temp.df)
}


write.csv(sum.out.df,paste(dump.ml.opti.fld,"S2_Opti_Vs_ML_StatSummary.csv",sep = "/"))


