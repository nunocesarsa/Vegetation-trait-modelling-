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
dump.esa.results <- "./Out_ESA"
dump.neon.results <- "./Out_Neon"
dump.ml.opti.fld <- "./Out_OptimVsML"
dump.ml.fld <- "./Out_MachL"
dump.fld <- "./Out_Optim"


out.df <- read.csv(paste(dump.esa.results,"ESA_Preds_new.csv",sep="/"))

names(out.df)
out.df <- out.df[,c(2,5,13)]
names(out.df)<- c("LAI","SCE_LAI","ANN_LAI")

library(ggplot2)
library(ggthemes)

p.sce <- ggplot(out.df, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_LAI), color = "orange",size=2) + 
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "orange",linetype = "dashed",size=.5)+
  ylab("Predicted LAI")+
  theme_bw()



p.mlr <- ggplot(out.df, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=1) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = ANN_LAI), color="steelblue",size=2)+
  geom_smooth(aes(x=LAI,y = ANN_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  ylab("")+
  theme_bw()

library(gridExtra)
grid.arrange(p.sce,p.mlr,
             nrow = 2,ncol=1)

tiff(paste(dump.esa.results,"OVP_LAI_Prediction_Horiz.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))


plot(grid.arrange(p.sce,p.mlr,
                  nrow = 2,ncol=1))

dev.off()


tiff(paste(dump.esa.results,"OVP_LAI_Prediction_SCE.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))


plot(p.sce)

dev.off()


tiff(paste(dump.esa.results,"OVP_LAI_Prediction_MLR.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))


plot(p.mlr)

dev.off()

############################
########## storing correlation statistics


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
for (i in 1:1){
  x = out.df[,i]
  y = out.df[,i+1]
  
  trait.pred <- names(out.df)[i+1]
  print(names(out.df)[i+1])
  
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
for (i in 1:1){
  x = out.df[,i]
  y = out.df[,i+2]
  
  trait.pred <- names(out.df)[i+2]
  print(names(out.df)[i+2])
  
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

write.csv(sum.out.df,paste(dump.esa.results,"S2_ESA_StatSummary.csv",sep = "/"))
