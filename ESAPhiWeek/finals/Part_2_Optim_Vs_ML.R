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

dir.create(dump.ml.opti.fld) #it gives out a warning if the folder exits
#setting seed
set.seed(2000)


#first we generate a dataset which we will divide in 3 parts:
#now we set up the upper and lower limits (important! to ensure convergence)
#param.maxmin <- matrix(c(10,80, #Cab
#                         5,40, #Car
#                         0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
#                         0.01,0.02, #Cm
#                         0.5,7),#LAI
#                       #0.05,0.1), #hotstop
#                       nrow=5,ncol = 2,byrow = T)

param.maxmin <- matrix(c(5,100, #Cab
                         5,50, #Car
                         0.005,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.03, #Cm
                         0.5,9.5),#LAI
                       #0.05,0.1), #hotstop
                       nrow=5,ncol = 2,byrow = T)

#now we create a set of initial parameters to be able to start prosail
LHS.3000 <- Latinhyper(param.maxmin,3000)
LHS.6000 <- Latinhyper(param.maxmin,6000)

#stting up data for ML model - ANN
train.par.df.3000 <- data.frame(Cab=LHS.3000[,1],
                                Car=LHS.3000[,2],
                                Cw= LHS.3000[,3],
                                Cm= LHS.3000[,4],
                                LAI=LHS.3000[,5])

train.par.df.6000 <- data.frame(Cab=LHS.3000[,1],
                           Car=LHS.3000[,2],
                           Cw= LHS.3000[,3],
                           Cm= LHS.3000[,4],
                           LAI=LHS.3000[,5])

spclib.plist.3000 <-PROSAIL(parameterList = train.par.df.3000 )
spclib.plist.6000 <-PROSAIL(parameterList = train.par.df.6000)

s2.spclib.3000 <- spectralResampling(spclib.plist.3000,"Sentinel2",response_function = T)
s2.spclib.6000 <- spectralResampling(spclib.plist.6000,"Sentinel2",response_function = T)

s2.spclib.3000.20m <- s2.spclib.3000[,c(2,3,4,5,6,7,9,12,13)]
s2.spclib.6000.20m <- s2.spclib.6000[,c(2,3,4,5,6,7,9,12,13)]

s2.spectra.3000 <- as.data.frame(spectra(s2.spclib.3000.20m))
s2.spectra.6000 <- as.data.frame(spectra(s2.spclib.6000.20m))

s2.band.names <- c("B02","B03","B04",
                   "B05","B06","B07",
                   "B8A","B11","B12")

names(s2.spectra.3000) <- s2.band.names
names(s2.spectra.6000) <- s2.band.names

df.3000 <- cbind(train.par.df.3000,s2.spectra.3000)
df.6000 <- cbind(train.par.df.6000,s2.spectra.6000)

#setting up for ANN
Y.mat.3000 <- as.matrix(df.3000[,c(1:5)])
X.mat.3000 <- as.matrix(df.3000[,c(6:ncol(df.3000))])

Y.mat.6000 <- as.matrix(df.6000[,c(1:5)])
X.mat.6000 <- as.matrix(df.6000[,c(6:ncol(df.6000))])

#ANN structure and hyperparameters
net.st <- c(10,6)
act.fn <- c("tanh","relu")
lrates <- 0.005
epochs <- 3000
optype <- "adam"

#creating and predicting the ANN
ann.3000 <- neuralnetwork(X=X.mat.3000,y=Y.mat.3000,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype)

ann.6000 <- neuralnetwork(X=X.mat.6000,y=Y.mat.6000,regression=T,
                          hidden.layers =net.st,
                          loss.type = "squared",
                          activ.functions = act.fn,
                          n.epochs = epochs,
                          standardize = T,
                          learn.rates = lrates,
                          optim.type = optype)

#once we have a neural network set, we can fetch the ouput of the ml
df.optim <- read.csv(paste(dump.fld,"Optim_S2_5Param_outdata_BroadL.csv",sep="/"))

head(df.optim)
#doing it in one go
df.optim.spectra <- as.data.frame(spectra(spectralResampling(PROSAIL(parameterList = df.optim[,c(2:6)]),
                                                             "Sentinel2",response_function = T)[,c(2,3,4,5,6,7,9,12,13)]))


names(df.optim.spectra) <- s2.band.names

X.mat.optim <- as.matrix(df.optim.spectra)

#predicting
ann.3000.pred <- predict(ann.3000,X.mat.optim)
ann.6000.pred <- predict(ann.6000,X.mat.optim)

#when i optimized i used a different order for variables so now we have to re-order before bringin it all together
colnames(ann.3000.pred$predictions)
ann.3000.pred.reorder <- ann.3000.pred$predictions[,c(1,3,5,4,2)]
ann.6000.pred.reorder <- ann.6000.pred$predictions[,c(1,3,5,4,2)]

df.ann.3000 <- as.data.frame(ann.3000.pred.reorder) 
df.ann.6000 <- as.data.frame(ann.6000.pred.reorder) 

names(df.ann.3000) <- paste("ann3",names(df.ann.3000),sep="_")
names(df.ann.6000) <- paste("ann6",names(df.ann.6000),sep="_")

#we will now bring it all together
names(df.optim)
df.final <- df.optim[,c(2:6,7:11)]
df.final <- cbind(df.final,df.ann.3000)
df.final <- cbind(df.final,df.ann.6000)

#in a massive splot
plot(df.final)

head(df.final)

#lets do a plt each time

plot.lst <- {}
#plot.lst <- vector('list', 5)

for (i in 1:5){
  out.df <- df.final[,c(i,i+5,i+10,i+15)]
  
  print(head(out.df))
  p <- ggplot(out.df, aes(x=out.df[,1])) + 
    #line of values
    geom_point(aes(y = out.df[,1]), color = "black",size=1.5) + 
    geom_smooth(aes(x=out.df[,1],y = out.df[,1]),method = "lm",color = "black",linetype = "dashed",size=.5)+
    #line of optim
    geom_point(aes(y = out.df[,2]), color = "darkred",size=1.5)+
    geom_smooth(aes(x=out.df[,1],y = out.df[,2]),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
    #line of ann3
    geom_point(aes(y = out.df[,3]), color = "blue4",size=1.5)+
    geom_smooth(aes(x=out.df[,1],y = out.df[,3]),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
    #line of ann4
    #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
    #geom_smooth(aes(x=out.df[,1],y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
    theme_bw()
  
  print(p)
  plot.lst[[i]] <- p
}

#its not workiing properly so we have to do one by one
head(df.final)

theme_update(plot.title = element_text(hjust = 0.5))
p.cab <- ggplot(df.final, aes(x=Cab)) + 
  #line of values
  geom_point(aes(y = Cab), color = "black",size=1.5) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_Cab), color = "darkred",size=1.5)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  #line of ann3
  geom_point(aes(y = ann3_Cab), color = "blue4",size=1.5)+
  geom_smooth(aes(x=Cab,y = ann3_Cab),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=Cab,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Cab")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p.cab

p.Cw <- ggplot(df.final, aes(x=Cw)) + 
  #line of values
  geom_point(aes(y = Cw), color = "black",size=1.5) + 
  geom_smooth(aes(x=Cw,y = Cw),method = "lm",color = "black",linetype = "dashed",size=.5)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_Cw), color = "darkred",size=1.5)+
  geom_smooth(aes(x=Cw,y = SCE_Cw),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  #line of ann3
  geom_point(aes(y = ann3_Cw), color = "blue4",size=1.5)+
  geom_smooth(aes(x=Cw,y = ann3_Cw),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=Cw,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Cw")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p.Cw

p.LAI <- ggplot(df.final, aes(x=LAI)) + 
  #line of values
  geom_point(aes(y = LAI), color = "black",size=1.5) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_LAI), color = "darkred",size=1.5)+
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  #line of ann3
  geom_point(aes(y = ann3_LAI), color = "blue4",size=1.5)+
  geom_smooth(aes(x=LAI,y = ann3_LAI),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=LAI,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("LAI")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p.LAI

p.Cm <- ggplot(df.final, aes(x=Cm)) + 
  #line of values
  geom_point(aes(y = Cm), color = "black",size=1.5) + 
  geom_smooth(aes(x=Cm,y = Cm),method = "lm",color = "black",linetype = "dashed",size=.5)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_Cm), color = "darkred",size=1.5)+
  geom_smooth(aes(x=Cm,y = SCE_Cm),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  #line of ann3
  geom_point(aes(y = ann3_Cm), color = "blue4",size=1.5)+
  geom_smooth(aes(x=Cm,y = ann3_Cm),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=Cm,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Cm")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p.Cm

p.Car <- ggplot(df.final, aes(x=Car)) + 
  #line of values
  geom_point(aes(y = Car), color = "black",size=1.5) + 
  geom_smooth(aes(x=Car,y = Car),method = "lm",color = "black",linetype = "dashed",size=.5)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_Car), color = "darkred",size=1.5)+
  geom_smooth(aes(x=Car,y = SCE_Car),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  #line of ann3
  geom_point(aes(y = ann3_Car), color = "blue4",size=1.5)+
  geom_smooth(aes(x=Car,y = ann3_Car),method = "lm",color = "blue4",linetype = "dashed",size=.5)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=Car,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Car")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p.Car

library(grid)
legend <- cowplot::get_legend(p.Car)
grid.newpage()
grid.draw(legend)

library(gridExtra)

tiff(paste(dump.ml.opti.fld,"S2_Opti_Vs_ML.tif",sep = "/"),
     units="px", width = 2048, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(p.cab,p.Cw,
             p.LAI,p.Cm,
             p.Car,
             nrow = 2,ncol=3,
             top="Shuffled Complex Evolution vs ANN")

dev.off()

write.csv(df.final,paste(dump.ml.opti.fld,"S2_Opti_Vs_ML.csv",sep = "/"))

#############################################################
############## Working on neon data stars here ##############
############### But for now only for ML #####################
#############################################################



