#library for PROSAIL and also Spectral Angle Mapper
library(hsdar)

#library for optimization procedure - you still should read what this actually does lol
#library(SoilHyP)

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

gc()

#lets generate a space of varying biophysical parameters 

#we can generate them all together but then we will feed the model 1 by 1 so we can see
#the impact of each trait on the prosail simulation

#limits taken from: https://www.mdpi.com/2072-4292/8/2/119

#N Structure parameter
#Cab chlorophyll content
#Car Carotenoid content
#Cbrown Brown pigment content
#Cw Equivalent water thickness
#Cm Dry matter content
#psoil Dry/Wet soil factor
#LAI Leaf area index

#TypeLidf Type of leaf angle distribution. See details section
#lidfa Leaf angle distribution. See details section
#lidfb Leaf angle distribution. See details section
#hspot Hotspot parameter
#tts Solar zenith angle
#tto Observer zenith angle
#psi Relative azimuth angle
#parameterList An optional object of class 'data.frame'. Function will iterate over rows of parameterList setting missing entries to default values. See examples section.
#rsoil background (soil) reflectance. N



param.maxmin <- matrix(c(0.8, 2.5, #leaf layers or leaf structure
                         0.2,77, #Cab
                         0,15, #Car
                         0.0043,0.0753, #Cw
                         0.0017,0.0331, #Cm
                         0,8), #LAI
                       nrow=6,ncol = 2,byrow = T)

#we can now generate thousands of sampling, lets do 500 per trait: 
prosail.runs <- nrow(param.maxmin)*1000
LHS <- Latinhyper(param.maxmin,prosail.runs)
head(LHS)
summary(LHS)


#for future use
names.traits <- c("N","Cab",
                  "Car","Cw",
                  "Cm","LAI")

#first we run prosail for all combinations to save running time
spclib.container <- list() #empty list to store all data

#now we will run prosail independently with all fixed except one trait and see what happens to the spectra
for (i in 1:ncol(LHS)){
  
  #we generate an empty df and assign the name from the previous list
  param.list <- data.frame(temp=LHS[,i])
  colnames(param.list)<-names.traits[i]
  
  print(paste("Generating set for:", names(param.list)))
  spclib.container <- c(spclib.container,PROSAIL(parameterList = param.list,
                                                 TypeLidf="Erectophile")) #i guess most grass is this - NOT SURE THIS IS PROPER YET
  
}

par(mfrow=c(2,3))
for (i in 1:ncol(LHS)){
  
  #the order matters for plotting because the main graphic should overlap
  
  #first the maximum
  plot(spclib.container[[i]],main=paste("Varying",names.traits[i],
                  ";",param.maxmin[i,1],
                  "to",param.maxmin[i,2]),
       FUN = "max",col="blue",lty=2)
  
  #then the minimum
  plot(spclib.container[[i]],#main=paste("Varying",names.traits[i],
                             #           ";",param.maxmin[i,1],
                             #           "to",param.maxmin[i,2]),
       FUN = "min",col="red",new=FALSE,lty=2)

  #then the average
  plot(spclib.container[[i]],#main=paste("Varying",names.traits[i],
                             #           ";",param.maxmin[i,1],
                             #           "to",param.maxmin[i,2]),
       FUN = "mean",new=FALSE,lty=2) 
  
}

#now lets plot the variance, its a much more direct visual of where the main differences occur.
for (i in 1:ncol(LHS)){
  
  plot(spclib.container[[i]],main=paste("Varying",names.traits[i],
                                        ";",param.maxmin[i,1],
                                        "to",param.maxmin[i,2]),
       FUN="var")
  
}

#this was while using Hyperspectral simulated data, lets see if it was a sentinel response

#first we generate a new list with sentinel responses
s2.spclib.container <- list() 
for (i in 1:ncol(LHS)){
  
  s2.spclib.container <- c(s2.spclib.container,spectralResampling(spclib.container[[i]],
                                                                  "Sentinel2",response_function = TRUE))
}

#now we can redo the same plots but with sentinel data

par(mfrow=c(2,3))
for (i in 1:ncol(LHS)){
  
  #the order matters for plotting because the main graphic should overlap
  
  #first the maximum
  plot(s2.spclib.container[[i]],main=paste("Varying",names.traits[i],
                                        ";",param.maxmin[i,1],
                                        "to",param.maxmin[i,2]),
       FUN = "max",col="blue",lty=2,type="b",pch=3)
  
  #then the minimum
  plot(s2.spclib.container[[i]],#main=paste("Varying",names.traits[i],
       #           ";",param.maxmin[i,1],
       #           "to",param.maxmin[i,2]),
       FUN = "min",col="red",new=FALSE,lty=2,type="b",pch=3)
  
  #then the average
  plot(s2.spclib.container[[i]],#main=paste("Varying",names.traits[i],
       #           ";",param.maxmin[i,1],
       #           "to",param.maxmin[i,2]),
       FUN = "mean",new=FALSE,lty=2,type="b",pch=19) 
  
}

#now lets plot the variance, its a much more direct visual of where the main differences occur.
for (i in 1:ncol(LHS)){
  
  plot(s2.spclib.container[[i]],main=paste("Varying",names.traits[i],
                                        ";",param.maxmin[i,1],
                                        "to",param.maxmin[i,2],type="b"),
       FUN="var")
  
}

#now lets traing a mRF on all these parameters at once and see if we are able to do good predictions

#first we must generate a new sampling of the trait space
param.list <- data.frame(C1=LHS[,1],C2=LHS[,2],
                         C3=LHS[,3],C4=LHS[,4],
                         C5=LHS[,5],C6=LHS[,6],
                         TypeLidf="Erectophile") #this might require change if you change something before
#now the names must be changed
colnames(param.list)<-c(names.traits,"TypeLidf")

mRF.spclib <- PROSAIL(parameterList = param.list)
mRF.s2.spclib <- spectralResampling(mRF.spclib,
                                    "Sentinel2",response_function = TRUE)


#commonly we dont use all of them due to their GSD:
dim(mRF.s2.spclib)

#lets pick only the bands no larger than 20m resolution:
#b 2,3,4,5,6,7,8,8a,11,12
mRF.s2.spc.sel <- mRF.s2.spclib[,c(2,3,4,5,6,7,8,9,12,13)]

head(mRF.s2.spc.sel)

#brings it all together
full.df <- cbind(param.list,as.data.frame(mRF.s2.spc.sel))
names(full.df)
#next step is optional
names(full.df) <- c(names.traits,"TypeLidf",
                    "B02","B03","B04",
                    "B05","B06","B07",
                    "B08","B8A","B11",
                    "B12")

names(full.df)
#before going into the model we should remove the type lidf since we are trying with
#just one type of leaf
full.df <- full.df[,-7]
names(full.df)

#training the mfr
mRF_all <- rfsrc(Multivar(N,Cab,Car,
                          Cw,Cm,LAI)~.,data=full.df,block.size = 500)

#creating a test dataset
#let's see if the prediction against a test data is good first
#these have not been bounded to the same bounds as the model and that is intentional to force errors
prosail.runs.test <- 5000

#param.maxmin <- matrix(c(0.8, 2.5, #leaf layers or leaf structure
#                         0.2,77, #Cab
#                         0,15, #Car
#                         0.0043,0.0753, #Cw
#                         0.0017,0.0331, #Cm
#                         0,8), #LAI
#                       nrow=6,ncol = 2,byrow = T)

param.maxmin.mean <- (param.maxmin[,2]-param.maxmin[,1])/2 + param.maxmin[,1]
param.maxmin.mean
param.maxmin


#this is giving waaaay too many points outside the training range..
df.param.list.test <- data.frame(N=rtnorm(prosail.runs.test,
                                         mean=param.maxmin.mean[1],
                                         sd=param.maxmin.mean[1]/2, #this is a random criteria..
                                         lower=param.maxmin[1,1],
                                         upper=param.maxmin[1,2]),
                                 Cab=rtnorm(prosail.runs.test,
                                            mean=param.maxmin.mean[2],
                                            sd=param.maxmin.mean[2]/2, #this is a random criteria..
                                            lower=param.maxmin[2,1],
                                            upper=param.maxmin[2,2]),
                                 Car=rtnorm(prosail.runs.test,
                                            mean=param.maxmin.mean[3],
                                            sd=param.maxmin.mean[3]/2, #this is a random criteria..
                                            lower=param.maxmin[3,1],
                                            upper=param.maxmin[3,2]),
                                 Cw=rtnorm(prosail.runs.test,
                                           mean=param.maxmin.mean[4],
                                           sd=param.maxmin.mean[4]/2, #this is a random criteria..
                                           lower=param.maxmin[4,1],
                                           upper=param.maxmin[4,2]),
                                 Cm=rtnorm(prosail.runs.test,mean=param.maxmin.mean[4],
                                           sd=param.maxmin.mean[5]/2, #this is a random criteria..
                                           lower=param.maxmin[5,1],
                                           upper=param.maxmin[5,2]),
                                 LAI=rtnorm(prosail.runs.test,mean=param.maxmin.mean[5],
                                            sd=param.maxmin.mean[6]/2, #this is a random criteria..
                                            lower=param.maxmin[6,1],
                                            upper=param.maxmin[6,2]),
                                 TypeLidf="Erectophile")




param.maxmin
param.maxmin.mean
head(df.param.list.test)

summary(df.param.list.test)
summary(LHS)

#something is wrong when using the rtnorm distribution for genreating random samples..
#LHS.2 <- Latinhyper(param.maxmin,prosail.runs.test)
#first we must generate a new sampling of the trait space
#param.list.2 <- data.frame(#C1=LHS[,1],
#  C2=LHS.2[,1],
#  C3=LHS.2[,2],C4=LHS.2[,3],
#  C5=LHS.2[,4],C6=LHS.2[,5],
#  TypeLidf="Erectophile") #this might require change if you change something before
#now the names must be changed
#colnames(param.list.2)<-c(names.traits,"TypeLidf")



#generating the responses in hyperspectral and sentinel 
test.prosail <- PROSAIL(parameterList = df.param.list.test)#,prosail.runs.test)
#using the LHS
#test.prosail <- PROSAIL(parameterList = df.param.list.2)#,prosail.runs.test)

test.prosail.s2 <- spectralResampling(test.prosail,
                                      "Sentinel2",response_function = TRUE)

par(mfrow=c(2,2))
plot(mRF.spclib,main="RF train spectra set")
plot(mRF.s2.spclib)

plot(test.prosail,main="test spectra set")
plot(test.prosail.s2)

#bringing it all together
#test.prosail.s2.spectra <-spectra(test.prosail.s2)
test.prosail.s2.spectra.subsel <- test.prosail.s2[,c(2,3,4,5,6,7,8,9,12,13)]

test.df <- cbind(df.param.list.test,as.data.frame(test.prosail.s2.spectra.subsel ))
names(test.df)
names(test.df) <- c(names.traits,"TypeLidf",
                    "B02","B03","B04",
                    "B05","B06","B07",
                    "B08","B8A","B11",
                    "B12")
#again, for now we dont care for the leaf
test.df <- test.df[,-7]


#now.. we must be able to predic the df.param.list.test
mRF.test.pred <- predict(mRF_all,newdata=test.df)
out.mRF.test.pred <- get.mv.predicted(mRF.test.pred)

#lets generate a out df
out.df <- df.param.list.test[,-7]

head(out.df)
head(out.mRF.test.pred)

out.df <- cbind(out.df,out.mRF.test.pred[,1:ncol(LHS)]) #CAREFUL THE NAMES ARE THE SAME
names(out.df)
names(out.df) <- c(names.traits,paste(names.traits,"mRF",sep = "_"))

head(out.df)

#time to plot
par(mfrow=c(2,3))
ncol(LHS)
for (i in 1:ncol(LHS)){
  print(i)
  
  plot(out.df[,i],
       out.df[,i+6],
       main=paste(names.traits[i]),
       xlab="Reference",
       ylab="Predicted")
  
}
