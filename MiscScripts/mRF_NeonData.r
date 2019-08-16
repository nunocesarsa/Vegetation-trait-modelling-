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

#investigate
library(randomForestSRC)

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

#these parameters com from here
#https://www.sciencedirect.com/science/article/pii/S0303243415000100
#also of interest: https://www.mdpi.com/2072-4292/11/13/1597/htm
#car limits taken from: #limits taken from: https://www.mdpi.com/2072-4292/8/2/119

param.maxmin <- matrix(c(1.5, 1.9, #leaf layers or leaf structure
                         15,55, #Cab
                         0,15, #Car
                         0.01,0.02, #Cw
                         0.005,0.01, #Cm
                         0.1,8,#LAI
                         0.05,0.1), #hotstop
                       nrow=7,ncol = 2,byrow = T)

#we can now generate thousands of sampling, lets do 500 per trait: 
prosail.runs <- nrow(param.maxmin)*1000
LHS <- Latinhyper(param.maxmin,prosail.runs)

#now lets traing a mRF on all these parameters at once and see if we are able to do good predictions

#from snap s2
#sun_zenith mean: 25.8734
#sun_azimuth mean: 141.9867
#view_zenith mean: 5.4748
#view_azimuth mean: 281.0692

sun_zenith = 25.8734
obs_zenith = 5.4748
rel_azimut = 281.0692 - 5.4748


#first we must generate a new sampling of the trait space
param.list <- data.frame(C1=LHS[,1],C2=LHS[,2],
                         C3=LHS[,3],C4=LHS[,4],
                         C5=LHS[,5],C6=LHS[,6],
                         C7=LHS[,7],
                         TypeLidf="Erectophile",
                         tts=sun_zenith, #fixed zenith, azimuth and rel azimuth
                         tto=obs_zenith,
                         psi=rel_azimut,
                         psoil=1,  #psoil - fixed this as the average that i saw from the paper [0.5 to 1.5]
                         lidfa=65) #fixed also the avreage leaf angle
#this might require change if you change something before
#now the names must be changed

names.traits <- c("N","Cab",
                   "Car","Cw",
                   "Cm","LAI","hspot")
c(names.traits,names(param.list)[8:13])


colnames(param.list)<-c(names.traits,names(param.list)[8:13])

mRF.spclib <- PROSAIL(parameterList = param.list)
mRF.s2.spclib <- spectralResampling(mRF.spclib,
                                    "Sentinel2",response_function = TRUE)
plot(mRF.spclib)
plot(mRF.s2.spclib)

#lets pick only the bands that are on our dataset
#b 2,3,4,5,6,7,8a,11,12
mRF.s2.spc.sel <- mRF.s2.spclib[,c(2,3,4,5,6,7,9,12,13)]

head(mRF.s2.spc.sel)

#brings it all together
full.df <- cbind(param.list,as.data.frame(mRF.s2.spc.sel))
names(full.df)
#next step is optional

names(full.df) <- c(c(names.traits,names(param.list)[8:13]),
                    "B02","B03","B04",
                    "B05","B06","B07",
                    "B8A","B11","B12")

names(full.df)
#before going into the model we should remove  everything that is constant... since we DONT want to predict that
#just one type of leaf
full.df <- full.df[,-c(8:13)]
names(full.df)


#training the mfr
mRF_all <- rfsrc(Multivar(N,Cab,Car,
                          Cw,Cm,LAI,hspot)~.,data=full.df,block.size = 50)


#now, we will get the samples from the point shapefile we created elsewhere - this case only grasslands
neon.test.shp <- readShapePoints("D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_small_pts_5kSample.shp")

#now lets extract the values from the sentinel image
s2.AOI.stack <- stack("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif")/10000 #to convert to reflectance

#the names of the layers come in different order, for simplicity lets make it all the same
names(s2.AOI.stack)<-c("B02","B03","B04",
                       "B05","B06","B07",
                       "B8A","B11","B12")

s2.AOI.refl.shp <- raster::extract(s2.AOI.stack,
                                   neon.test.shp,
                                   sp=T)

#now lets load all the data from the traits of NEON and have that also in our shape
#loading the raster
path2neon <- list.files("D:/NEON_Data/CLBJ/Processed/",
                        pattern=".tif",
                        full.names=T)
name2neon <- list.files("D:/NEON_Data/CLBJ/Processed/",
                        pattern=".tif")

#now we remove the one that doesn't matter.. BEWARE - if you change anything .tif in th folder, this might change
path2neon <- path2neon[-1]
name2neon <- name2neon[-1]
path2neon
name2neon

#they also need to be cropped to the same extents or we cant stack them
AOI_small.shp <- readShapePoly("D:/NEON_Data/CLBJ/qgis_shape_bound/qgis_shape_bound_small_Area_diss.shp")
for (i in 1:length(path2neon)){
  
  out.file <- paste("D:/NEON_Data/CLBJ/Processed_cropped/Small_",name2neon[i],sep="")
  print(out.file)
  crop(x=raster(path2neon[i]), y=AOI_small.shp, 
       filename=out.file,
       options=c("COMPRESS=LZW"),
       overwrite=TRUE)
}

#loading the raster
path2neon.small <- list.files("D:/NEON_Data/CLBJ/Processed_cropped/",
                        pattern=".tif",
                        full.names=T)
name2neon.small <- list.files("D:/NEON_Data/CLBJ/Processed_cropped/",
                        pattern=".tif")

#now we can stack and extract everything in one go
valid.shp <- raster::extract(stack(path2neon.small ),
                             s2.AOI.refl.shp,
                             sp=T)

head(valid.shp)
#now we want proper names... these names are so freaking bad
names(valid.shp) <- c(names(valid.shp)[1:12],
                      c("biomass","chm","fpar",
                        "lai","msi","ndii",
                        "ndwi","ndmi","wbi"))
head(valid.shp)

#no we have everything to make some prediction, we can go use our trained mRF to predict some outputs
valid.df <- as.data.frame(valid.shp)
head(valid.df)

#only the bands to avoid confusion
valid.df.bandsOnly <- valid.df[,c(4:12)]
head(valid.df.bandsOnly)

#now lets do some predictions
mRF.test.pred <- predict(mRF_all,newdata=valid.df.bandsOnly)
out.mRF.test.pred <- get.mv.predicted(mRF.test.pred)

#lets wrap it up in a dataframe
test.df <- as.data.frame(out.mRF.test.pred)
names(test.df)
valid.df.neonOnly <- valid.df[,c(13:21)]
names(valid.df.neonOnly)
test.df <- cbind(test.df,valid.df.neonOnly)

#using the function of https://www.mdpi.com/2072-4292/11/13/1597/htm - AGB as a function of LAI*cm
par(mfrow=c(1,1))
plot(test.df$lai,test.df$LAI,xlab="NEON data",ylab="Prediction (mRF): LAI")

#canopy water content indices
par(mfrow=c(3,2))

plot(test.df$msi,test.df$Cw,xlab="NEON data: MSI",ylab="Prediction (mRF): Cw")
plot(test.df$ndii,test.df$Cw,xlab="NEON data: ndii",ylab="Prediction (mRF): Cw")
plot(test.df$ndwi,test.df$Cw,xlab="NEON data: ndwi",ylab="Prediction (mRF): Cw")
plot(test.df$ndmi,test.df$Cw,xlab="NEON data: ndmi",ylab="Prediction (mRF): Cw")
plot(test.df$wbi,test.df$Cw,xlab="NEON data: wbi",ylab="Prediction (mRF): Cw")

#using the function of https://www.mdpi.com/2072-4292/11/13/1597/htm - AGB as a function of LAI*cm
par(mfrow=c(1,1))
plot(test.df$biomass,(test.df$Cm*test.df$LAI),xlab="NEON data: biomass",ylab="Prediction (mRF): Cm*LAI")

#testing to fpar
plot(test.df$fpar,(test.df$LAI),xlab="NEON data: fpar",ylab="Prediction (mRF): LAI")

#check this
names(test.df)
paste("mRF_",names(test.df)[1:7],sep="")
paste("NEON_",names(test.df)[8:16],sep="")
names(test.df) <- c(paste("mRF_",names(test.df)[1:7],sep=""),
                    paste("NEON_",names(test.df)[8:16],sep=""))
names(test.df)

library(corrplot)
M<-cor(test.df)
corrplot(M, type="upper")
