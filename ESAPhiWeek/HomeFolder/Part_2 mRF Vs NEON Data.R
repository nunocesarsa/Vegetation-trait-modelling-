
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

#GIS/RS
library(maptools)
library(rgeos)
library(rgdal)
library(raster)
library(sp)

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
                         0.005,0.03, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
                         0.005,0.03, #Cm
                         0.1,8,#LAI
                         0.05,0.1), #hotstop
                       nrow=7,ncol = 2,byrow = T)

#we can now generate thousands of sampling, lets do 500 per trait: 
prosail.runs <- nrow(param.maxmin)*1500
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
names(valid.shp)
names(valid.shp)[1:12]
#now we want proper names... these names are so freaking bad
names(valid.shp) <- c(names(valid.shp)[1:12],
                      'ARVI','BIOm','CHM',
                      'EVI','fPAR','LAI',
                      'MSI','NDII','NDLI',
                      'NDNI','NDVI','NDWI',
                      'NMDI','PRI','SAVI',
                      'WBI')
                      #c("biomass","chm","fpar",
                      #  "lai","msi","ndii",
                      #  "ndwi","ndmi","wbi"))
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


##recreate the reflectances to see if we can calculate the indices
new.param.list <- as.data.frame(out.mRF.test.pred )
head(new.param.list)
new.param.list$TypeLidf="Erectophile"
new.param.list$tts=sun_zenith #fixed zenith, azimuth and rel azimuth
new.param.list$tto=obs_zenith
new.param.list$psi=rel_azimut
new.param.list$psoil=1  #psoil - fixed this as the average that i saw from the paper [0.5 to 1.5]
new.param.list$lidfa=65
head(new.param.list)
pred.prosail.scplib <- PROSAIL(parameterList = new.param.list)

par(mfrow=c(1,1))
plot(pred.prosail.scplib)

#now we create a dataframe to receive everything
pred.spectra <- as.data.frame(spectra(pred.prosail.scplib))
names(pred.spectra) <- paste("band_",pred.prosail.scplib@wavelength,sep="")
names(pred.spectra)

#now we calculate the indices
pred.indices <- data.frame(linenr=seq(1:nrow(pred.spectra)))
#sources: https://www.harrisgeospatial.com/docs/CanopyWaterContent.html
pred.indices$MSI  <- pred.spectra$band_1599/pred.spectra$band_819
pred.indices$NDII <- (pred.spectra$band_819-pred.spectra$band_1649)/(pred.spectra$band_819+pred.spectra$band_1649)
pred.indices$NDWI <- (pred.spectra$band_857-pred.spectra$band_1241)/(pred.spectra$band_857+pred.spectra$band_1241)
pred.indices$NDMI <- (pred.spectra$band_860-(pred.spectra$band_1640-pred.spectra$band_2130))/(pred.spectra$band_860+(pred.spectra$band_1640-pred.spectra$band_2130))
pred.indices$WBI  <- pred.spectra$band_970/pred.spectra$band_900

#canopy water content indices
par(mfrow=c(3,2))
plot(test.df$NEON_msi ,pred.indices$MSI,xlab="NEON data: MSI",ylab="Predicted MSI")
plot(test.df$NEON_ndii,pred.indices$NDII,xlab="NEON data: ndii",ylab="Predicted NDII")
plot(test.df$NEON_ndwi,pred.indices$NDWI,xlab="NEON data: ndwi",ylab="Predicted NDWI")
plot(test.df$NEON_ndmi,pred.indices$NDMI,xlab="NEON data: ndmi",ylab="Predicted NDMI")
plot(test.df$NEON_wbi ,pred.indices$WBI,xlab="NEON data: wbi",ylab="PredictedWBI")


#can it be an effect of outliers?
par(mfrow=c(1,1))

ncol(test.df)
par(mfrow=c(3,3))
for (i in 8:ncol(test.df)){
  boxplot(test.df[,i])
}


test.df.redux <- test.df
for (i in 8:ncol(test.df)){
  outlier.redux <- boxplot(test.df.redux[,i], plot=FALSE)$out
  print(length(outlier.redux))
  
  if (length(outlier.redux)!=0){
    test.df.redux<-test.df.redux[-which(test.df.redux[,i] %in% outlier.redux),]
  }
  boxplot(test.df.redux[,i],main=names(test.df.redux)[i])
}

#so now, i removed all the outliers from the validation data. How many points i have left?
nrow(test.df.redux)
#still over 4k, lets repeat the previous
par(mfrow=c(1,1))
plot(test.df.redux$NEON_lai,test.df.redux$mRF_LAI,xlab="NEON data (no outliers)",ylab="Prediction (mRF): LAI")

#canopy water content indices
par(mfrow=c(3,2))

plot(test.df.redux$NEON_msi ,test.df.redux$mRF_Cw,xlab="NEON data: MSI",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NEON_ndii,test.df.redux$mRF_Cw,xlab="NEON data: ndii",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NEON_ndwi,test.df.redux$mRF_Cw,xlab="NEON data: ndwi",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NEON_ndmi,test.df.redux$mRF_Cw,xlab="NEON data: ndmi",ylab="Prediction (mRF): Cw")
plot(test.df.redux$NEON_wbi ,test.df.redux$mRF_Cw,xlab="NEON data: wbi",ylab="Prediction (mRF): Cw")


##recreate the reflectances to see if we can calculate the indices - redux dataset
names(test.df.redux)
new.param.list.redux <- test.df.redux[,1:7]
head(new.param.list.redux)
names(new.param.list.redux) <- names(as.data.frame(out.mRF.test.pred ))
head(new.param.list.redux)

new.param.list.redux$TypeLidf="Erectophile"
new.param.list.redux$tts=sun_zenith #fixed zenith, azimuth and rel azimuth
new.param.list.redux$tto=obs_zenith
new.param.list.redux$psi=rel_azimut
new.param.list.redux$psoil=1  #psoil - fixed this as the average that i saw from the paper [0.5 to 1.5]
new.param.list.redux$lidfa=65
head(new.param.list.redux)
pred.prosail.scplib.redux <- PROSAIL(parameterList = new.param.list.redux)


par(mfrow=c(1,1))
plot(pred.prosail.scplib.redux)

#now we create a dataframe to receive everything
pred.spectra.redux <- as.data.frame(spectra(pred.prosail.scplib.redux))
names(pred.spectra.redux) <- paste("band_",pred.prosail.scplib.redux@wavelength,sep="")
names(pred.spectra.redux)

#now we calculate the indices
pred.indices.redux <- data.frame(linenr=seq(1:nrow(pred.spectra.redux)))
#sources: https://www.harrisgeospatial.com/docs/CanopyWaterContent.html
pred.indices.redux$MSI  <- pred.spectra.redux$band_1599/pred.spectra.redux$band_819
pred.indices.redux$NDII <- (pred.spectra.redux$band_819-pred.spectra.redux$band_1649)/(pred.spectra.redux$band_819+pred.spectra.redux$band_1649)
pred.indices.redux$NDWI <- (pred.spectra.redux$band_857-pred.spectra.redux$band_1241)/(pred.spectra.redux$band_857+pred.spectra.redux$band_1241)
pred.indices.redux$NDMI <- (pred.spectra.redux$band_860-(pred.spectra.redux$band_1640-pred.spectra.redux$band_2130))/(pred.spectra.redux$band_860+(pred.spectra.redux$band_1640-pred.spectra.redux$band_2130))
pred.indices.redux$WBI  <- pred.spectra.redux$band_970/pred.spectra.redux$band_900

#canopy water content indices
par(mfrow=c(3,2))
plot(test.df.redux$NEON_msi ,pred.indices.redux$MSI,xlab="NEON data: MSI",ylab="Predicted MSI")
plot(test.df.redux$NEON_ndii,pred.indices.redux$NDII,xlab="NEON data: ndii",ylab="Predicted NDII")
plot(test.df.redux$NEON_ndwi,pred.indices.redux$NDWI,xlab="NEON data: ndwi",ylab="Predicted NDWI")
plot(test.df.redux$NEON_ndmi,pred.indices.redux$NDMI,xlab="NEON data: ndmi",ylab="Predicted NDMI")
plot(test.df.redux$NEON_wbi ,pred.indices.redux$WBI,xlab="NEON data: wbi",ylab="PredictedWBI")

#test r squared
summary(lm(pred.indices.redux$MSI~test.df.redux$NEON_msi))$r.squared
names(test.df.redux)
names(pred.indices.redux)
redux.indices.table <- cbind(test.df.redux[,12:16],pred.indices.redux[,c(2:6)])
head(redux.indices.table)

library(corrplot)
par(mfrow=c(1,1))
M.indices.redux <-cor(redux.indices.table)
corrplot(M.indices.redux, type="upper")

#lets bring the unreduced
nrow(test.df)
nrow(pred.indices)
indices.table <- cbind(test.df[,c(12:16)],pred.indices[,c(2:6)])
names(indices.table)
M.indices <-cor(indices.table)
corrplot(M.indices, type="upper")

par(mfrow=c(1,2))
corrplot(M.indices, type="upper",title="With outliers",mar=c(0,0,1,0),addCoef.col = "black",diag=FALSE)
corrplot(M.indices.redux, type="upper",main="Without outliers",mar=c(0,0,1,0),addCoef.col = "black",diag=FALSE)
