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


#for std col stats
library(matrixStats)

library(keras)

#setwd
setwd("D:/ESAPhiWeek/")
gc()
dump.fld <- "./dumps"
dir.create(dump.fld) #it gives out a warning if the folder exits
set.seed(1000)


param.maxmin <- matrix(c(#1.5, 1.9, #leaf layers or leaf structure
  10,60, #Cab
  5,25, #Car
  0.01,0.02, #Cw #original it was from [0.01 to 0.02] but i suspect saturation
  0.01,0.02, #Cm
  0.5,6),#LAI
  #0.05,0.1), #hotstop
  nrow=5,ncol = 2,byrow = T)

#creating a training space
train.n <- 500 * nrow(param.maxmin) #this represents the number of runs that prosail will be, x pts per Trait
train.LHS <- Latinhyper(param.maxmin,train.n)

#to force having the limits we will add something oo the train
train.grid <- Grid(param.maxmin,train.n)
train.LHS <- rbind(train.LHS,train.grid)

valid.n <- 2* nrow(train.LHS) #this represents the number of runs that prosail will be, x pts per Trait
valid.LHS <- Latinhyper(param.maxmin,valid.n)

#checking the stats
summary(train.LHS)
summary(valid.LHS)

#sentinel position
sun_zenith = 25.8734
obs_zenith = 5.4748
rel_azimut = 281.0692 - 5.4748


#lets build the dataset and train the model
train.trait.df <- data.frame(#N=train.LHS[,1],
  Cab=train.LHS[,1],
  Car=train.LHS[,2],
  Cw=train.LHS[,3],
  Cm=train.LHS[,4],
  LAI=train.LHS[,5],
  #hspot=train.LHS[,7],
  #TypeLidf = 0 ,
  tts = sun_zenith,
  tto = obs_zenith,
  psi = rel_azimut)

valid.trait.df <- data.frame(#N=valid.LHS[,1],
  Cab=valid.LHS[,2],
  Car=valid.LHS[,2],
  Cw=valid.LHS[,3],
  Cm=valid.LHS[,4],
  LAI=valid.LHS[,5],
  #hspot=valid.LHS[,7],
  #TypeLidf = 0 ,
  tts = sun_zenith,
  tto = obs_zenith,
  psi = rel_azimut)

#there are no NA here 
summary(train.trait.df)
summary(valid.trait.df)

#creating the spectral libraries from PROSAIL RTM
train.spclib <- PROSAIL(parameterList = train.trait.df)
train.spectr <- spectra(train.spclib)
train.spectr.df <- as.data.frame(train.spectr)

valid.spclib <- PROSAIL(parameterList = valid.trait.df)
valid.spectr <- spectra(valid.spclib)
valid.spectr.df <- as.data.frame(valid.spectr)

#checking out the training and validation spectral spaces
par(mfrow=c(1,2))
plot(train.spclib)
plot(valid.spclib)

#ressampling to sentinel
s2.train.spclib <- spectralResampling(train.spclib,"Sentinel2",response_function = T)
s2.train.spectr <- spectra(s2.train.spclib)

s2.train.spclib.20m <- s2.train.spclib[,c(2,3,4,5,6,7,9,12,13)] #only the 20m bands
s2.train.spectr.20m <- spectra(s2.train.spclib.20m)

s2.valid.spclib <- spectralResampling(valid.spclib,"Sentinel2",response_function = TRUE)
s2.valid.spectr <- spectra(s2.valid.spclib)

s2.valid.spclib.20m <- s2.valid.spclib[,c(2,3,4,5,6,7,9,12,13)]
s2.valid.spectr.20m <- spectra(s2.valid.spclib.20m)

#normalizing the variables
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

tr.trait.df <- train.trait.df[,-c(6,7,8)]
vl.trait.df <- valid.trait.df[,-c(6,7,8)]

tr.trait.df.norm <- as.data.frame(lapply(tr.trait.df,normalize))
vl.trait.df.norm <- as.data.frame(lapply(vl.trait.df,normalize))

# Random sample indexes
train_index <- sample(1:nrow(tr.trait.df.norm), 0.8 * nrow(tr.trait.df.norm))
test_index <- setdiff(1:nrow(tr.trait.df.norm), train_index)

# Build X_train, y_train, X_test, y_test
X_train <- as.matrix(maxmindf[train_index, -15])
y_train <- as.matrix(maxmindf[train_index, "sales"])

X_test <- as.matrix(maxmindf[test_index, -15])
y_test <- as.matrix(maxmindf[test_index, "sales"])

model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 12, activation = 'relu', kernel_initializer='RandomNormal', input_shape = c(6)) %>% 
  layer_dense(units = 8, activation = 'relu') %>%
  layer_dense(units = 1, activation = 'linear')
summary(model)


#testing
setwd("D:/ESAPhiWeek/kerastest/datasets")
cars<-read.csv("cars.csv")

#Max-Min Normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
maxmindf <- as.data.frame(lapply(cars, normalize))
attach(maxmindf)

# Random sample indexes
train_index <- sample(1:nrow(maxmindf), 0.8 * nrow(maxmindf))
test_index <- setdiff(1:nrow(maxmindf), train_index)
# Build X_train, y_train, X_test, y_test
X_train <- as.matrix(maxmindf[train_index, -15])
y_train <- as.matrix(maxmindf[train_index, "sales"])
X_test <- as.matrix(maxmindf[test_index, -15])
y_test <- as.matrix(maxmindf[test_index, "sales"])

model %>% compile(
  loss = 'mean_squared_error',
  optimizer = 'adam',
  metrics = c('mae')
)

history <- model %>% fit(
  X_train, y_train, 
  epochs = 150, batch_size = 50, 
  validation_split = 0.2
)

model %>% evaluate(X_test, y_test)

library(tensorflow)
install_tensorflow(version="GPU")
##################### IMPLEMENT ME

library(keras)
#https://stackoverflow.com/questions/48472646/r-keras-declare-multiple-outputs


# placeholder data
Y <- data.frame(y1=1:100,y2=1:100)
X <- data.frame(x1=1:100,x2=1:00,x3=1:100)

# add covariates
input <- layer_input(shape=dim(X)[2],name="covars")

# add hidden layers
base_model <- input  %>%
  layer_dense(units = 3, activation='relu') %>%
  layer_dense(units = 2, activation='relu') 

# add outputs
y1 <- base_model %>% 
  layer_dense(units = 1, name="y1") 

y2 <- base_model %>% 
  layer_dense(units = 1, name="y2") 

# combine
model <- keras_model(input,list(y1,y2))


