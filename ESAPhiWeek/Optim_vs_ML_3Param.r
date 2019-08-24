#we convert it to a spectral object
spclib.plist.300  <-PROSAIL(parameterList = train.par.df.300 )
spclib.plist.900  <-PROSAIL(parameterList = train.par.df.900 )
spclib.plist.1500 <-PROSAIL(parameterList = train.par.df.1500 )

spectra.300  <- as.data.frame(spectra(spclib.plist.300))
spectra.900  <- as.data.frame(spectra(spclib.plist.900))
spectra.1500 <- as.data.frame(spectra(spclib.plist.1500))


#lets combine all the data
df.300  <-  cbind(train.par.df.300,spectra.300)
df.900  <-  cbind(train.par.df.900,spectra.900)
df.1500 <-  cbind(train.par.df.1500,spectra.1500)

#we are training with hyperspectral data
mRF.300 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.300,block.size = 10)
mRF.900 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.900,block.size = 10)
mRF.1500 <- rfsrc(Multivar(Cab,Cw,LAI)~.,data=df.1500,block.size = 10)

#we will train the best with sentinel data.. or we can even try to predict the same traits when trained with
#sentinel data

#now lets generate a spectra of our original traind ata
df.param.list.spclib <- PROSAIL(parameterList=df.param.list)
df.param.list.spectr <- as.data.frame(spectra(df.param.list.spclib))

#we can now try to predict based on each model
mRF.pred.300  <- predict(mRF.300,newdata=df.param.list.spectr )
mRF.pred.900  <- predict(mRF.900,newdata=df.param.list.spectr )
mRF.pred.1500 <- predict(mRF.1500,newdata=df.param.list.spectr )

#now we have to extrac the tables from the predicted model object
mv.mRF.pred.300 <- get.mv.predicted(mRF.pred.300)
mv.mRF.pred.900 <- get.mv.predicted(mRF.pred.900)
mv.mRF.pred.1500 <- get.mv.predicted(mRF.pred.1500)


#and finally we can bring it to the original df and add it to the plot
mv.mRF.pred.300
#just to keep a table of the previous results
out.df.stored <- out.df

out.df$mRF_Cab300 <- mv.mRF.pred.300[,1]
out.df$mRF_Cw300  <- mv.mRF.pred.300[,2]
out.df$mRF_LAI300 <- mv.mRF.pred.300[,3]
out.df$mRF_Cab900 <- mv.mRF.pred.900[,1]
out.df$mRF_Cw900  <- mv.mRF.pred.900[,2]
out.df$mRF_LAI900 <- mv.mRF.pred.900[,3]
out.df$mRF_Cab1500 <- mv.mRF.pred.1500[,1]
out.df$mRF_Cw1500  <- mv.mRF.pred.1500[,2]
out.df$mRF_LAI1500 <- mv.mRF.pred.1500[,3]


par(mfrow=c(2,4))
plot(out.df$Cab,out.df$Cab,pch=19,xlab="Reference",ylab="Prediction",main="Cab",cex=1.5)
points(out.df$Cab,out.df$SCE_Cab,pch=2,col="blue",cex=.75)
points(out.df$Cab,out.df$LBFGSB_Cab,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_Cab~out.df$Cab),col="blue")
abline(lm(out.df$LBFGSB_Cab~out.df$Cab),col="red")

plot(out.df$Cw,out.df$Cw,pch=19,xlab="Reference",ylab="Prediction",main="Cw",cex=1.5)
points(out.df$Cw,out.df$SCE_Cw,pch=2,col="blue",cex=.75)
points(out.df$Cw,out.df$LBFGSB_Cw,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_Cw~out.df$Cw),col="blue")
abline(lm(out.df$LBFGSB_Cw~out.df$Cw),col="red")

plot(out.df$LAI,out.df$SCE_LAI,pch=19,xlab="Reference",ylab="Prediction",main="LAI",cex=1.5)
points(out.df$LAI,out.df$SCE_LAI,pch=2,col="blue",cex=.75)
points(out.df$LAI,out.df$LBFGSB_LAI,pch=3,col="red",cex=.75)
abline(lm(out.df$SCE_LAI~out.df$LAI),col="blue")
abline(lm(out.df$LBFGSB_LAI~out.df$LAI),col="red")

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','SCE','L-BFGS-B'), 
       pch=c(19,2,3), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red'))


plot(out.df$Cab,out.df$Cab,pch=19,xlab="Reference",ylab="Prediction",main="Cab",cex=1.5)
points(out.df$Cab,out.df$mRF_Cab300,pch=13,col="blue",cex=.75)
points(out.df$Cab,out.df$mRF_Cab900,pch=15,col="red",cex=.75)
points(out.df$Cab,out.df$mRF_Cab1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_Cab300~out.df$Cab),col="blue")
abline(lm(out.df$mRF_Cab900~out.df$Cab),col="red")
abline(lm(out.df$mRF_Cab1500~out.df$Cab),col="green")

plot(out.df$Cw,out.df$Cw,pch=19,xlab="Reference",ylab="Prediction",main="Cw",cex=1.5)
points(out.df$Cw,out.df$mRF_Cw300,pch=13,col="blue",cex=.75)
points(out.df$Cw,out.df$mRF_Cw900,pch=15,col="red",cex=.75)
points(out.df$Cw,out.df$mRF_Cw1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_Cw300~out.df$Cw),col="blue")
abline(lm(out.df$mRF_Cw900~out.df$Cw),col="red")
abline(lm(out.df$mRF_Cw1500~out.df$Cw),col="green")

plot(out.df$LAI,out.df$SCE_LAI,pch=19,xlab="Reference",ylab="Prediction",main="LAI",cex=1.5)
points(out.df$LAI,out.df$mRF_LAI300,pch=13,col="blue",cex=.75)
points(out.df$LAI,out.df$mRF_LAI900,pch=15,col="red",cex=.75)
points(out.df$LAI,out.df$mRF_LAI1500,pch=23,col="green",cex=.75)
abline(lm(out.df$mRF_LAI300~out.df$LAI),col="blue")
abline(lm(out.df$mRF_LAI900~out.df$LAI),col="red")
abline(lm(out.df$mRF_LAI1500~out.df$LAI),col="green")

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Target value','mRF (n=300)','mRF (n=900)','mRF (n=1500)'), 
       pch=c(19,13,15,23), pt.cex=1.5, cex=1.5, bty='n',
       col = c('black','blue','red','green'))
