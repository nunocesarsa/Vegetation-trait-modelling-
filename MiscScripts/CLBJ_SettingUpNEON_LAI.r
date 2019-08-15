
library(raster)
library(rgdal)
library(maptools)

gc()
#you could even set up an unzip, process and delete procedure but cba for that now

#path2LAI
path2LAI <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_lai-spectrometer-mosaic/",
                       pattern=".tif",
                       full.names = T)

name2LAI <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_lai-spectrometer-mosaic/",
                       pattern=".tif")

#now we first just make a big mosaic of all the data and for that, we loop.

#outfld
out.fld <- "D:/NEON_Data/CLBJ/Processed/"

#without ressampling it simply takes too much time, we need to ressample first
tempdir() #go check if there are forgotten processes here and clean it

s2.AOI.rst <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif",
                     band=1)

#folder to store all resampled data
out.fld.res <- "D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m/"


for (i in 1:length(path2LAI)){
 #focal(raster(i),
 #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2LAI[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2LAI[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
}

#now that we ressampled, lets aggreagte them
path2LAI_agg <- list.files(out.fld.res,
                           pattern=".tif",
                           full.names = T)
path2LAI_agg

temp.rst <- raster(path2LAI_agg[1])
for (i in 2:length(path2LAI_agg)){
  print(paste("processing image:", i,"of",length(path2LAI_agg)))
  
  temp.rst <- mosaic(temp.rst,raster(path2LAI_agg[i]),fun="max")
}

plot(temp.rst)
writeRaster(temp.rst,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_LAI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)

