library(raster)
library(rgdal)
library(maptools)

gc()
#you could even set up an unzip, process and delete procedure but cba for that now

#path2biomass
path2BIO <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_biomass-spectrometer-mosaic/",
                       pattern=".tif",
                       full.names = T)

name2BIO <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_biomass-spectrometer-mosaic/",
                       pattern=".tif")

#outfld
out.fld <- "D:/NEON_Data/CLBJ/Processed/"

#without ressampling it simply takes too much time, we need to ressample first
tempdir() #go check if there are forgotten processes here and clean it

s2.AOI.rst <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif",
                     band=1)

#folder to store all resampled data
out.fld.res <- "D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m/"

#FIRST LETS CHECK THE RESOLUTION
temp.rst <- raster(path2BIO[1])
temp.rst

for (i in 1:length(path2BIO)){
  #focal(raster(i),
  #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2BIO[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2BIO[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
}

#now that we ressampled, lets aggreagte them
path2BIO_agg <- list.files(out.fld.res,
                           pattern="Biomass.tif",
                           full.names = T)
path2BIO_agg

temp.rst <- raster(path2BIO_agg[1])
for (i in 2:length(path2BIO_agg)){
  print(paste("processing image:", i,"of",length(path2BIO_agg)))
  
  temp.rst <- mosaic(temp.rst,raster(path2BIO_agg[i]),fun="max")
}

plot(temp.rst)
writeRaster(temp.rst,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_BIO.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
