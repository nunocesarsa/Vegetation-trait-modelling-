library(raster)
library(rgdal)
library(maptools)

gc()
#you could even set up an unzip, process and delete procedure but cba for that now

#path2biomass
path2fpar <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_fpar-spectrometer-mosaic/",
                       pattern=".tif",
                       full.names = T)

name2fpar <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_fpar-spectrometer-mosaic/",
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
temp.rst <- raster(path2fpar[1])
temp.rst

for (i in 1:length(path2fpar)){
  #focal(raster(i),
  #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2fpar[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2fpar[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
}

#now that we ressampled, lets aggreagte them
path2fpar_agg <- list.files(out.fld.res,
                           pattern="fPAR.tif",
                           full.names = T)
path2fpar_agg

temp.rst <- raster(path2fpar_agg[1])
for (i in 2:length(path2fpar_agg)){
  print(paste("processing image:", i,"of",length(path2fpar_agg)))
  
  temp.rst <- mosaic(temp.rst,raster(path2fpar_agg[i]),fun="max")
}

plot(temp.rst)
writeRaster(temp.rst,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_fPAR.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)