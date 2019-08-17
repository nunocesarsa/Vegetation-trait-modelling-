library(raster)
library(rgdal)
library(maptools)

gc()
#you could even set up an unzip, process and delete procedure but cba for that now

#path2biomass
path2XAN <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_xanthophyll-canopy-mosaic/",
                        pattern=".tif",
                        full.names = T)

name2XAN <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_xanthophyll-canopy-mosaic/",
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
temp.rst <- raster(path2XAN[1])
temp.rst

for (i in 1:length(path2XAN)){
  #focal(raster(i),
  #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2XAN[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2XAN[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
}

#now we have multiple indices, we have to run each independant
path2XAN.ARVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="ARVI.tif",
                           full.names = T)

name2XAN.ARVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="ARVI.tif")

path2XAN.EVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="EVI.tif",
                            full.names = T)

name2XAN.EVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="EVI.tif")

path2XAN.NDLI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="NDLI.tif",
                           full.names = T)

name2XAN.NDLI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="NDLI.tif")

path2XAN.NDNI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDNI.tif",
                            full.names = T)

name2XAN.NDNI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDNI.tif")

path2XAN.NDVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDVI.tif",
                            full.names = T)

name2XAN.NDVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDVI.tif")

path2XAN.PRI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="PRI.tif",
                            full.names = T)

name2XAN.PRI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="PRI.tif")

path2XAN.SAVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="SAVI.tif",
                           full.names = T)

name2XAN.SAVI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="SAVI.tif")


#now that we ressampled, lets aggreagte them
temp.ARVI <- raster(path2XAN.ARVI[1])
temp.EVI  <- raster(path2XAN.EVI[1])
temp.NDLI <- raster(path2XAN.NDLI[1])
temp.NDNI <- raster(path2XAN.NDNI[1])
temp.NDVI <- raster(path2XAN.NDVI[1])
temp.PRI  <- raster(path2XAN.PRI[1])
temp.SAVI <- raster(path2XAN.SAVI[1])


for (i in 2:length(path2XAN.ARVI)){
  print(paste("processing image set:", i,"of",length(path2XAN.ARVI)))
  
  #now that we ressampled, lets aggreagte them
  temp.ARVI <- mosaic(temp.ARVI,raster(path2XAN.ARVI[i]),fun="max")
  temp.EVI  <- mosaic(temp.EVI, raster(path2XAN.EVI[i]) ,fun="max")
  temp.NDLI <- mosaic(temp.NDLI,raster(path2XAN.NDLI[i]),fun="max")
  temp.NDNI <- mosaic(temp.NDNI,raster(path2XAN.NDNI[i]),fun="max")
  temp.NDVI <- mosaic(temp.NDVI,raster(path2XAN.NDVI[i]),fun="max")
  temp.PRI  <- mosaic(temp.PRI, raster(path2XAN.PRI[i]) ,fun="max")
  temp.SAVI <- mosaic(temp.SAVI,raster(path2XAN.SAVI[i]),fun="max")
  
}


writeRaster(temp.ARVI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_ARVI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.EVI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_EVI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.NDLI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NDLI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.NDNI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NDNI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.NDVI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NDVI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.PRI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_PRI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
writeRaster(temp.SAVI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_SAVI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
