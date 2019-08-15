library(raster)
library(rgdal)
library(maptools)

gc()

#you could even set up an unzip, process and delete procedure but cba for that now

#path2LAI
path2h2o <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_conc-h2o-canopy-mosaic/WaterIndices",
                       pattern=".tif",
                       full.names = T)

name2h2o <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_conc-h2o-canopy-mosaic/WaterIndices",
                       pattern=".tif")


#now we can iterate the same as before
#outfld
out.fld <- "D:/NEON_Data/CLBJ/Processed/"
#without ressampling it simply takes too much time, we need to ressample first
tempdir() #go check if there are forgotten processes here and clean it
s2.AOI.rst <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif",
                     band=1)
#folder to store all resampled data
out.fld.res <- "D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m/"

for (i in 1:length(path2h2o)){
  #focal(raster(i),
  #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2h2o[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2h2o[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
} #no


#now we have multiple indices, we have to run each independant
path2h2o.MSI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="MSI.tif",
                           full.names = T)

name2h2o.MSI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="MSI.tif")

path2h2o.NDII <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDII.tif",
                            full.names = T)

name2h2o.NDII <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDII.tif")

path2h2o.NDWI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDWI.tif",
                            full.names = T)

name2h2o.NDWI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NDWI.tif")

path2h2o.NMDI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NMDI.tif",
                            full.names = T)

name2h2o.NMDI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                            pattern="NMDI.tif")

path2h2o.WBI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="WBI.tif",
                           full.names = T)

name2h2o.WBI <- list.files("D:/NEON_Data/CLBJ/BiophysicalParamsNEON_20m",
                           pattern="WBI.tif")


#now that we ressampled, lets aggreagte them
name2h2o.WBI

temp.MSI <- raster(path2h2o.MSI[1])
temp.NDII <- raster(path2h2o.NDII[1])
temp.NDWI <- raster(path2h2o.NDWI[1])
temp.NMDI <- raster(path2h2o.NMDI[1])
temp.WBI <- raster(path2h2o.WBI[1])


for (i in 2:length(path2h2o.MSI)){
  print(paste("processing image set:", i,"of",length(path2h2o.WBI)))
  
  temp.MSI  <- mosaic(temp.MSI,raster(path2h2o.MSI[i]),fun="max")
  temp.NDII  <- mosaic(temp.NDII,raster(path2h2o.NDII[i]),fun="max")
  temp.NDWI  <- mosaic(temp.NDWI,raster(path2h2o.NDWI[i]),fun="max")
  temp.NMDI  <- mosaic(temp.NMDI,raster(path2h2o.NMDI[i]),fun="max")
  temp.WBI  <- mosaic(temp.WBI,raster(path2h2o.WBI[i]),fun="max")
}


writeRaster(temp.MSI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_MSI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)

writeRaster(temp.NDII,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NDII.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)

writeRaster(temp.NDWI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NDWI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)

writeRaster(temp.NMDI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_NMDI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)

writeRaster(temp.WBI,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_WBI.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
