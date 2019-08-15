library(raster)
library(rgdal)
library(maptools)

gc()
#you could even set up an unzip, process and delete procedure but cba for that now

#path2biomass
path2CHM <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_struct-ecosystem/",
                        pattern=".tif",
                        full.names = T)

name2CHM <- list.files("D:/NEON_Data/CLBJ/Mosaic/NEON_struct-ecosystem/",
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
temp.rst <- raster(path2CHM[1])
temp.rst

for (i in 1:length(path2CHM)){
  #focal(raster(i),
  #     w=matrix(1,20,20)) #to aggregate
  out.name <- paste(out.fld.res,
                    "NotA_20m_",
                    name2CHM[i],
                    sep="")
  print(out.name)
  aggregate(raster(path2CHM[i]), 
            fact=20, 
            fun=mean, expand=FALSE, na.rm=TRUE, 
            filename=out.name,
            options=c("COMPRESS=LZW"),overwrite=TRUE)
  
  
}

#now that we ressampled, lets aggreagte them
path2CHM_agg <- list.files(out.fld.res,
                            pattern="CHM.tif",
                            full.names = T)
path2CHM_agg

temp.rst <- raster(path2CHM_agg[1])
for (i in 2:length(path2CHM_agg)){
  print(paste("processing image:", i,"of",length(path2CHM_agg)))
  
  temp.rst <- mosaic(temp.rst,raster(path2CHM_agg[i]),fun="max")
}

plot(temp.rst)
writeRaster(temp.rst,
            filename = paste(out.fld,
                             "Final_20m_NEON_D11_CLBJ_CHM.tif",sep=""),
            options=c("COMPRESS=LZW"),overwrite=TRUE)
