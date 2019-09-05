library(raster)
library(rgdal)

#test.img <- readGDAL("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/T14SPB_20190420T170851_B02_20m.jp2")
#test.img <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/T14SPB_20190420T170851_B02_20m.jp2")
#test.img

gc()

s2.file.list <- list.files("D:/OVP_EOLDAS_R/ESA/S2B_MSIL2A_20190717T104029_N0213_R008_T31UFU_20190717T151314.SAFE/GRANULE/L2A_T31UFU_A012332_20190717T104353/IMG_DATA/R20m/",
                           full.names = T)
s2.name.list <- list.files("D:/OVP_EOLDAS_R/ESA/S2B_MSIL2A_20190717T104029_N0213_R008_T31UFU_20190717T151314.SAFE/GRANULE/L2A_T31UFU_A012332_20190717T104353/IMG_DATA/R20m")

s2.name.list
#lets load it in order
s2.file.list.ordered <- s2.file.list[c(2,3,4,5,6,7,10,8,9)]
s2.name.list.ordered <- s2.name.list[c(2,3,4,5,6,7,10,8,9)]

s2.file.list.ordered
s2.name.list.ordered 

#lets stack
s2.full.img.stack<- stack(s2.file.list.ordered)

writeRaster(s2.full.img.stack,
            filename="D:/ESAPhiWeek/ESA_data/S2B_MSIL2A_20190717_T31UFU.tif",
            options=c("COMPRESS=LZW"),
            overwrite=TRUE)


#lets crop
library(sp)
library(maptools)
shp.rst <- raster("D:/OVP_EOLDAS_R/TestData/S2_20190401_OVPMask_RGB8A_UTM31N.tif")
#shp.aoi <- readShapePoly("D:/NEON_Data/CLBJ/qgis_shape_bound/qgis_shape_bound.shp")

s2.crop.img.stack <- crop(x = s2.full.img.stack,
                          y = shp.rst )

#shp.aoi.rst <- rasterize(shp.aoi,s2.crop.img.stack)

s2.crop.mask.img.stack <- mask(x = s2.crop.img.stack,
                               mask = shp.rst,
                               filename="D:/ESAPhiWeek/ESA_data/S2B_MSIL2A_20190717_T31UFU_Small.tif",
                               options=c("COMPRESS=LZW"),
                               overwrite=TRUE) 
