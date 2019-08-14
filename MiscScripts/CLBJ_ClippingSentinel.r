ibrary(raster)
library(rgdal)

#test.img <- readGDAL("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/T14SPB_20190420T170851_B02_20m.jp2")
#test.img <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/T14SPB_20190420T170851_B02_20m.jp2")
#test.img

gc()

s2.file.list <- list.files("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/",
                           full.names = T)
s2.name.list <- list.files("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302/S2A_MSIL2A_20190420T170851_N0211_R112_T14SPB_20190420T213302.SAFE/GRANULE/L2A_T14SPB_A019986_20190420T171657/IMG_DATA/R20m/")

s2.name.list
#lets load it in order
s2.file.list.ordered <- s2.file.list[c(2,3,4,5,6,7,10,8,9)]
s2.name.list.ordered <- s2.name.list[c(2,3,4,5,6,7,10,8,9)]

s2.name.list
s2.name.list.ordered 

#lets stack
s2.crop.img.stack<- stack(s2.file.list.ordered)
writeRaster(s2.full.img.stack,
            filename="D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB.tif",
            options=c("COMPRESS=LZW"),
            overwrite=TRUE)

#lets crop
library(sp)
library(maptools)
shp.aoi <- readShapePoly("D:/NEON_Data/CLBJ/qgis_shape_bound/qgis_shape_bound.shp")

s2.crop.img.stack <- crop(x = s2.full.img.stack,
                          y = shp.aoi)

shp.aoi.rst <- rasterize(shp.aoi,s2.crop.img.stack)

s2.crop.mask.img.stack <- mask(x = s2.crop.img.stack,
                               mask = shp.aoi,
                               filename="D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif",
                               options=c("COMPRESS=LZW"),
                               overwrite=TRUE) 



#lets save it to a folder - BIG FILE
