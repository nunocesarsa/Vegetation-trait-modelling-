library(raster)
library(rgdal)
library(maptools)


s2.AOI.rst <- raster("D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_20190420_T14SPB_AOI.tif",
                     band=1)

s2.AOI.rst

s2.AOI.20m.shp <- rasterToPolygons(s2.AOI.rst)
#its a big file!
writePolyShape(s2.AOI.20m.shp,
               "D:/NEON_Data/CLBJ/S2A_MSIL2A_20190420/SL2A_T14SPB_AOI_F20mGrid.shp")


#Clipping the US land cover map (30m resolution, CONUS)

list.files("D:/NEON_Data/CLBJ/NLCD_subsection/")
LC.large.rst <- raster("D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_2016_Land_Cover_L48_20190424_sDHeFWVTqu1fNbMojiDY.tiff")

#first we ressample to 20x20m as in the sentinel image
CRS(projection(s2.AOI.rst))
LC.large.rst.UTM14N <- projectRaster(from=LC.large.rst,
                                     crs=CRS(projection(s2.AOI.rst)),
                                     method="ngb",
                                     res=30)
LC.large.rst.UTM14N
LC.large.rst
                               
LC.large.rst.20m <- resample(x=LC.large.rst.UTM14N,
                             y=s2.AOI.rst,
                             method="ngb")
#now we crop it and mask it
LC.large.rst.20m.crop.msk <- mask(crop(x=LC.large.rst.20m,
                                       y=s2.AOI.rst),
                                  mask=s2.AOI.rst)

LC.large.rst.20m.crop.msk
s2.AOI.rst
writeRaster(LC.large.rst.20m.crop.msk,
            "D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_2016_Land_Cover_AOI.tif",
            options=c("COMPRESS=LZW"),
            overwrite=TRUE)


#lets seelect the grasslands
names(LC.large.rst.20m.crop.msk) <- "ClassNr"
NLCD_grass <- LC.large.rst.20m.crop.msk[LC.large.rst.20m.crop.msk$ClassNr == 71]
plot(LC.large.rst.20m.crop.msk)
plot(LC.large.rst.20m.crop.msk==71)

NLCD_grass <- LC.large.rst.20m.crop.msk==71
plot(NLCD_grass)
NLCD_grass[NLCD_grass!=1]<- NA
writeRaster(NLCD_grass,
            "D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI.tif",
            options=c("COMPRESS=LZW"),
            overwrite=TRUE)

#converting the grass herbs layer to a point and polygon shapefile

NLCD_grass.poly.shp <- rasterToPolygons(NLCD_grass)
NLCD_grass.pts.shp <- rasterToPoints(NLCD_grass,spatial = T)

#its a big file!
nrow(NLCD_grass.poly.shp)
head(NLCD_grass.poly.shp)
writePolyShape(NLCD_grass.poly.shp,
               "D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_Poly.shp")
writePointsShape(NLCD_grass.pts.shp,
               "D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_pts.shp")

#generting a random sample of points (lets make it 5000)
row.index <- c(1:nrow(NLCD_grass.pts.shp))
head(row.index)
row.index.sp <- sample(row.index,5000)
head(row.index.sp)

NLCD_grass.pts.shp_SAMPLE <- NLCD_grass.pts.shp[row.index.sp,]
writePointsShape(NLCD_grass.pts.shp_SAMPLE,
                 "D:/NEON_Data/CLBJ/NLCD_subsection/NLCD_grass_AOI_pts_5kSample.shp")
