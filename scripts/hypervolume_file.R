# File for Hypervolume file creation

library(raster)
library(tidyverse)


# read in bioclim and occurrence file

bioclim <- brick("D:/endophyteDimensions/data/bioclim_5km.tif")
names(bioclim) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14",
                    "bio15", "bio16", "bio17", "bio18", "bio19", "bio2",
                    "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")


load("D:/endophyteDimensions/data/plant.occ.final.RData")

# turn points into spatial object
occ.sp <-SpatialPointsDataFrame(plant.occ[,c('decimallongitude','decimallatitude'),], 
                                data = plant.occ) 

proj4string(occ.sp) <- CRS("+init=epsg:4326")


occ.test <- occ.sp[1:10,]


occ.test$bio1 <- raster::extract(bioclim[["bio1"]], occ.test)
occ.test$bio2 <- raster::extract(bioclim[["bio2"]], occ.test)
occ.test$bio10 <- raster::extract(bioclim[["bio10"]], occ.test)


# then to make final dataset
occ.test[, c("species", "bio1","bio10")] %>% as.data.frame()


#final set, I know there is a cleaner way to do this but I don't feel like figuring it out

occ.sp$bio1 <- raster::extract(bioclim[["bio1"]], occ.sp)
occ.sp$bio2 <- raster::extract(bioclim[["bio2"]], occ.sp)
occ.sp$bio3 <- raster::extract(bioclim[["bio3"]], occ.sp)
occ.sp$bio4 <- raster::extract(bioclim[["bio4"]], occ.sp)
occ.sp$bio5 <- raster::extract(bioclim[["bio5"]], occ.sp)
occ.sp$bio6 <- raster::extract(bioclim[["bio6"]], occ.sp)
occ.sp$bio7 <- raster::extract(bioclim[["bio7"]], occ.sp)
occ.sp$bio8 <- raster::extract(bioclim[["bio8"]], occ.sp)
occ.sp$bio9 <- raster::extract(bioclim[["bio9"]], occ.sp)
occ.sp$bio10 <- raster::extract(bioclim[["bio10"]], occ.sp)
occ.sp$bio11 <- raster::extract(bioclim[["bio11"]], occ.sp)
occ.sp$bio12 <- raster::extract(bioclim[["bio12"]], occ.sp)
occ.sp$bio13 <- raster::extract(bioclim[["bio13"]], occ.sp)
occ.sp$bio14 <- raster::extract(bioclim[["bio14"]], occ.sp)
occ.sp$bio15 <- raster::extract(bioclim[["bio15"]], occ.sp)
occ.sp$bio16 <- raster::extract(bioclim[["bio16"]], occ.sp)
occ.sp$bio17 <- raster::extract(bioclim[["bio17"]], occ.sp)
occ.sp$bio18 <- raster::extract(bioclim[["bio18"]], occ.sp)
occ.sp$bio19 <- raster::extract(bioclim[["bio19"]], occ.sp)


occ.hyper <- occ.sp[, c("species", "bio1", "bio2", "bio3", "bio4",
                        "bio5", "bio6", "bio7", "bio8", "bio9", "bio10",
                        "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                        "bio17", "bio18", "bio19")] %>% as.data.frame()


write.csv(occ.hyper, "D:/endophyteDimensions/data/hypervolume_data.csv")
