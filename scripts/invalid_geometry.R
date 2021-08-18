# Fixing invalid geometries

## 

library(rgdal)
library(rgeos)
library(rJava)
options(java.parameters = "-Xmx20g" ) #asking for 23gb mem for jobs
library(dismo)
library(dplyr) 
library(readr)
#library(janitor) #need this?
library(stringr)
library(sf)

rasterOptions(tmpdir = 'temp/', progress = "text", maxmemory = 12e9) #edit this


library(sf)

#load in all bioclim and occ data

bioclim <- brick("data/bioclim_5km.tif")
# need to add layer names, not preserved when raster was written
names(bioclim) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14",
                    "bio15", "bio16", "bio17", "bio18", "bio19", "bio2",
                    "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")


load("data/plant.occ.final.RData") # this loads object 'plant.occ'

specs <- read.csv("data/species_list.csv") %>% pull(species)

buffer_dist <- read.csv("data/buffer_dist_model.csv")



plant_sf <- st_as_sf(plant.occ, coords = c("decimallongitude", "decimallatitude"),
                       crs = 4326, dim = "XY")

#convert to projected coordinates (this fixed invalid geometries)
plant_proj <- st_transform(plant_sf, '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')



occ <- subset(plant_proj, plant_proj$species == specs[i])

occ.sp <- as(occ, "Spatial")



# make range polygon
dist <- as.numeric(buffer_dist[buffer_dist$species == specs[i], "distance_km"])


range <- buffer(occ.sp, width = (dist*1000), dissolve = TRUE) #uses meters

#clip predictors to range polygon

## need to project raster first, save this to use in models

bioclim.prj <- projectRaster(bioclim, crs = crs(plant_proj))
writeRaster(bioclim.prj, filename = "D:/endophyteDimensions/data/bioclim_prj_5km.tif", options="INTERLEAVE=BAND")

preds <- crop(bioclim, extent(range)) %>% mask(range) #this all works now!


#ISSUE raster calculations do not work in projected coordinates, need to be lat/long. 

## try making range polygon then reprojecting

occ.sp <-SpatialPointsDataFrame(plant.occ[,c('decimallongitude','decimallatitude'),], 
                                data = plant.occ) 
proj4string(occ.sp) <- CRS("+proj=longlat + ellps=WGS84")


# now just project for range polygon creation, then reproject polygon for predictor crop
occ.planar <- spTransform(occ.sp, CRSobj = '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')



occ1 <- subset(occ.sp, occ.sp$species == specs[23])


occ2 <- subset(occ.planar, occ.planar$species == specs[23])


# make range polygon
dist <- as.numeric(buffer_dist[buffer_dist$species == specs[i], "distance_km"])


range1 <- buffer(occ1, width = (dist*1000), dissolve = TRUE) #uses meters

range2 <- buffer(occ2, width = (dist*1000), dissolve = TRUE)

#reproject range2
range2_sp <- spTransform(range2, crs(occ.sp)) #this lookds weird...


#try reprojection raster object? 

map_outputs <- brick("data/Acacia mearnsii_map_outputs.tif")

crs(map_outputs) #in planar coordinates


map_outputs_prj <- projectRaster(map_outputs, crs = crs(occ.sp))


# test extent and suitable area calc
