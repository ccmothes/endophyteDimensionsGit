# calculate buffer distance for each species

# Set Up -----------------------------------------------------------
library(raster)
library(geodist)
#library(geosphere)
#library(rdist)
library(cluster)
library(tidyverse)



# Read in occurrence file

load("D:/endophyteDimensions/data/thinned_occ_all.RData")
spec_occ <- read_csv("data/species_occ_totalUpdate.csv")



# Buffer calculation ----------------------------------------

## geodist will not work on specs with > 14,000 occ (275 specs)
## first round do just those < 14,000 and remove <1 cause error

specs <- spec_occ %>% filter(num_occ < 14000 & num_occ > 1) %>% pull(species) %>% unique()

##rename columns for geodist
names(thin.occ.final.fix)[8:9] <- c("latitude", "longitude")

dist <- vector("list", length = length(specs))

for (i in 75:length(specs)){
  
  s <- dplyr::filter(thin.occ.final.fix, species == specs[i])
  
  #d <- distm(as.matrix(s[,c("longitude","latitude")]), fun = distHaversine)
  #try rdist
  #d <- pdist(as.matrix(s[,c("decimallongitude","decimallatitude")]), metric = "euclidean")
  
  #try geodist
  d <- geodist(as.matrix(s[,c("longitude","latitude")]), sequential = FALSE, measure = "haversine") 
  #this seems faster, returns exact same values

  ag <-
    agnes(d,
          diss = TRUE,
          metric = "euclidean",
          method = "single")
  dist[[i]] <- max(ag$height)/1000
  
  print(i)
  
}

dist.all <- dist[1:length(specs)] %>% unlist()

## save distances with species names
buffer_dist <- data.frame(species = specs, distance_km = dist.all)
# left join with number of occurrences and cut distances in half to make range polygons

buffer_dist <- buffer_dist %>% as_tibble() %>% left_join(spec_occ, by = "species") %>% select(-X1) %>% 
  mutate(distance_km_half = distance_km/2)

write_csv(buffer_dist, "data/buffer_distances.csv")


# see how many have really large distances

buffer_dist %>% filter(num_occ >= 10) # 7233 total to run models on

# how many are likely in US and China (~11,000 km apart)
buffer_dist %>% filter(num_occ >= 10) %>% 
  filter(distance_km > 10000) # 401

# how many likely china and Us OR Europe (~6,000 km)
buffer_dist %>% filter(num_occ >= 10) %>% 
  filter(distance_km > 6000) # 1,797

# investigate max distances
max(buffer_dist$distance_km) #17604
buffer_dist %>% filter(distance_km > 17000)


# investigate some of these on a map
thin.occ.final.fix %>% filter(species == "Cremanthodium campanulatum") %>% 
  select(longitude, latitude) %>% write.csv("data/test_buffer.csv")


##only use half the distance to make geo extent so make new column with distance/2
buffer_dist <- mutate(buffer_dist, distance_km_half = distance_km/2)
write.csv(buffer_dist, file = "buffer_distances.csv")



# INVESTIGATE buffer threshold ---------------------------------------

load("D:/endophyteDimensions/data/thinned_occ_all_update.RData")
buffer_dist <- read_csv("data/buffer_distances_final.csv")

hist(buffer_dist$distance_km)
hist(buffer_dist$distance_km_half)

## separate those >1650 and <1650
greater <- buffer_dist %>% filter(distance_km > 1650)
less <- buffer_dist %>% filter(distance_km < 1650)

## choose 15 random species to investigate in each
greater_id <- sample.int(n = nrow(greater), size = 15)
less_id <- sample.int(n = nrow(less), size = 15)


## filter species from random IDs
greater_sp <- greater[greater_id,]
less_sp <- less[less_id,]

bind_rows(greater_sp, less_sp) %>% write.csv("data/buffer_test.csv")

## UPDATE it is 1650 HALF km distance
## most still fit the above/below group, but need to redo a couple

# need 4 more above species
greater <- buffer_dist %>% filter(distance_km_half > 1650)

## choose 4 random species to investigate in each
greater_id <- sample.int(n = nrow(greater), size = 4)

## filter species from random IDs
greater_sp <- greater[greater_id,]


## now go through and vizualize all these sp. points 
## and record native/inv range
thin.occ.update %>% filter(species == paste(less_sp[3,1])) %>% 
  dplyr::select(decimallongitude, decimallatitude) %>% 
  write.table(., "clipboard", sep="\t", col.names = TRUE)

## when occ too big for clipboard
thin.occ.update %>% filter(species == paste(greater_sp[4,1])) %>% 
  dplyr::select(decimallongitude, decimallatitude) %>% 
  write.csv("data/test_points.csv")


# REDO SPECIES ABOVE THRESHOLD ---------------------------------------

# read in buffer distance file (NOTE this does not include large occ sp (see below))

buffer_dist <- read_csv("data/buffer_distances_final.csv")

## filter buffers above threshold
buffer_redo <- buffer_dist %>% filter(distance_km_half > 1650) #2715 species


##rename columns of updated file for geodist
names(thin.occ.update)[8:9] <- c("latitude", "longitude")

dist <- vector("list", length = nrow(buffer_redo))

for (i in 1:nrow(buffer_redo)){
  
  s <- dplyr::filter(thin.occ.update, species == paste(buffer_redo[i,1]))
  
  d <- geodist(as.matrix(s[,c("longitude","latitude")]), sequential = FALSE, measure = "haversine") 
  
  ag <-
    agnes(d,
          diss = TRUE,
          metric = "euclidean",
          method = "single")
  
  if((max(ag$height)/1000) > 3330){
    dist[[i]] <- max(ag$height[ag$height < (3330*1000)]/1000) 
    
  }
  else{
    dist[[i]] <- max(ag$height)/1000
    
  }
  
  print(i)
  
}

dist.redo <- dist[1:nrow(buffer_redo)] %>% unlist()

## save distances with species names
buffer_redo_dist <- data.frame(species = buffer_redo$species, distance_km = dist.redo)
# left join with number of occurrences and cut distances in half to make range polygons

buffer_redo_occ <- thin.occ.update %>% group_by(species) %>% count() %>% 
  filter(species %in% buffer_redo$species)

buffer_redo_dist <- buffer_redo_dist %>% as_tibble() %>% left_join(buffer_redo_occ, by = "species") %>% 
  mutate(distance_km_half = distance_km/2)


#now fix final buffer distance file

buffer_good <- buffer_dist %>% filter(distance_km_half < 1650)

buffer_dist_update <- bind_rows(buffer_redo_dist, buffer_good)


## load in large occ dist buffers (calculated on peggy)
buffer_dist_large <- list.files(path = getwd(), pattern = "buffer_dist", full.names = TRUE) %>% 
  lapply(read_csv) %>% bind_rows() %>% dplyr::select(-X1)

buffer_dist_update <- bind_rows(buffer_dist_update, buffer_dist_large)
write.csv(buffer_dist_update, "data/buffer_distances_update.csv")


# REDO thinned large occ --------------------------------------------
load("D:/endophyteDimensions/data/thinned_occ_all_update.RData")

buffer_dist <- read_csv("data/buffer_distances_final.csv")

# filter species that are not in buffer distances file
specs_big <-  thin.occ.update %>% filter(!(species %in% buffer_dist$species)) %>% 
  group_by(species) %>% count() %>% filter(n >= 10) %>% 
  filter(species != 'Cerastium glomeratum') #remove error species

## create occ file with just large occ sp for buffer calc in peggy

thin.occ.large <- thin.occ.update %>% filter(species %in% specs_big$species)
# nearly 1/3 of total occ
save(thin.occ.large, file = "D:/endophyteDimensions/data/thinned_occ_large.RData")

##rename columns of updated file for geodist
names(thin.occ.update)[8:9] <- c("latitude", "longitude")

dist.big <- vector("list", length = length(specs_big))

for (i in 1:length(specs_big)){
  
  s <- dplyr::filter(thin.occ.update, species == specs_big[i])
  
  d <- geodist(as.matrix(s[,c("longitude","latitude")]), sequential = FALSE, measure = "haversine") 
  
  ag <-
    agnes(d,
          diss = TRUE,
          metric = "euclidean",
          method = "single")
  
  dist.big[[i]] <- max(ag$height)/1000
  
  print(i)
  
}

dist.big.all <- dist.big[1:length(specs)] %>% unlist()

## save distances with species names
buffer_dist <- data.frame(species = specs, distance_km = dist.all)
# left join with number of occurrences and cut distances in half to make range polygons

buffer_dist <- buffer_dist %>% as_tibble() %>% left_join(spec_occ, by = "species") %>% select(-X1) %>% 
  mutate(distance_km_half = distance_km/2)

# DID THIS IN PEGGY
buffer_dist_large <- lapply(list.files(pattern = "buffer_dist"), read_csv) %>% bind_rows()

write.csv(buffer_dist_large, "data/buffer_distances_large.csv")
# these are added to 'buffer_distances_final.csv"


##create buffer around native range localities for each species to make geographic extent polygons

range.output = list()
for (i in 1:length(unique(data.occ.sp$Species))){
  s <- subset(data.occ.sp, data.occ.sp$Species == unique(data.occ.sp$Species)[i])
  x <- circles(s@coords, d = buffer_dist$distance_km_half[i], lonlat = T)
  s.range <- gUnionCascaded(x@polygons)
  range.output[[i]]<- s.range
  print(i)
}


# ADD SYNONYM SPECIES ---------------------------------------------------
load("D:/endophyteDimensions/data/syn.thin.final.RData")


##rename columns of updated file for geodist
names(syn.thin)[8:9] <- c("latitude", "longitude")

#filter out species > 14000, just 2 sp
syn.large <- syn.thin %>% group_by(species) %>% count() %>% filter(n > 14000)

specs <- syn.thin %>% filter(!(species %in% syn.large$species)) %>%
  group_by(species) %>% count() %>% filter(n > 1) %>% 
  pull(species) %>% unique()

dist <- vector("list", length = length(specs))

for (i in 1:length(specs)){

  s <- dplyr::filter(syn.thin, species == specs[i])
  
  d <- geodist(as.matrix(s[,c("longitude","latitude")]), sequential = FALSE, measure = "haversine") 
  
  ag <-
    agnes(d,
          diss = TRUE,
          metric = "euclidean",
          method = "single")
  
  if((max(ag$height)/1000) > 3330){
    dist[[i]] <- max(ag$height[ag$height < (3330*1000)]/1000) 
    
  }
  else{
    dist[[i]] <- max(ag$height)/1000
    
  }
  
  print(i)
  
  
}


dist.all <- unlist(dist)

## save distances with species names
syn_buffer_dist <- data.frame(species = specs, distance_km = dist.all)
# left join with number of occurrences and cut distances in half to make range polygons

syn_spec_occ <- syn.thin %>% group_by(species) %>% count() %>% rename(num_occ = n)

syn_buffer_dist_final <- syn_buffer_dist %>% as_tibble() %>% left_join(syn_spec_occ, by = "species") %>% 
  mutate(distance_km_half = distance_km/2)

write.csv(syn_buffer_dist_final, "data/syn_buffer_dist.csv")
