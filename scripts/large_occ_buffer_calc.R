# Large OCC buffer calculation

# Set Up -----------------------------------------------------------
library(raster)
library(geodist)
library(cluster)
library(readr)
library(dplyr)

setwd("/scratch/projects/endophytedim/")

# read in occ
load("data/thinned_occ_large.RData")

specs <- thin.occ.large %>% pull(species) %>% unique()

##rename columns of updated file for geodist
names(thin.occ.large)[8:9] <- c("latitude", "longitude")

dist <- vector("list", length = length(specs))

for (i in 1:4){
  
  s <- dplyr::filter(thin.occ.large, species == paste(specs[i]))
  
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

dist <- unlist(dist)

## save distances with species names
buffer_dist <- data.frame(species = specs[1:4], distance_km = dist)

# left join with number of occurrences and cut distances in half to make range polygons
buffer_dist_final <- thin.occ.large %>% group_by(species) %>% count() %>% 
  filter(species %in% buffer_dist$species) %>% left_join(buffer_dist, by = "species") %>% 
  mutate(distance_km_half = distance_km/2) %>% rename(num_occ = n)

write.csv(buffer_dist_final, "data/buffer_dist_test.csv")

print("done")

