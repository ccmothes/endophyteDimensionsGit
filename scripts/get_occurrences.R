# Retrieve plant occurence data

library(tidyverse)
library(dismo)
library(taxize)


# FILE PREP -------------------------------

# read in files

original <- read_csv("data/species_original_combined.csv") %>% 
  dplyr::select(family = 'Family Name', species = "Species Name") %>% 
  separate(species, c('genus', 'species'), sep = ' ', extra = 'drop') %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))

new <- read_csv("data/species_new_combined.csv") %>% 
  dplyr::select(family = Family, genus = Genus, species = Species) %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))


species_list <- rbind(original, new)


# download occurrences

## make empty character vector to store species names without records
no_occ <- vector(mode = "character", length = nrow(species_list))

## make empty list to store occ in
occ <- vector("list", length = nrow(species_list))


for (i in 482:nrow(species_list)){
  genus <- species_list[i, 2]
  species <- species_list[i, 3]

  download <-
    tryCatch({
      gbif(
        genus = genus,
        species = species,
        geo = TRUE,
        removeZeros = TRUE,
        download = TRUE,
        ntries = 3,
        nrecs = 300,
        end = 190000
      )}, error = function(msg){
        return(NA)
      })
  
  
  
  
  if(is.null(download)){
    no_occ[[i]] <- paste(genus, species, sep = " ")
  } else {
    
    if(is.na(download)){
      no_occ[[i]] <- paste(genus, species, sep = " ")
      
    } else {
      
      if(!("lat" %in% names(download))){
        no_occ[[i]] <- paste(genus, species, sep = " ")
        
      } else {
        
        occ[[i]] <- download %>% 
          dplyr::select(acceptedScientificName, acceptedTaxonKey, basisOfRecord,
                        country, ISO2, lat, lon, year, family, genus, species) %>% 
          filter(!(is.na(lat)) | !(is.na(lon))) %>% distinct(lat, lon, .keep_all = TRUE)
        
        
      }
      
    }
    
    
    
  }
  
 print(i) 
  
}

occ_1 <- do.call("rbind", occ)
write.csv(occ_1, "data/occ_1.csv")

occ_2 <- do.call("rbind", occ)
write.csv(occ_2, "data/occ_2.csv")
#occ_1 = 1:112
#occ_2 = 113:312
#314 had too many species, had to skip
#339 too many
# 438 too many




no_occ <- no_occ[no_occ != ""]



download <-
  tryCatch({
    gbif(
      genus = genus,
      species = species,
      geo = TRUE,
      removeZeros = TRUE,
      download = TRUE,
      ntries = 3,
      nrecs = 300,
      end = 190000
    )}, error = function(msg){
      return(NA)
    })

