# clean occurrences

library(tidyverse)
library(vroom)
library(CoordinateCleaner)
library(countrycode)
#library(spThin)
library(raster)
#remotes::install_github('adamlilith/enmSdm', dependencies=TRUE)
library(enmSdm)

# read in occurrence dataframes -----------------------------------------


occ_1 <- vroom("D:/endophyteDimensions/data/rgbif_occ_1.csv", 
               col_select = c(datasetKey, scientificName, taxonKey, speciesKey,
                              family, genus, species,
                              decimalLatitude, decimalLongitude, basisOfRecord,
                              countryCode, eventDate, year))

## remove duplicates, zeros and NAs
occ_1_filtered <- occ_1 %>% distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>% 
  filter(decimalLatitude != 0) %>% filter(decimalLongitude != 0) %>% 
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))


## coordinate cleaner, remove all flagged records ------------------------

### to use country - coordinate mismatch flag, convert country code
### convert country code from ISO2c to ISO3c
occ_1_filtered$countryCode <-  countrycode(occ_1_filtered$countryCode, origin =  'iso2c', destination = 'iso3c')


#to avoid specifying it in each function
names(occ_1_filtered)[8:9] <- c("decimallatitude", "decimallongitude")

clean <- occ_1_filtered%>%
  cc_val()%>%
  cc_equ()%>%
  cc_cap()%>%
  cc_cen()%>%
  cc_coun(iso3 = "countryCode")%>%
  cc_gbif()%>%
  cc_inst()%>%
  cc_sea(ref = buffland)%>%
  cc_zero()%>%
  cc_outl()%>%
  cc_dupl()

## spatial thin with 'thin' to 5km, 1 rep

occ_clean_thin <- thin(clean, lat.col = "decimallatitude",
                       long.col = "decimallongitude",
                       spec.col = "species",
                       thin.par = 5,
                       reps = 0,
                       locs.thinned.list.return = TRUE,
                       write.files = TRUE,
                       out.dir = "")

## ERROR file size too large, try loop over species

specs <- clean %>% pull(species)
occ_clean_thin <- vector("list", length = length(specs))

for (i in 2:length(specs)){
  
  dat <- filter(clean, species == specs[i])
  
  occ_clean_thin[[i]] <- thin(dat, lat.col = "decimallatitude",
                         long.col = "decimallongitude",
                         spec.col = "species",
                         thin.par = 5,
                         reps = 0,
                         locs.thinned.list.return = TRUE)
  
  print(i)
  
}

# still will not work...


## summarize # occ per species

species_occ_1 <- clean %>% group_by(species) %>% count()

## save files
write_csv(species_occ_1, "data/species_occ_1.csv")
write_csv(clean, "D:/endophyteDimensions/clean_occ_1.csv")


## clean all occurrences -----------------------------------------------------


# make list of file names

files <- list.files(path = "D:/endophyteDimensions/data/", pattern = "rgbif_occ", full.names = TRUE)


# remove occ_1 cause already did that and 10 cause too large (R crashes)
files <- files[-c(1,2)]

#re-do occ_10, just keep 10 1:3
files <- files[c(3,4,5)]


for (i in 1:length(files)){
  
  
  # read in occurrence dataframes
  
  occ <- vroom(files[i], 
                 col_select = c(datasetKey, scientificName, taxonKey, speciesKey,
                                family, genus, species,
                                decimalLatitude, decimalLongitude, basisOfRecord,
                                countryCode, eventDate, year))
  
  ## remove duplicates, zeros and NAs
  occ_filtered <- occ %>% distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>% 
    filter(decimalLatitude != 0) %>% filter(decimalLongitude != 0) %>% 
    filter(!is.na(decimalLongitude)) %>%
    filter(!is.na(decimalLatitude))
  
  ## run through coordinate cleaner, remove all flagged records
  
  ### to use country - coordinate mismatch flag, convert country code
  ### convert country code from ISO2c to ISO3c
  occ_filtered$countryCode <-  countrycode(occ_filtered$countryCode, origin =  'iso2c', destination = 'iso3c')
  
  
  #to avoid specifying it in each function
  names(occ_filtered)[8:9] <- c("decimallatitude", "decimallongitude")
  
  clean <- occ_filtered%>%
    cc_val()%>%
    cc_equ()%>%
    cc_cap()%>%
    cc_cen()%>%
    cc_coun(iso3 = "countryCode")%>%
    cc_gbif()%>%
    cc_inst()%>%
    cc_sea(ref = buffland)%>%
    cc_zero()%>%
    cc_outl()%>%
    cc_dupl()
  

  ## summarize # occ per species
  
  species_occ <- clean %>% group_by(species) %>% count()
  
  ## save files
  write_csv(species_occ, paste0("data/species", str_sub(files[i], start = 34L)))
  write_csv(clean, paste0("D:/endophyteDimensions/data/clean", str_sub(files[i], start = 34L)))
  
  
  print(i)
  
}

## ADD A.TRICHOPODA ----------------------------

occ <- vroom("D:/endophyteDimensions/data/rgbif_a.trichopoda.csv", 
             col_select = c(datasetKey, scientificName, taxonKey, speciesKey,
                            family, genus, species,
                            decimalLatitude, decimalLongitude, basisOfRecord,
                            countryCode, eventDate, year))

## remove duplicates, zeros and NAs
occ_filtered <- occ %>% distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>% 
  filter(decimalLatitude != 0) %>% filter(decimalLongitude != 0) %>% 
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

## run through coordinate cleaner, remove all flagged records

### to use country - coordinate mismatch flag, convert country code
### convert country code from ISO2c to ISO3c
occ_filtered$countryCode <-  countrycode(occ_filtered$countryCode, origin =  'iso2c', destination = 'iso3c')


#to avoid specifying it in each function
names(occ_filtered)[8:9] <- c("decimallatitude", "decimallongitude")

clean <- occ_filtered%>%
  cc_val()%>%
  cc_equ()%>%
  cc_cap()%>%
  cc_cen()%>%
  cc_coun(iso3 = "countryCode")%>%
  cc_gbif()%>%
  cc_inst()%>%
  cc_sea(ref = buffland)%>%
  cc_zero()%>%
  cc_outl()%>%
  cc_dupl()

#save file as occ_11

write.csv(clean, "D:/endophyteDimensions/data/clean_occ_11.csv")



# combine species files ---------------------------------------------


species_occ <- list.files(path = "data/", pattern = "species_occ", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows() %>% 
  distinct(species, .keep_all = TRUE) #for some reason there were duplicates...

# add a.trich
a.trich.occ <- a.trich %>% group_by(species) %>% count() %>% rename(num_occ = n)

species_occ <- bind_rows(a.trich.occ, species_occ)

write.csv(species_occ, "data/species_occ_total.csv")

#See which ones are duplicated

species_occ %>% group_by(species) %>% count %>% filter(n > 1)

# check to see if these are duplicate # occ or should be summed?

species_occ %>% filter(species == "Abelmoschus manihot")
species_occ %>% filter(species == "Acer amplum")
species_occ %>% filter(species == "Aconitum kusnezoffii")
# seem similar...just remove duplicates. BUT in final df will need to remove duplicate
# points again


# analyze # occurrences

## number less than 10
species_occ %>% filter(num_occ < 10) # 2,508

## number greater than 50,000
species_occ %>% filter(num_occ > 50000) # 179

## number greater than 20,000
species_occ %>% filter(num_occ > 20000) # 299

species_occ %>% filter(num_occ > 100000) # 92 > 100,000




# combine cleaned occurrences ------------------------------------------------

cleaned_occ <- list.files(path = "D:/endophyteDimensions/data/", pattern = "clean_occ", full.names = TRUE) %>% 
  lapply(vroom) %>% 
  bind_rows() 

# remove duplicate species coordinates

cleaned_occ <- cleaned_occ %>% distinct(species, decimallatitude, decimallongitude, .keep_all = TRUE)

write.csv(cleaned_occ, "D:/endophyteDimensions/data/cleaned_occ_all.csv")
# ERROR file too big, save as R object?
save(cleaned_occ, file = "D:/endophyteDimensions/data/cleaned_occ_all.RData")



# missing species -------------------------------------------------

## get list of original species

original <- read_csv("data/species_original_combined.csv") %>% 
  dplyr::select(family = 'Family Name', species = "Species Name") %>% 
  separate(species, c('genus', 'species'), sep = ' ', extra = 'drop') %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))

new <- read_csv("data/species_new_combined.csv") %>% 
  dplyr::select(family = Family, genus = Genus, species = Species) %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))


species_list <- rbind(original, new) %>% unite(taxon, c("genus", "species"), sep = " ")

# get missing gbif taxonomy list

missing_taxon <- read_csv("data/missing_gbif_species.csv")


no_occur <-
  species_list %>% filter(!(taxon %in% cleaned_occ$species)) %>%
  filter(!(taxon %in% missing_taxon$species))

write.csv(no_occur, "data/missing_occurrences.csv")


# thin from clim raster -------------------------------------------------------

load("D:/endophyteDimensions/data/cleaned_occ_all.RData")

#rbind with a.trich
a.trich <- read_csv("D:/endophyteDimensions/data/clean_occ_11.csv") %>% dplyr::select(-X1)

cleaned_occ <- bind_rows(cleaned_occ, a.trich)

# read in raster after first download
#clim <- getData(name = "worldclim", var = 'bio', res = 2.5, download = TRUE, path = "D:/endophyteDimensions/data")

specs <- cleaned_occ %>% pull(species) %>% unique()


thin.occ <- vector("list", length = length(specs))

for (i in 9376:length(specs)){
  
  occ <- filter(cleaned_occ, species == specs[i])
  
  thin.occ[[i]] <- elimCellDups(as.data.frame(occ), clim, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

# something wrong with specs[954], skipped for now
# specs[1005]
# 2417 too big?
#4930 too big
# 4969
#5853
# 7564
#8466
# 9353
#9375

## skipping these above...bind all the onces that did work

thin.occ.all <- bind_rows(thin.occ)

save(thin.occ.all, file = "D:/endophyteDimensions/data/thinned_occ_all.RData")


# investigate missing ones

missing <- c(954, 1005, 2417, 4930, 4969, 5853, 7564, 8466, 9353, 9375)

specs[missing]

# in order of missing
filter(cleaned_occ, species %in% specs[missing]) %>% group_by(species) %>% count()


#try re-doing  thin with missing species

specs_missing <- specs[missing]

thin.occ.missing <- vector("list", length = length(specs_missing))

for (i in 3:length(specs_missing)){
  
  occ <- filter(cleaned_occ, species == specs_missing[i])
  
  thin.occ.missing[[i]] <- elimCellDups(as.data.frame(occ), clim, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

# 3:10 worked, must have been a memory issue during loop.
# 1:2 might just not have any cell duplicates? Not sure what's going on there

## bind everyting to thin.occ

thin.occ.missing[[1]] <- filter(cleaned_occ, species == specs_missing[1])
thin.occ.missing[[2]] <- filter(cleaned_occ, species == specs_missing[2])


thin.occ.final <- bind_rows(thin.occ.all, thin.occ.missing)


thin.occ.final %>% pull(species) %>% unique() %>% length() # only 10124
# missing one species...see which one

#try error species using spThin Cerastium glomeratum'

c.glom.thin <- spThin::thin(as.data.frame(c.glom), lat.col = "decimallatitude",
                            long.col = "decimallongitude", spec.col = "species",
                            thin.par = 5,
                            reps=1, write.files = FALSE)
## still doesn't work....have to just remove



# there was NAs in the original cleaned occ file. Thinning has removed those


# make new spec occ files

species_occ_update <- thin.occ.final %>% group_by(species) %>% count() %>% rename(num_occ = n)

## remove the weird species combos
filter(species_occ_update, !(species %in% species_occ$species))

thin.occ.final.fix <- as_tibble(thin.occ.final) %>% 
  dplyr::filter(species != "Glyceria declinata × fluitans") %>% 
  dplyr::filter(species != "Glyceria fluitans × notata")

## save thinned and cleaned final dataset (updating old thinned_occ_all file)
save(thin.occ.final.fix, file = "D:/endophyteDimensions/data/thinned_occ_all.RData")


## new species occ summary (after removing two species above)

species_occ_update <- thin.occ.final.fix %>% group_by(species) %>% count() %>% rename(num_occ = n)

write.csv(species_occ_update, "data/species_occ_totalUpdate.csv")

# analyze # occurrences

## number less than 10
species_occ_update %>% filter(num_occ < 10) # 2,508 -> 2,614. Thinning removed 106 more sp.

## number greater than 50,000
species_occ_update %>% filter(num_occ > 50000) # 179 -> 80

## number greater than 20,000
species_occ_update %>% filter(num_occ > 20000) # 299 -> 222

## number greater than 14 k (geo extent cut-off)
species_occ_update %>% filter(num_occ > 14000) # 299 -> 222


species_occ_update %>% filter(num_occ > 100000) # 275


## thin big occ ------------------------------------------------

clim <- raster("D:/endophyteDimensions/data/wc2-5/bio1.bil")

# aggregate

clim_ag <- aggregate(clim, fact = 2, fun = max)


## get species with > 14000 occ

specs_big <- spec_occ %>% filter(num_occ >= 14000) %>% pull(species) %>% unique()


# write loop for thinning with bigger pixel size


thin.occ.big <- vector("list", length = length(specs_big))

for (i in 58:length(specs_big)){
  
  occ <- filter(thin.occ.final.fix, species == specs_big[i])
  
  thin.occ.big[[i]] <- elimCellDups(as.data.frame(occ), clim_ag, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

#57 didn't work
#Cerastium glomeratum
filter(spec_occ, species == specs_big[57]) # has 78,000 occ

thin.occ.big <- bind_rows(thin.occ.big)

# check new occurrence #s

specs_occ_big <- thin.occ.big %>% group_by(species) %>% count()

## how many under 14k now
specs_occ_big %>% filter(n < 14000) #only 59


## take remaining species and thin by ~20km res (factor of 4)

clim_20 <- aggregate(clim_ag, fact = 2, fun = max)

specs_remain <- specs_occ_big %>% filter(n >= 14000) %>% pull(species)


thin.occ.big2 <- vector("list", length = length(specs_remain))

for (i in 1:length(specs_remain)){
  
  occ <- filter(thin.occ.final.fix, species == specs_remain[i])
  
  thin.occ.big2[[i]] <- elimCellDups(as.data.frame(occ), clim_20, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

thin.occ.big2 <- bind_rows(thin.occ.big2)

specs_occ_big2 <- thin.occ.big2 %>% group_by(species) %>% count()

## how many under 14k now
specs_occ_big2 %>% filter(n < 14000) # 147


## filter remaining points to 40k

clim_40 <- aggregate(clim_20, fact = 2, fun = max)


specs_remain2 <-  specs_occ_big2 %>% filter(n >= 14000) %>% pull(species)


thin.occ.big3 <- vector("list", length = length(specs_remain2))

for (i in 1:length(specs_remain2)){
  
  occ <- filter(thin.occ.final.fix, species == specs_remain2[i])
  
  thin.occ.big3[[i]] <- elimCellDups(as.data.frame(occ), clim_40, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

thin.occ.big3 <- bind_rows(thin.occ.big3)

specs_occ_big3 <- thin.occ.big3 %>% group_by(species) %>% count()

## how many under 14k now
specs_occ_big3 %>% filter(n < 14000) # 67

# any over 14000?
specs_occ_big3 %>% filter(n >= 14000) # one species is 14026...see if distance function will run


library(geodist)

s <- dplyr::filter(thin.occ.big3, species == "Achillea millefolium")
names(s)[8:9] <- c("latitude", "longitude")

d <- geodist(as.matrix(s[,c("longitude","latitude")]), sequential = FALSE, measure = "haversine") 
# Yay it works!


#now save final file, and summary file of the filter used for each species

round1 <- filter(specs_occ_big, n < 14000) # 59
round2 <- filter(specs_occ_big2, n < 14000) #147
round3 <- filter(specs_occ_big3, n < 14100) #68, the rest of them

thin_summary <- bind_rows(round1, round2, round3) %>% 
  rename(num_occ = n) %>% 
  mutate(thin_distance = case_when(species %in% round1$species ~ 10, 
                                   species %in% round2$species ~ 20,
                                   species %in% round3$species ~ 40))

write.csv(thin_summary, "data/species_large_thin_summary.csv")


round1.occ <- thin.occ.big %>% filter(species %in% round1$species)
round2.occ <- thin.occ.big2 %>% filter(species %in% round2$species)
round3.occ <- thin.occ.big3 %>% filter(species %in% round3$species) #formality

thin.occ.big.all <- bind_rows(round1.occ, round2.occ, round3.occ)

# now remove previous occurrence sets of these species in original file
data1 <- thin.occ.final.fix %>% filter(!(species %in% thin_summary$species))

# and bind with new data

thin.occ.update <- bind_rows(data1, thin.occ.big.all)

# save as most recent occurrence and species file
## but note that one species still has 78k occ....

save(thin.occ.update, file = "D:/endophyteDimensions/data/thinned_occ_all_update.RData")

specs_occ_total_update <- thin.occ.update %>% group_by(species) %>% count() %>% 
  rename(num_occ = n)

write.csv(specs_occ_total_update, "data/species_occ_totalUpdate2.csv")


# CLEAN SYNONYM SPECS --------------------------------------------------

syn.occ <- vroom("D:/endophyteDimensions/data/syn_occ.csv", 
             col_select = c(datasetKey, scientificName, taxonKey, speciesKey,
                            family, genus, species,
                            decimalLatitude, decimalLongitude, basisOfRecord,
                            countryCode, eventDate, year))

## remove duplicates, zeros and NAs
syn.occ_filtered <- syn.occ %>% distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>% 
  filter(decimalLatitude != 0) %>% filter(decimalLongitude != 0) %>% 
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

## run through coordinate cleaner, remove all flagged records

### to use country - coordinate mismatch flag, convert country code
### convert country code from ISO2c to ISO3c
syn.occ_filtered$countryCode <-  countrycode(syn.occ_filtered$countryCode, origin =  'iso2c', destination = 'iso3c')


#to avoid specifying it in each function
names(syn.occ_filtered)[8:9] <- c("decimallatitude", "decimallongitude")

syn.clean <- syn.occ_filtered%>%
  cc_val()%>%
  cc_equ()%>%
  cc_cap()%>%
  cc_cen()%>%
  cc_coun(iso3 = "countryCode")%>%
  cc_gbif()%>%
  cc_inst()%>%
  cc_sea(ref = buffland)%>%
  cc_zero()%>%
  cc_outl()%>%
  cc_dupl()


## summarize # occ per species

syn_species_occ <- syn.clean %>% group_by(species) %>% count() #702 species left

## save files
write_csv(species_occ, paste0("data/species", str_sub(files[i], start = 34L)))
write_csv(clean, paste0("D:/endophyteDimensions/data/clean", str_sub(files[i], start = 34L)))



## thin occurrences

clim <- raster("D:/endophyteDimensions/data/wc2-5/bio1.bil")


specs <- syn.clean %>% pull(species) %>% unique()


syn.thin.occ <- vector("list", length = length(specs))

for (i in 1:length(specs)){
  
  occ <- filter(syn.clean, species == specs[i])
  
  syn.thin.occ[[i]] <- elimCellDups(as.data.frame(occ), clim, longLat = c("decimallongitude", "decimallatitude"))
  
  
  print(i)
  
}

syn.thin <- bind_rows(syn.thin.occ)

#see if any of these species are found in original data (plant.occ)
syn.thin %>% group_by(species) %>% count() %>% filter(n >= 10) #462 species qualify for model run

syn.thin %>% filter(species %in% plant.occ$species) %>% pull(species) %>% unique()
# 4 rhododendron species here were also in original/final occ file
# are these found in the original sp list?
syn.thin %>% filter(species %in% plant.occ$species) %>% 
  filter(species %in% species_list$taxon) %>% pull(species) %>% unique()
# no...so there are 458 NEW species

# see if there are any new occ for these 4 dup species
dups <- syn.thin %>% filter(species %in% plant.occ$species) %>% pull(species) %>% unique()

syn.thin.dups <- as_tibble(syn.thin) %>% filter(species %in% dups)
plant.occ.dups <- plant.occ %>% filter(species %in% dups)

# any syn.thin.dups occ unique?

syn.thin.dups %>% filter(!(decimallatitude %in% plant.occ.dups$decimallatitude) & !(decimallongitude %in% plant.occ.dups$decimallongitude))
# 14 unique points for 2 species (7 and 7)

# lets just remove all dup species from syn.occ dataset, already handled in final dataset
syn.thin <- syn.thin %>% filter(!(species %in% dups))

syn.thin %>% group_by(species) %>% count()

# save final cleaned/thinned dataset

save(syn.thin, file = "D:/endophyteDimensions/data/syn.thin.final.RData")
