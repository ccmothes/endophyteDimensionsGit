# fixing final dataset


library(tidyverse)


#load in files

#previous final dataset (plant.occ)
load("D:/endophyteDimensions/data/thinned_occ_final.RData")

#added species (syn.thin)
load("D:/endophyteDimensions/data/syn.thin.final.RData")

#first round of removable species
load("data/remove_sp.RData") #455

# load in original species list
original <- read_csv("data/species_original_combined.csv") %>% 
  dplyr::select(family = 'Family Name', species = "Species Name") %>% 
  separate(species, c('genus', 'species'), sep = ' ', extra = 'drop') %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))

new <- read_csv("data/species_new_combined.csv") %>% 
  dplyr::select(family = Family, genus = Genus, species = Species) %>% 
  filter(!(is.na(species))) %>% filter(!(str_detect(species, "sp.$")))


species_list <- rbind(original, new) %>% unite(taxon, c("genus", "species"), sep = " ")

# add new species
new_sp <- data.frame(family = "Amborellaceae", taxon = "Amborella trichopoda")

species_list <- rbind(species_list, new_sp) #THERE ARE 1304 DUPLICATE NAMES HERE



#species key
spec_key <- read_csv("data/species_name_key.csv") %>% dplyr::select(-X1)

#get list of duplicated original species (i.e., originals with multiple synonym match in dataset)
original_dup <- spec_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  arrange(original_sciname) %>% group_by(original_sciname) %>% count() # 51 original names have multiple model names
# are all of these in one of the final datasets?
original_dup %>% filter(original_sciname %in% plant.occ$species) #42 in original file
original_dup %>% filter(original_sciname %in% syn.thin$species) #none here...so 9 of the 51 must have been dropped in cleaning process



#get list of duplicated model name species (i.e., duplicate model species from combining datasets)
model_dup <- spec_key %>% filter(duplicated(model_name) | duplicated(model_name, fromLast = TRUE)) %>% 
  arrange(model_name) %>% group_by(model_name) %>% count() #26....
# do these match multiple original names?
spec_key %>% filter(duplicated(model_name) | duplicated(model_name, fromLast = TRUE)) %>% 
  arrange(model_name) %>% View()
# yes some of them...see where these model names originated from

model_dup %>% filter(model_name %in% plant.occ$species) # 5 in original dataset
model_dup %>% filter(model_name %in% syn.thin$species) #17 in new synonym dataset
model_dup %>% filter(model_name %in% plant.occ$species & model_name %in% syn.thin$species) #no dups



#read in files used to make species key

## this one was not added to species key because names were going to be changed...
name.switch1 <- read_csv("data/name.switch1.csv") %>% dplyr::select(-X1)


# FINAL DATASET ----------------------------------------------------------------

## double check, any dups between syn.thin and plant occ?

syn.thin %>% filter(species %in% plant.occ$species)
plant.occ %>% filter(species %in% syn.thin$species)
#nope!


#but do any syn.thin match to plant.occ sp?
spec_key %>% filter(model_name %in% syn.thin$species) %>% filter(original_sciname %in% plant.occ$species)
# nope!

#how about the opposite
spec_key %>% filter(model_name %in% plant.occ$species) %>% filter(original_sciname %in% syn.thin$species)
# nope! all good


# remove spec. from plant.occ that do not have match to original name


plant.occ.update <- plant.occ %>% filter(!(species %in% remove_sp))


syn.thin %>% filter(!(species %in% spec_key$model_name)) #all in spec.key
plant.occ.update %>% filter(!(species %in% spec_key$model_name)) %>% pull(species) %>% unique()
#need to add "amborella trichopoda" to key. 
spec_key <- spec_key %>% bind_rows(data.frame(model_name = "Amborella trichopoda", original_sciname = "Amborella trichopoda"))

#4 other species in dataset but not in key..these are name.switch1 species

# remove for now, will be put in file for later
plant.occ.later <- plant.occ.update %>% filter(species %in% name.switch1$species)

#update final file
plant.occ.update <- plant.occ.update %>% filter(!(species %in% name.switch1$species))


#now fix key, filter to just model names in plant.occ.update and syn.thin and see if still duplicates
spec_key_update <- spec_key %>% filter(model_name %in% plant.occ.update$species | model_name %in% syn.thin$species)

original_dup_updated <- spec_key_update %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  arrange(original_sciname)
# 47 original names have multiple model names
# save all of this with the plant.occ.later file
spec_key_later <- name.switch1 %>% dplyr::select(model_name = species, original_sciname) %>% bind_rows(original_dup_updated)
# now remove these from species key
spec_key_update <-  spec_key_update %>% filter(!(original_sciname %in% spec_key_later$original_sciname))

plant.occ.later <- plant.occ.update %>% filter(species %in% spec_key_later$model_name | species %in% spec_key_later$original_sciname) %>% 
  bind_rows(plant.occ.later)

plant.occ.update <- plant.occ.update %>% filter(!(species %in% spec_key_later$model_name)) %>% filter(!(species %in% spec_key_later$original_sciname))

# OKAY PLANT.OCC FINAL ALL SET, SAVE THIS AND FILE FOR LATER
write.csv(spec_key_later, "data/species_key_later.csv")
save(plant.occ.later, file = "D:/endophyteDimensions/data/plant.occ.later.RData")


#NOW FIX MODEL DUPLICATES, REMOVE THOSE MATCHING TWO DIFFERENT ORIGINAL NAMES BUT KEEP EXACT MATCHES (assuming that was from original search)
model_dup_updated <- spec_key_update %>% filter(duplicated(model_name) | duplicated(model_name, fromLast = TRUE)) %>% 
  arrange(model_name) 
#22 model names have multiple original names...
# keep cases where duplicates are identical
# find cases where model != any o
model_dup_updated_keep <- model_dup_updated %>% filter(model_name == original_sciname)

model_dup_updated_remove <- model_dup_updated %>% 
  filter(model_name != original_sciname & !(model_name %in% model_dup_updated_keep$model_name)) %>% distinct(model_name)


spec_key_update <- spec_key_update %>% filter(!(model_name %in% model_dup_updated_remove$model_name))

#Use this to remove incorrect matchings (rows)
incorrect_match <- spec_key_update %>% filter(duplicated(model_name) | duplicated(model_name, fromLast = TRUE)) %>% 
  arrange(model_name) %>% filter(model_name != original_sciname)

spec_key_update <- spec_key_update %>% 
  filter(!(model_name %in% incorrect_match$model_name) & !(original_sciname %in% incorrect_match$original_sciname))

syn.thin %>% filter(!(species %in% spec_key_update$model_name)) %>% distinct(species)
#17 syn.thin species now not in key
syn.thin %>% filter(!(species %in% spec_key_update$model_name)) %>% 
  distinct(species) %>% filter(species %in% model_dup_updated_remove$model_name) #13 species

syn.thin %>% filter(!(species %in% spec_key_update$model_name)) %>% 
  distinct(species) %>% filter(species %in% incorrect_match$model_name) #4 species

#remove these from syn.thin
syn.thin.update <- syn.thin %>% filter(species %in% spec_key_update$model_name)


## double check for duplicates
spec_key_update %>% filter(duplicated(model_name) | duplicated(model_name, fromLast = TRUE)) %>% 
  arrange(model_name)  #good here!
spec_key_update %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) #same!


#make sure all model names in final occ dataset are in species key
plant.occ.update %>% filter(!(species %in% spec_key_update$model_name)) #Linum salcutum removed from key for some reason...
syn.thin.update %>% filter(!(species %in% spec_key_update$model_name)) #good here, all in dataset

spec_key_update <- spec_key_update %>% bind_rows(data.frame(model_name = "Linum sulcatum", original_sciname = "Linum sulcatum"))


write.csv(spec_key_update, "data/species_key_update.csv")

plant.occ <- bind_rows(plant.occ.update, syn.thin.update) #woohoo same # rows at key! 7640
save(plant.occ, file = "D:/endophyteDimensions/data/plant.occ.final.RData")



# save final buffer distances file
buffer_dist <- read_csv("data/buffer_distances_update.csv") # this includes all large buffers
buffer_dist_syn <- read_csv("data/syn_buffer_dist.csv")


buffer_dist_all <- bind_rows(buffer_dist, buffer_dist_syn)

#filter to just species in final dataset
buffer_dist_all <- buffer_dist_all %>% filter(species %in% plant.occ$species) #only 7591

# oh no, some species not in here?
plant.occ.update %>% filter(!(species %in% buffer_dist_all$species)) %>% distinct(species) #all included
syn.thin.update %>% filter(!(species %in% buffer_dist_all$species)) %>% distinct(species) # all missing are in here
# oh, forgot to remove species with < 10 occ in syn file...

syn.thin.remove <- syn.thin.update %>% group_by(species) %>% count() %>% filter(n < 10)

plant.occ <- plant.occ %>% filter(!(species %in% syn.thin.remove$species))
plant.occ %>% filter(!(species %in% buffer_dist_all$species)) %>% distinct(species) #two species still missing buffer dist..these were the large syn I forgot to do. Remove for now

buffer_dist_all <- buffer_dist_all %>% filter(species %in% plant.occ$species)

plant.occ <- plant.occ %>% filter(species %in% buffer_dist_all$species)
spec_key_update <- spec_key_update %>% filter(model_name %in% buffer_dist_all$species)

save(plant.occ, file = "D:/endophyteDimensions/data/plant.occ.final.RData")
write.csv(spec_key_update, "data/species_key_update.csv")
write.csv(buffer_dist_all, "data/buffer_dist_model.csv")

# get vector of species
specs <- plant.occ %>% distinct(species) %>% arrange(species) %>% pull(species)
write.csv(specs, "data/species_list.csv")


#FIX SPECIES NAMES FOR FILE SAVING

specs <- read_csv("data/species_list.csv")

plant.occ <- plant.occ %>% mutate(species = str_replace(species, " ", "_"))
specs <- specs %>% mutate(species = str_replace(species, " ", "_"))
buffer_dist <- buffer_dist %>% mutate(species = str_replace(species, " ", "_"))

specs_key <- read_csv("data/species_key_update.csv")
specs_key <- specs_key %>% mutate(model_name = str_replace(model_name, " ", "_"),
                                  original_sciname = str_replace(original_sciname, " ", "_"))


save(plant.occ, file = "D:/endophyteDimensions/data/plant.occ.final.RData")
write.csv(specs, "data/species_list.csv")
write.csv(buffer_dist, "data/buffer_dist_model.csv")
write.csv(specs_key, "data/species_key_update.csv")

