# download occurences from rgbif (gbif recommends)

# pulled from this example: https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

library(tidyverse)
library(rgbif)
library(taxize)
library(vroom)

#fill in gbif user credentials

user <- "ccmothes"
pwd <- "maxNsunny1!"
email <- "ccmothes@gmail.com"


# import data -----------------------------------------------------------


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

#match the taxon names -------------------------------------------------------

gbif_taxon_keys <- species_list %>% 
  pull("taxon") %>% 
  taxize::get_gbifid_(method = "backbone") %>% 
  imap(~ .x %>% mutate(original_sciname = .y)) %>% 
  bind_rows() %T>%
  readr::write_tsv(path = "data/all_matches.tsv") %>% 
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  filter(kingdom == "Plantae") %>% 
  pull(usagekey)

#get taxon key for added sp A.trich

a.trich.key <- new_sp %>% pull(taxon) %>%
  taxize::get_gbifid_(method = "backbone") %>% 
  imap(~ .x %>% mutate(original_sciname = .y)) %>% 
  bind_rows() %T>%
  readr::write_tsv(path = "data/a.trichopoda_matches.tsv") %>% 
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  filter(kingdom == "Plantae") %>% 
  pull(usagekey)
  
# only 10330 (10330 with A.trich) species with taxon keys. With other method, at 3500 species 1746
# did not have occurrences, so so far this is looking better

#download occurrences ------------------------------------------------------------

occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

# get occ for A.trichopoda
occ_download(
  pred_in("taxonKey", a.trich.key),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)


# file is now in gbif account 8gb with 77 MILLION occurrences

# read in new file


rgbif_occ <- read_csv("D:/endophyteDimensions/rgbif_occurrences/rgbif_occ.csv")
# quit after an hour and only 70% loaded

rgbif_occ <- vroom("D:/endophyteDimensions/rgbif_occurrences/rgbif_occ.csv")
# still takes forever.....


# UPDATE download occ -----------------------------------------------------------


## split up search into 12 groups to read and clean separately then combine summary stats

## create list of split taxon keys to run multiple occurrences searches

taxon_keys_split <- split(gbif_taxon_keys, ceiling(seq_along(gbif_taxon_keys)/1033))


## download occ for each chunk of taxon keys

## need to do in sests of 3 (only 3 download allowed at a time)

for (i in 10:10){
  
  occ_download(
    pred_in("taxonKey", taxon_keys_split[[i]]),
    pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email
  )
  
  print(i)
  
}

#last file (occ_10) is too large to clean, split that up into three sets

gbif_taxa <- read_tsv("data/all_matches.tsv") %>% 
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  filter(kingdom == "Plantae") %>% 
  pull(usagekey)


taxon_keys_split <- split(gbif_taxa, ceiling(seq_along(gbif_taxa)/1033))

#split last section into 3's
split_10 <- split(taxon_keys_split[[10]], ceiling(seq_along(taxon_keys_split[[10]])/345))

# re-download occ_10 in 3 sections

for (i in 1:length(split_10)){
  
  occ_download(
    pred_in("taxonKey", split_10[[i]]),
    pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email
  )
  
  print(i)
  
}



# Missing taxa -----------------------------------------------------------

## finall dataset used to search occ
gbif_taxa <- read_tsv("data/all_matches.tsv") %>% filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  filter(kingdom == "Plantae")

## full dataset returned from taxize
all_matches <- read_tsv("data/all_matches.tsv")


missing <- species_list %>% filter(!(taxon %in% gbif_taxa$original_sciname)) %>% rename(species = taxon)
# 1745 -> 1739 (previous version searched those not in gbif_taxa$species)
write.csv(missing, "data/missing_gbif_species_update.csv")


#19 species not found anywhere (NOTE THERE ARE OVER 1,000 DUPLICATE NAMES IN SP FILE)
# TOTAL # SPECIES = 10765
species_list %>% filter(!(taxon %in% all_matches$original_sciname))


# get comparison file (original vs accepted/gbif name)
name_compare <- gbif_taxa %>% dplyr::select(original_sciname, species) #all have match except one (capitalization)


## taxa problem #1 ----------------------------------------------
# uh-oh, there are some species names in occ file that are not in the original species list (ignore a. trichopoda, in separate file)
# gbif search return synonym names? Investigate...

#using final occ dataset
load("D:/endophyteDimensions/data/thinned_occ_final.RData")

#load taxonomy file
all_matches <- read_tsv("data/all_matches.tsv")

x <- plant.occ %>% filter(!(species %in% species_list$taxon)) %>% pull(species) %>% unique()
# 510 missing names

#see how many occ for each species
plant.occ %>% filter(!(species %in% name_compare$species)) %>% group_by(species) %>% count() %>% View()
  

all_matches %>% filter(matchtype == "EXACT" & status == "SYNONYM") %>% 
  filter(species %in% x) %>% distinct(species) #14 were synonym match of original species names

#are these original species names in species list?
all_matches %>% filter(matchtype == "EXACT" & status == "SYNONYM") %>% 
  filter(species %in% x) %>% filter(original_sciname %in% species_list$taxon) %>% 
  distinct(species)
#yes all of them....

syn_dup <- all_matches %>% filter(matchtype == "EXACT" & status == "SYNONYM") %>% 
  filter(species %in% x) %>% distinct(species, original_sciname, .keep_all = TRUE)

#BUT only 6 have the original name in occ dataset (panicum coloratum has 2 syn spec included)
plant.occ %>% filter(species %in% syn_dup$original_sciname) %>% distinct(species)

# so refine syn. dup 
syn.occ <- plant.occ %>% filter(species %in% syn_dup$original_sciname) %>% pull(species) %>% unique()
syn_dup2 <- syn_dup %>% filter(original_sciname %in% syn.occ)

#are the localities the same??
occ.syn <- plant.occ %>% filter(species %in% syn_dup2$species)

occ.original <- plant.occ %>% filter(species %in% syn_dup2$original_sciname)

occ.syn %>% group_by(species) %>% 
  filter(!(decimallatitude %in% occ.original$decimallatitude) & !(decimallongitude %in% occ.original$decimallongitude)) %>% 
  count() # 5 syn species have (a lot) of unique occ
# species                    n
# <chr>                  <int>
#   1 Panicum antidotale       213
# 2 Panicum atrosanguineum    23
# 3 Panicum repens          1245
# 4 Panicum virgatum        2562
# 5 Rhododendron ponticum   7066

# okay for these 5 species get matching original name to change in occ dataset
## WILL NEED TO REDO BUFFER CALC FOR THESE
name.switch1 <- occ.syn %>% group_by(species) %>% 
  filter(!(decimallatitude %in% occ.original$decimallatitude) & !(decimallongitude %in% occ.original$decimallongitude)) %>% 
  pull(species) %>% unique()

name.switch1 <- filter(syn_dup, species %in% name.switch1)
#WHAT TO DO WITH P. ANTIDOTALE matches with 2 original names. Remove? YEs
name.switch1 <- name.switch1 %>% filter(species != "Panicum antidotale")
write.csv(name.switch1, "data/name.switch1.csv")

#Now get second list of dups that had matching original name but not found in occ set
# i.e., just switch name NO NEED TO RECALCULATE ANYTHING
name.switch2 <- syn_dup %>% filter(!(species %in% name.switch1$species))
write.csv(name.switch2, "data/name.switch2.csv")

#okay now working with remaining unknown species
occ_unknown_sp <- x[!(x %in% syn_dup$species)]


# retrieve gbif synonyms for each species
unknown_gbif_syn <- occ_unknown_sp %>% 
  taxize::get_gbifid_(method = "backbone") %>% 
  imap(~ .x %>% mutate(original_sciname = .y)) %>% 
  bind_rows() #%T>%

write_tsv(unknown_gbif_syn, "data/all_unknown_matches.tsv")
  #readr::write_tsv(path = "data/all_unknown_matches.tsv") 
  # filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  # filter(kingdom == "Plantae") %>% 
  # pull(usagekey)

unknown_gbif_syn_exact <- unknown_gbif_syn %>% filter(matchtype == "EXACT" & status != "DOUBTFUL")

#see how many have exact synonym match
unknown_gbif_syn_exact %>% group_by(original_sciname) %>% 
  filter(species %in% species_list$taxon) 
# only 13...some match with multiple (19 total rows)

#see how many have any sort of match with original names (might be doubtful/fuzzy match)
unknown_gbif_syn %>% group_by(original_sciname) %>% 
  filter(species %in% species_list$taxon) 
# 91...with a lot having multple matches (148 rows...)

#just keep list of species names that have single orignal match
unknown_with_match <- unknown_gbif_syn %>% group_by(original_sciname) %>% 
  filter(species %in% species_list$taxon) %>% 
  count() %>% 
  filter(n == 1) #52 total species

name.switch3 <- unknown_gbif_syn %>% filter(original_sciname %in% unknown_with_match$original_sciname) %>% 
  filter(species %in% species_list$taxon)
write.csv(name.switch3, "data/name.switch3.csv")
  
#GET LIST OF SPECIES TO REMOVE

## these are species not listed in any of the name.switches
remove_sp <- x[!(x %in% name.switch1$species)] %>% .[!(. %in% name.switch2$species)] %>% 
  .[!(. %in% name.switch3$original_sciname)]
# 445 remove
save(remove_sp, file = "data/remove_sp.RData")

unknown_sp_name_match <- bind_rows(name.switch2, name.switch3)
write.csv(unknown_sp_name_match, "data/unknown_sp_name_match.csv")

## taxa problem #2 -------------------------------------------------------

#...there were a bunch of species kicked out because they were 'synonym' not 'accepted'...
all_matches %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(matchtype == "EXACT" & status == "SYNONYM") %>% distinct(original_sciname) #1265 species...


# see how many synonyms per species...
all_matches %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(matchtype == "EXACT") %>% group_by(original_sciname) %>% count() %>% View()

all_matches %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(matchtype == "EXACT") %>% group_by(original_sciname) %>% count() %>%
  filter(n > 1) #352 more than 1 synonym, how to choose?

#inspect one with 18 synonyms
all_matches %>% filter(original_sciname == "Kyllinga brevifolia") %>% filter(matchtype == "EXACT")

# see how many of these had a synonym included in the final dataframe... 269
all_matches %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(species %in% plant.occ$species) %>% distinct(original_sciname)

# filter missing species out that have same synonym as another in the occ dataset (n = 263)
missing_syn_included <- all_matches %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(species %in% plant.occ$species) %>% filter(species %in% species_list$taxon)


#okay remove those from the dataset and see what remaining 'EXACT" species working with..1067
all_matches %>% 
  filter(matchtype == "EXACT") %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(!(original_sciname %in% missing_syn_included$original_sciname)) %>% distinct(original_sciname)

#how many have a single synonym we could search (not really a way to choose among multiple and don't want to search the wrong species)
missing_single_syn <- all_matches %>% 
  filter(matchtype == "EXACT") %>% 
  filter(!(original_sciname %in% gbif_taxa$original_sciname)) %>% 
  filter(!(original_sciname %in% missing_syn_included$original_sciname)) %>% 
  group_by(original_sciname) %>% count() %>% filter(n <= 1)
#788....worth re-doing search??

syn_search_sp <- all_matches %>% filter(original_sciname %in% missing_single_syn$original_sciname) %>% 
  filter(matchtype == "EXACT" & status == "SYNONYM")
# save to match synonym names to later
write.csv(syn_search_sp, "data/synonym_species.csv")

# NEW OCC DOWNLOAD FOR THESE SPECIES
taxon_keys_new <- syn_search_sp %>% pull(usagekey)

occ_download(
  pred_in("taxonKey", taxon_keys_new),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)


# make species name key file ----------------------


## first/original occ search, but filter only those that returned occ so far
gbif_taxa <- read_tsv("data/all_matches.tsv") %>% filter(matchtype == "EXACT" & status == "ACCEPTED") %>% #get only accepted and matched names
  filter(kingdom == "Plantae") %>% select(species, original_sciname) %>% distinct() %>% 
  rename(model_name = species) %>% 
  filter(model_name %in% plant.occ$species)

## add name_switch1? because there are also species in other datasets that are dups of already included original species



##name switch 2 = species was the synonym name found in dataset, original
## was the original sp name
name_switch2 <- read.csv("data/name.switch2.csv") %>% select(species, original_sciname) %>% 
  rename(model_name = species) # problem, two models per one original name match...special case

#example of how to filter out the duplicated rows
name_switch2 %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE))


name_switch3 <- read.csv("data/name.switch3.csv") %>% select(species, original_sciname) 

name_switch3 %>% filter(species %in% species_list$taxon) 
# all species names are original names, original_sci name is the dataset synonym that was
# previously unmatched
name_switch3 <- name_switch3 %>% rename(model_name = original_sciname, original_sciname = species)


# read in the additional synonym searched species (those dropped the first round because no exact name match)
syn_search_sp <- read.csv("data/synonym_species.csv") %>% select(species, original_sciname) %>% 
  rename(model_name = species) %>% as_tibble() #no duplicates, all names in original list


#now bind all to make species name key
species_name_key <- bind_rows(gbif_taxa, name_switch2, name_switch3, syn_search_sp)
write.csv(species_name_key, "data/species_name_key.csv")

# search for duplicates
species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  arrange(original_sciname) %>% group_by(original_sciname) %>% count() # 47 original names have multiple model names


# see where these model names are coming from
species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  filter(model_name %in% syn_search_sp$model_name)
# 4 here, but these were the ones getting removed from syn dataset

species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  filter(original_sciname %in% name_switch2$original_sciname)
# 6 here...

species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  filter(original_sciname %in% gbif_taxa$original_sciname)
# 42 here...
# are any of these in another dataset? is that where the duplicates come from?
species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  filter(original_sciname %in% gbif_taxa$original_sciname) %>% filter(original_sciname %in% name_switch2$original_sciname)
# all 42 of them in switch3

species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  View()

#name_switch2 all in switch3?
species_name_key %>% filter(duplicated(original_sciname) | duplicated(original_sciname, fromLast = TRUE)) %>% 
  filter(original_sciname %in% name_switch2$model_name) %>% filter()



