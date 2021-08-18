# removing extra plots in model folders

#cd to model test folder

specs <- read.csv("data/species_list.csv") %>% pull(species)

baseDir <- "/scratch/projects/endophytedim/batch_21_100/results/"

for (i in 21:100){
  
  setwd(baseDir)
  
  setwd(paste0(baseDir, specs[i], "_model/plots"))
  
  files <- list.files(path = getwd(), full.names = TRUE)
  
  files.remove <- grep(pattern = "species_bio", x = files, value = TRUE, invert = TRUE)
  
  file.remove(files.remove)
  
  
  
}


# test first
setwd("model.test/Dubautia ciliolata/bestbeta/results/plots/")

files <- list.files(path = getwd(), full.names = TRUE)

files.remove <- grep(pattern = "species_bio", x = files, value = TRUE, invert = TRUE)

file.remove(files.remove)
