#google drive connection test

library(googledrive)

drive_find(n_max = 10)


#now authorized, time to test

test <- data.frame(site = 1:5, name = c("a", "b", "c", "d", "e"))

write.csv(test, "test.csv")
write.csv(test, "species_1_model/test2.csv")


# set path to google drive folder
drive <- "https://drive.google.com/drive/u/0/folders/0AGb7RkShCg1VUk9PVA/"

#set path to drive model results folder
drive_results <- "https://drive.google.com/drive/u/0/folders/1osfjlJP5dCoBp6COXMveidomPZPD_xGi"

drive_upload("results/test.csv", path = drive)

#create a species folder
folder <- drive_mkdir(name = "species_test", path = drive_results)

#list all files in local results folder (set working directory to results)
local_files <- list.files("results/", recursive = TRUE, full.names = TRUE)
 
# upload all files to species folder (this takes a while)
with_drive_quiet(
  files <- purrr::map(local_files, ~ drive_upload(.x, path = folder))
)
                     