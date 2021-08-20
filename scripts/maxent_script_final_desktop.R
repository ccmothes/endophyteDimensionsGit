# Script for running models on Desktop comp

# SET UP --------------------------------------------------------------------

#create package load function to install (if necessary) and load all required packages
packageLoad <-
  function(x) {
    for (i in 1:length(x)) {
      if (!x[i] %in% installed.packages()) {
        install.packages(x[i])
      }
      library(x[i], character.only = TRUE)
    }
  }

#install and/or load packages
packageLoad(
  c(
    "rgdal",
    "rgeos",
    "rJava",
    "dismo",
    "dplyr",
    "purrr",
    "readr",
    "stringr",
    "googledrive"
  )
)

options(java.parameters = "-Xmx12g" ) #may want to change this

#if changing java params above, match maxmemory here too
rasterOptions(progress = "text", maxmemory = 12e9)  


#set path to google drive folder for file saving later
drive_results <- "https://drive.google.com/drive/u/0/folders/1osfjlJP5dCoBp6COXMveidomPZPD_xGi"


#create raster averaging function (at least twice as fast as 'calc')
rasterstack_mean <- function(x) {
  s0 <- nlayers(x)
  s1 <- raster(x, layer=1)
  for(ri in 2:s0) {
    r <- raster(x, layer=ri)
    s1 <- s1 + r
  }
  return(s1/s0)
}




# FILE PREP ----------------------------------------------

# read in global predictors
bioclim <- brick("data/bioclim_prj_5km.tif")

# need to add layer names, not preserved when raster was written
names(bioclim) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14",
                    "bio15", "bio16", "bio17", "bio18", "bio19", "bio2",
                    "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")


# read in occurrences (and convert to spatial based on bioclim params)
load("data/plant.occ.final.RData") # this loads object 'plant.occ'

occ.sp <-SpatialPointsDataFrame(plant.occ[,c('decimallongitude','decimallatitude'),], 
                                data = plant.occ) 

proj4string(occ.sp) <- CRS("+init=epsg:4326")

#create set of planar coordinates for buffer calculation
occ.planar <-
  spTransform(occ.sp, CRSobj = '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


#get list of species names
specs <- read.csv("data/species_list.csv") %>% pull(species)


# read in buffer distance file
buffer_dist <- read.csv("data/buffer_dist_model.csv")


##make sequence of betas to use
betas <- array(dim = c(24, 2))
betas[,1] <- "betamultiplier="
betas[,2] <- c(seq(0.2,1, 0.2), seq(2, 20, 1))
betas <- data.frame(betas)
betas[,3] <- paste(betas[,1], betas[,2],sep = '')
names(betas)[3] <- 'setbeta'
betas[,2] <- as.character(betas[,2])

#make output data frame for beta results
bestbetas <- data.frame(array(dim = c(length(specs), 5)))
names(bestbetas) <- c('species','no.pres','no.params default','no.params best','bestbeta')


#make output data frame for all summary results
final.summary <- data.frame(array(dim = c(length(specs), 94)))
names(final.summary) <-
  c(
    'species',
    'no.pres',
    'no.background',
    'train.AUC',
    'test.AUC',
    'geo.ext.area.sqkm',
    'maxss.threshold',
    'MTP.threshold',
    'MTP.10.percent.threshold',
    'maxss.range.size.sqkm',
    'MTP.range.size.sqkm',
    'MTP10.range.size.sqkm',
    'top.var.contribution',
    'max.top.contribution',
    'min.top.contribution',
    'top.var.importance',
    'max.top.importance',
    'min.top.importance',
    paste(paste0("bio", c(1, 10:19, 2:9)), "contribution", sep = "."), #to match order in maxentResults files
    paste(paste0("bio", c(1, 10:19, 2:9)), "permutation.importance", sep = "."),
    paste("max.suit.bio", 1:19, sep = "."),
    paste("min.suit.bio", 1:19, sep = ".")
  )


##### EDIT THIS ############ 
#batch number is the range of index #'s in species file you want to run
batch <- 


# create batch specific folder (in case running multiple sessions on same comp)
dir.create(paste("batch", batch[1], batch[length(batch)], sep = "_"))  

baseDir <-
  paste0(getwd(),
    paste("/batch",
    batch[1],
    batch[length(batch)],
    sep = "_"
  ))
  
setwd(baseDir) #now working in batch specific directory
dir.create("results")
  
# for each species:
for (i in batch[1]:batch[length(batch)]) {
  tryCatch({
  
  setwd(baseDir)
  
  #remove folder from previous species run
 
  unlink("allbetas", recursive = TRUE)
  
  
  dir.create("allbetas")
  
  #folder <- list.dirs(path = getwd(), pattern = paste(specs[i]))
  
  setwd("allbetas") #now files save in this folder
  
  
  occ <- subset(occ.sp, occ.sp$species == specs[i])
  
  occ.p <- subset(occ.planar, occ.planar$species == specs[i])
  
  # make range polygon
  dist <- as.numeric(buffer_dist[buffer_dist$species == specs[i], "distance_km"])
  
  
  #switch to planar to avoid invalid geometry
  range <- buffer(occ.p, width = (dist*1000), dissolve = TRUE) #uses meters

  
  #clip predictors to range polygon
  preds <- crop(bioclim, extent(range)) %>% mask(range)
  
  
  #convert back to lat/long
  preds <- projectRaster(preds, crs = crs(occ.sp))
  
  
  # clip all occurrences to range polygon for background points
  background <- occ.planar[!is.na(over(occ.planar, range)),]
  
  
  #sample just 10,000 occ if over
  if (nrow(background) > 10000)
    background <- background[sample(nrow(background),10000),]
  
  #convert background back to lat/long
  background <- spTransform(background, CRSobj = crs(occ.sp))
  
  
  allbetas <- data.frame(array(dim = c(24, 7)))
  names(allbetas) <-
    c('species',
      'beta',
      'no.params',
      'no.pres',
      'loglik',
      'AIC',
      'AICc')
  allbetas[, 1] <- specs[i]
  
  
  for (k in 1:nrow(betas)) {
    
    xm <-
      maxent(
        x = preds,
        p = occ,
        a = background,
        args = c(
          betas[k, "setbeta"],
          "responsecurves=TRUE",
          "writebackgroundpredictions=TRUE"
        ),
        path = getwd()
      )
    
    #save(xm, file = paste(specs[i], betas[k, "X2"], "xm.RData", sep = '_'))
    
    #get lambdas info out of xm
    lambdas <- data.frame(strsplit(xm@lambdas, ','))
    lambdas <- t(lambdas)
    lambdas <- data.frame(lambdas)
    rownames(lambdas) <- 1:nrow(lambdas)
    #vitalconstants <- lambdas[(nrow(lambdas)-3):nrow(lambdas),1:2]
    lambdas <- lambdas[-((nrow(lambdas) - 3):nrow(lambdas)), ]
    colnames(lambdas) <- c('variable', 'lambda_estimate', 'min', 'max')
    lambdas$lambda_estimate <-
      as.numeric(as.character(lambdas$lambda_estimate))
    
    #read and write predictions at presences
    pres.raw <- read.csv('species_samplePredictions.csv')
    write.csv(pres.raw,
              paste(specs[i], betas[k, "X2"], "pres.predictions.csv", sep = '_'))
    
    #calculate AICc
    no.params <- length(which(lambdas$lambda_estimate != 0))
    no.pres <- xm@results[1]
    loglik <- sum(log(pres.raw$Raw.prediction))
    AIC <- (2 * no.params) - (2 * loglik)
    AICc <-
      AIC + (((2 * no.params) * (no.params + 1)) / (no.pres - no.params - 1))
    
    
    allbetas[k, 2:7] <-
      c(betas[k, "X2"], no.params, no.pres, loglik, AIC, AICc)
    

  }
  
  # need to save this?
  write.csv(allbetas, file = paste(specs[i], "allbetas.csv", sep = "_"))
  
  
  # find best beta and save in bestbetas file
  
  bestbetas[i, 1] <- as.character(allbetas$species[1])
  bestbetas[i, 2] <- allbetas$no.pres[1]
  bestbetas[i, 3] <- allbetas[which(allbetas$beta == 1), 'no.params']
  
  if (max(allbetas$no.pres) > min(allbetas$no.params))
    bestbetas[i, 5] <-
    min(allbetas[which(allbetas$AICc == min(allbetas[which(allbetas$no.params < allbetas$no.pres &
                                                             allbetas$no.params > 0), 'AICc'])), 'beta'])
  
  if ((max(allbetas$no.pres) > min(allbetas$no.params)))
    bestbetas[i, 4] <-
    allbetas[which(allbetas$beta == bestbetas[i, 5]), 'no.params']
  
  
  # BEST BETA MODEL
  
  #set new directory for current best beta
  dir.create("bestbeta") #still working in allbetas folder
  setwd("bestbeta")
  
  
  xm <- maxent(x = preds, p = occ, a = background, args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"outputformat=logistic", "responsecurves=true", "replicates=10"), path = getwd())
  
  
  save(xm, file = paste(specs[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  #write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  #predict model onto native range with logistic output and save raster file
  suitability <- predict(xm, preds)
  
  # #this is a stack for each rep, get single average raster
  # native.range.average <- calc(native.range, fun = mean, na.rm = T)
  # #try faster function
  # rasterstack_mean <- function(x) {
  #   s0 <- nlayers(x)
  #   s1 <- raster(x, layer=1)
  #   for(ri in 2:s0) {
  #     r <- raster(x, layer=ri)
  #     s1 <- s1 + r
  #   }
  #   return(s1/s0)
  # }
  # system.time({rasterstack_mean(native.range)}) #this is about twice as fast
  # system.time({calc(native.range, fun = mean, na.rm = T)})
  
  # Need average of all 10 files (for each cv run)
  suit.average <- rasterstack_mean(suitability)
  
  names(suit.average) <- paste(specs[i], "habitat_suitability", sep = "_")
  
  #stack with predictors
  maps <- stack(preds, suit.average)
  
  
  # start filling in final summary -----------------------------------------------
  
  # read in maxent results file
  
  max.results <- read.csv("maxentResults.csv") %>% filter(Species == "species (average)")
  
  
  final.summary[i, 1:2] <- bestbetas[i, 1:2]
  final.summary[i,3:5] <- max.results[c("X.Background.points", "Training.AUC", "Test.AUC")]
  
  #calculate area of geo extent (using preds raster to remove water area)
  extent_area <- area(preds[[1]], na.rm = TRUE)
  extent_area <- extent_area[!is.na(extent_area)]
  final.summary[i,6] <- sum(extent_area)
  
  final.summary[i,7:9] <- max.results[c("Maximum.test.sensitivity.plus.specificity.Logistic.threshold",
                                      "Minimum.training.presence.Logistic.threshold",
                                      "X10.percentile.training.presence.Logistic.threshold")]
 
  
  
  # get area of suitable habitat (i.e., range size) using the three different thresholds
  suit.average.maxss <- suit.average # to preserve elements of suitability map
  
  suit.average.maxss[suit.average.maxss < final.summary[i, "maxss.threshold"]] <- NA
  suit_area <- area(suit.average.maxss, na.rm = TRUE)
  suit_area <- suit_area[!is.na(suit_area)]
  final.summary[i,10] <- sum(suit_area)
  
  suit.average.mtp <- suit.average # to preserve elements of suitability map
  
  suit.average.mtp[suit.average.mtp < final.summary[i, "MTP.threshold"]] <- NA
  suit_area <- area(suit.average.mtp, na.rm = TRUE)
  suit_area <- suit_area[!is.na(suit_area)]
  final.summary[i,11] <- sum(suit_area)
  
  suit.average.mtp10 <- suit.average # to preserve elements of suitability map
  
  suit.average.mtp10[suit.average.mtp10 < final.summary[i, "MTP.10.percent.threshold"]] <- NA
  suit_area <- area(suit.average.mtp10, na.rm = TRUE)
  suit_area <- suit_area[!is.na(suit_area)]
  final.summary[i,12] <- sum(suit_area)
    
    
  #get most important variable and extract values @ all points
  
    
    ## contribution  
    final.summary[i,13] <- max.results %>% 
      dplyr::select(Species, contains("contribution")) %>% 
      filter(Species == "species (average)") %>% 
      data.frame(variable = names(.[2:length(.)]), average = as.numeric(.[1,2:length(.)])) %>% 
      dplyr::select(variable, average) %>% 
      filter(average == max(average)) %>% 
      pull(variable) %>% 
      str_sub(end = -14)
    
    top.cont.values <- raster::extract(preds[[final.summary[i,13]]], occ)
    final.summary[i,14] <- max(top.cont.values, na.rm = TRUE)
    final.summary[i,15] <- min(top.cont.values, na.rm = TRUE)
    
    ## permutation imp  
    final.summary[i, 16] <- max.results %>% 
      dplyr::select(Species, contains("permutation")) %>% 
      filter(Species == "species (average)") %>% 
      data.frame(variable = names(.[2:length(.)]), average = as.numeric(.[1,2:length(.)])) %>% 
      dplyr::select(variable, average) %>% 
      filter(average == max(average)) %>% 
      pull(variable) %>% 
      str_sub(end = -24)
    
    top.perm.values <- raster::extract(preds[[final.summary[i, 16]]], occ)
    final.summary[i, 17] <- max(top.perm.values, na.rm = TRUE)
    final.summary[i, 18] <- min(top.perm.values, na.rm = TRUE)
    
    #now paste raw contribution/importance for all vars
    final.summary[i, 19:56] <- max.results[12:49]
    
    #now get max/min of all predictor vars within suitable habitat using maxss
    preds.suit <- preds[!is.na(suit.average.maxss)] %>% as.data.frame()
    for (k in 1:19){
      
      final.summary[i, paste0("max.suit.bio.", k)] <- max(preds.suit[paste0("bio", k)], na.rm = TRUE)
      final.summary[i, paste0("min.suit.bio.", k)] <- min(preds.suit[paste0("bio", k)], na.rm = TRUE)

    }
    
  
  ## FILE SAVING --------------------------------------
   
    #save model outputs in species folder
    dir.create(paste0(baseDir,"/results/", specs[i])) 
    
    file.rename(from = "maxentResults.csv", to = paste0(baseDir, "/results/", specs[i], "/", specs[i], "_maxentResults.csv"))
    write.csv(allbetas, paste0(baseDir, "/results/",  specs[i], "/", specs[i], "_allbetas.csv"))
    
    
    
    file.rename(from ="maxent.html", to = paste0(baseDir, "/results/", specs[i], "/", specs[i], "_maxent.html"))
    file.rename(from = "plots/", to = paste0(baseDir, "/results/", specs[i], "/", "plots/"))
    # remove all the extra plots
    files <- list.files(path = paste0(baseDir, "/results/", specs[i], "/", "plots/"), full.names = TRUE)
    files.remove <- grep(pattern = "species_bio", x = files, value = TRUE, invert = TRUE)
    file.remove(files.remove)
    
    save(xm, file = paste(paste0(baseDir, "/results/", specs[i], "/", specs[i]), bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
    
    #save the next two in case something happens and final file is not created at end of species loop
    write.csv(bestbetas[i,], paste0(baseDir, "/results/",  specs[i], "/", specs[i], "_beta_params.csv"))
    write.csv(final.summary[i,], paste0(baseDir, "/results/",  specs[i], "/", specs[i], "_final_summary.csv"))
    
    writeRaster(x = maps, filename = paste0(baseDir, "/results/",  specs[i], "/", specs[i], '_map_outputs.tif'), 
                options="INTERLEAVE=BAND", overwrite = TRUE)
    
    
   # upload to species folder on googledrive
    #create a species folder
    folder <- drive_mkdir(name = specs[i], path = drive_results)
    
    #list all files in local species results folder
    local_files <- list.files(paste0(baseDir,"/results/", specs[i], "/"), recursive = TRUE, full.names = TRUE)
    
    # upload all files to species folder
    with_drive_quiet(
      files <- purrr::map(local_files, ~ drive_upload(.x, path = folder))
    )
    
    
    print(i)
    
  }, error=function(e){cat("Species=", specs[i], "Index =", i, 
                           "ERROR :", conditionMessage(e), "\n")})
  
}

## AFTER LOOP file saving, working in endophytedim folder
setwd(paste0(baseDir, "/results/"))

#bestbetas and model params
write.csv(bestbetas, file =  paste(batch[1], batch[length(batch)], 'bestbetas.csv', sep = "_")) # give specific batch name

#final summary file
write.csv(final.summary, paste(batch[1], batch[length(batch)], "final_summary.csv", sep = "_")) #add specific batch name

#upload to google drive
drive_upload(paste0(
  baseDir,
  "/results/",
  paste(batch[1], batch[length(batch)], "final_summary.csv", sep = "_")
), path = drive_results)

drive_upload(paste0(
  baseDir,
  "/results/",
  paste(batch[1], batch[length(batch)], 'bestbetas.csv', sep = "_")
), path = drive_results)



# remove temp files
removeTmpFiles(h=0)

