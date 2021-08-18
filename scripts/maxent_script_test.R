# Test run script for running mass Maxent models on Pegasus

# SET UP --------------------------------------------------------------------

library(rgdal)
library(rJava)
options(java.parameters = "-Xmx12g" ) #increase this? depending on how much memory requested
library(dismo)
library(tidyverse) 
library(janitor) #need this?
library(stringr)

rasterOptions(tmpdir = 'temp/', progress = "text", maxmemory = 12e9) #edit this

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
bioclim <- stack(lapply(list.files("D:/endophyteDimensions/data/wc2-5/", pattern = ".bil$", full.names = TRUE), raster))
# writeRaster(bioclim, filename = "D:/endophyteDimensions/data/bioclim_5km.tif",
#             options="INTERLEAVE=BAND", overwrite = TRUE)
save(bioclim, file = "D:/endophyteDimensions/data/bioclim.RData")


load("D:/endophyteDimensions/data/bioclim.RData")


# read in occurrences (convert to spatial based on clim params)

load("path to occ Rdata file") # this is called thin.occ.update

#filter to those with #occ >= 10

model.species <- thin.occ.update %>% group_by(species) %>% count() %>% 
  filter(n >= 10)

model.occ <- thin.occ.update %>% filter(species %in% model.species$species) %>% 
  filter(species != "Cerastium glomeratum")
#remove one sp. that couldn't thin? Cerastium glomeratum'
## save this! then don't need to do this filtering step
plant.occ <- model.occ
save(plant.occ, file = "D:/endophyteDimensions/data/thinned_occ_final.RData")

occ.sp <-SpatialPointsDataFrame(model.occ[,c('decimallongitude','decimallatitude'),], 
                                data = model.occ) 
proj4string(occ.sp) <- CRS("+proj=longlat + ellps=WGS84")

#get list of species names
specs <- plant.occ %>% pull(species) %>% unique()

#test species = Panicum chillagoanum and Dubautia ciliolata
specs.test <- specs[specs %in% c('Panicum chillagoanum', "Dubautia ciliolata")]

# read in buffer distance file
## make sure all specs are in this file

buffer_dist <- read_csv("data/buffer_distances_update.csv")


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
final.summary <- data.frame(array(dim = c(length(specs), 90)))
names(final.summary) <-
  c(
    'species',
    'no.pres',
    'no.background',
    'train.AUC',
    'test.AUC',
    'maxss.threshold',
    'geo.ext.area.sqkm',
    'range.size.sqkm',
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


baseDir <- "/scratch/projects/endophytedim/model_results"

  
# for each species:
for (i in 1:length(specs)) {
  
  ## make species folder (say working in folder called "model results")
  
  setwd(baseDir)
  
  dir.create(paste(specs[i]))
  
  #folder <- list.dirs(path = getwd(), pattern = paste(specs[i]))
  
  setwd(paste(getwd(), specs[i], sep = "/")) #now files save to species folder
  
  
  occ <- subset(occ.sp, occ.sp$species == specs[i])
  
  
  # make range polygon
  dist <- as.numeric(buffer_dist[buffer_dist$species == specs[i], "distance_km"])
  
  
  range <- buffer(occ, width = (dist*1000), dissolve = TRUE) #uses meters

  #clip predictors to range polygon
  preds <- crop(bioclim, extent(range)) %>% mask(range)

  # clip all occurrences to range polygon for background points
  background <- occ.sp[!is.na(over(occ.sp, range)),]

  
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
        p = coordinates(occ),
        a = coordinates(background),
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
    
    
    
    print(k)
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
  dir.create("bestbeta") #still working in species folder
  setwd("bestbeta")
 
  
  xm <- maxent(x = preds, p = occ, a = background, args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"outputformat=logistic", "responsecurves=true", "replicates=10"), path = getwd())
  
  
  save(xm, file = paste(specs.test[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  #write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  #predict model onto native range with logistic output and save raster file
  suitability <- predict(xm, preds, filename = paste(specs[i],'suitability.tif',sep = '_'), overwrite =T)
  
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
  
  # read in maxent results file, still in bestbeta folder
  
  max.results <- read.csv("maxentResults.csv") %>% filter(Species == "species (average)")
  
  
  final.summary[i, 1:2] <- bestbetas[i, 1:2]
  final.summary[i,3:5] <- max.results[c("X.Background.points", "Training.AUC", "Test.AUC")]
  final.summary[i,6] <- max.results["Maximum.test.sensitivity.plus.specificity.Logistic.threshold"]
 
  #calculate area of geo extent (using preds raster to remove water area)
  extent_area <- area(preds[[1]], na.rm = TRUE)
  extent_area <- extent_area[!is.na(extent_area)]
  final.summary[i,7] <- sum(extent_area)
  
  # get area of suitable habitat (i.e., range size) using maxss threshold
  native.range.suit <- native.range.average # to preserve elements of native.range.average map
  
  native.range.suit[native.range.suit < final.summary[i, "maxss.threshold"]] <- NA
  suit_area <- area(native.range.suit, na.rm = TRUE)
  suit_area <- suit_area[!is.na(suit_area)]
  final.summary[i,8] <- sum(suit_area)
    
    #get most important variable and extract values @ all points
  
    
    ## contribution  
    final.summary[i,9] <- max.results %>% 
      dplyr::select(Species, contains("contribution")) %>% 
      filter(Species == "species (average)") %>% 
      data.frame(variable = names(.[2:length(.)]), average = as.numeric(.[1,2:length(.)])) %>% 
      dplyr::select(variable, average) %>% 
      filter(average == max(average)) %>% 
      pull(variable) %>% 
      str_sub(end = -14)
    
    top.cont.values <- raster::extract(preds[[final.summary[i,9]]], occ)
    final.summary[i,10] <- max(top.cont.values)
    final.summary[i,11] <- min(top.cont.values)
    
    ## permutation imp  
    final.summary[i, 12] <- max.results %>% 
      dplyr::select(Species, contains("permutation")) %>% 
      filter(Species == "species (average)") %>% 
      data.frame(variable = names(.[2:length(.)]), average = as.numeric(.[1,2:length(.)])) %>% 
      dplyr::select(variable, average) %>% 
      filter(average == max(average)) %>% 
      pull(variable) %>% 
      str_sub(end = -24)
    
    top.perm.values <- raster::extract(preds[[final.summary[i, 12]]], occ)
    final.summary[i, 13] <- max(top.perm.values)
    final.summary[i, 14] <- min(top.perm.values)
    
    #now paste raw contribution/importance for all vars
    final.summary[i, 15:52] <- max.results[12:49]
    
    #now get max/min of all predictor vars within suitable habitat
    preds.suit <- preds[!is.na(native.range.suit)] %>% as.data.frame()
    for (k in 1:19){
      
      final.summary[i, paste0("max.suit.bio.", k)] <- max(preds.suit[paste0("bio", k)])
      final.summary[i, paste0("min.suit.bio.", k)] <- min(preds.suit[paste0("bio", k)])

    }
    
  
  # make folder with species specific results want to keep
    dir.create("results")
    setwd("results")
    
    file.copy("maxentResults.csv", "results/")
    save(xm, file = paste(paste0("results/", specs.test[i]), bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
    write.csv(allbetas, paste0("results/", specs.test[i], "_allbetas.csv"))
    #save the next two in case something happens and final file is not created at end of species loop
    write.csv(bestbetas[i,], paste0("results/", specs.test[i], "_beta_params.csv"))
    write.csv(final.summary[i,], paste0("results/", specs.test[i], "final_summary.csv"))
    
    writeRaster(x = maps, filename = paste0("results/", specs.test[i], '_map_outputs.tif'), 
                options="INTERLEAVE=BAND", overwrite = TRUE)
    
    
    system("rclone copy [path to results folder] endodim:")
  
  
    #remove folder
    rm("path of folder name")
  
    
    print(i)
    
  
}

## AFTER LOOP file saving


#bestbetas and model params
write.csv(bestbetas, file = 'bestbetas.csv') # give specific batch name

#final summary file
write.csv(final.summary, "path to folder/final.summary") #add specific batch name

# remove files
removeTmpFiles(h=0)
