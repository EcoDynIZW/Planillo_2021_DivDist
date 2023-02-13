## 26 10 2018

## JSDM invertebrates. ENVIRONMENTAL MODEL. species in 3 or more sites

# CREATING MAPS FOR VISUALIZING OBSERVED AND PREDICTED DATA

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Model thin = 100
## Samples: 30000
## Chains: 3

rm(list = ls())

getMKLthreads()
setMKLthreads(4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(devtools)
#install Rtools to be able to install BayesLogit unless you have it already
# library(httr)
# set_config(use_proxy(url="http://192.168.2.2", port=3128, username="user",password="password"))
# 
# install_url('https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz')
# library(BayesLogit)
# install_github("gtikhonov/HMSC")
library(Hmsc)
# library(stats)
library(abind)
# library(corrplot)
library(ggplot2)
# library(nnet)
# library(coda)
# library(ape)
# library(phytools)
# library(fields)
# library(dplyr)
# MARKOV CHAINS
# devtools::install_github("njtierney/mmcc")
# library(mmcc)
# library(MCMCvis)
# MAPPING
library (sf)
library(tmap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workspace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#WorkDir <- "D:/Dropbox/_IZW/APlanilloIn2018_BIBSJSDMs"

WorkDir <- getwd()

Data_wd <- file.path(WorkDir, "2_Raw_data")
output_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_3sites/")
maps_wd <- file.path(WorkDir, "2_Raw_data/ForPrediction/maps")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(WorkDir, "/3_Scripts/Functions_MapModels.R"))


#####################
# Mapping

## Plotting observed and predicted data in bibs plots

# Survey data
plots <- read_sf(paste0(maps_wd, "/Dry_Grassland_Plots_Berlin_25833.gpkg"))
plot(st_geometry(plots), main = "BIBS plots")

# Environment
## berlin and land use
berlin <- read_sf(paste0(maps_wd, "/berlin_city_border_25833.gpkg"))
berlin.out <- read_sf(paste0(maps_wd,"/Berlin_border.gpkg"))

berlin.out <- berlin.out %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin.out <- st_transform(berlin.out, crs = 25833)

water <- read_sf(paste0(maps_wd, "/waterbodies_Berlin_25833.gpkg"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Formulas for each model

# regression formulas for env. covariates

# XFormula.bees1 = ~ patch.area + imperv.500m + temp.500m + open.green.500m
# XFormula.bees2 = ~ imperv.500m + open.green.500m
# 
# XFormula.carabids1 = ~ plant.abund + imperv.100m + noise.100m + dist.water.100m + temp.100m + open.green.100m
# XFormula.carabids2 = ~ imperv.100m + dist.water.100m + temp.100m
# 
# XFormula.grasshop1 = ~ imperv.100m + noise.100m + temp.500m + open.green.500m
# XFormula.grasshop2 = ~ imperv.100m + open.green.500m
# XFormula.grasshop3 = ~ imperv.100m + noise.100m + temp.100m + open.green.100m
# XFormula.grasshop4 = ~ imperv.100m + temp.100m

# XFormula.spiders1 = ~ patch.area + imperv.100m + noise.100m + temp.100m + open.green.100m
# XFormula.spiders2 = ~ imperv.100m + noise.100m + temp.100m + open.green.100m



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIBS PLOTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Observations in the bibs plots

# Load coordinates of plots
my.xy.tmp <- read.csv(paste0(Data_wd, '/env_42plots.csv'), sep = ';')
my.xy <- my.xy.tmp[,c(1, 9:10)]

# Load species observation and plot them in a map

model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4",
                 "spiders1", "spiders2")
species <- c("bees", "bees",
             "carabids", "carabids",
             "grasshoppers", "grasshoppers", "grasshoppers", "grasshoppers",
             "spiders", "spiders")

# Loop
for (i in 1:length(species)){
  obs.sp <- species[i]
  modelType <- model.names[i]
  # Load observation matrix and remove sp only in 1 or 2 sites
  observed.sp.tmp <- read.csv(paste0(Data_wd, "/invertebrates/", obs.sp, "_sitexsp_namesphylo_42plots.csv"))
  observed.sp.tmp2 <- as.matrix(observed.sp.tmp[, -1])  
  observed.sp <- observed.sp.tmp2[,colSums(observed.sp.tmp2>0)>2] # Removing species present only in 1 or 2 sites
  ncol(observed.sp)
  
  observed.abu <- rowSums(observed.sp) # Abundance per site
  observed.s <- rowSums(observed.sp > 1) # Richness per site
  
  # Create dataframe with all info
  
  observations <- cbind.data.frame(my.xy, observed.abu, observed.s)
  write.csv(observations, paste0(output_wd, "/", modelType, '/Observed_abu_richness_', 
                   obs.sp, '.csv'), row.names = FALSE)
  # Make them spatial object
  
  observ.spatial <- st_as_sf(observations, coords = c("X", "Y"), crs = 25833)
  
  # Make and save map
  default.breaks <- c(0, round(quantile(observations$observed.abu)[2:4], 0), round(max(observations$observed.abu), 0))
  default.title <- paste0(modelType, ". \nObserved Abundance and Richness in BIBS plots")
  output <- paste0(output_wd, "/", modelType, '/Map_observed_bibsplots_Abu_S_', 
                   obs.sp, '.pdf')
  see_save_pred_map(prediction.file = observ.spatial, bubbles.size = "observed.s", bubbles.col = "observed.abu",  
                    breaks = default.breaks, main.title = default.title, 
                    print.map = TRUE, save.mymap = TRUE, output = output)
}



## Predictions for the observed data

plant.tmp <- read.csv(paste0(Data_wd, "/EnvValues/plants_42plots.csv"))
env100m.tmp <- read.csv(paste0(Data_wd, "/EnvValues/bibsplots_env_centered_42plots_100m.csv"))
env500m.tmp <- read.csv(paste0(Data_wd, "/EnvValues/bibsplots_env_centered_42plots_500m.csv"))

Xdata <- cbind.data.frame(plant.abund = plant.tmp$plants.abund,
                          patch.area = env100m.tmp$patch.area,
                          imperv.100m = env100m.tmp$impervious_surface,
                          noise.100m = env100m.tmp$noise,
                          dist.water.100m = env100m.tmp$distance_water,
                          temp.100m = env100m.tmp$temp_day,
                          open.green.100m = env100m.tmp$open_green,
                          imperv.500m = env500m.tmp$impervious_surface,
                          temp.500m = env500m.tmp$temp_day,
                          open.green.500m = env500m.tmp$open_green)
summary(Xdata)

## spatial coordinates
xyData.tmp <- read.csv(paste0(Data_wd, '/env_42plots.csv'), sep = ';')
xyData <- cbind(x = xyData.tmp$X, y = xyData.tmp$Y)


model.names <- c("bees1", "bees2",
  "carabids1", "carabids2",
  "grasshop1", "grasshop2", "grasshop3", "grasshop4")
#"spiders1", "spiders2"

for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata", sep = ""))
  load(file = filename)
  
  ## Predict in the observed plots
  predYR.all <- m$predict()  
  predYR <- apply(abind(predYR.all,along=3),c(1,2),mean)

  # Get the predicted abundance and richness
  Abu.pred <- rowSums(predYR) # In this case, abundance. Richness with pa data
  S.pred <- rowSums(predYR > 1) # richness taken as species presence when predicted abundance is higher than 1

  # SAVE PREDICTIONS
  predictions.plots <- cbind.data.frame(xyData.tmp$Plot, xyData,Abu.pred,S.pred)
  write.csv(predictions.plots, file=paste0(output_wd, "/", modelType, '/Predictions_bibsplots_', 
                                     modelType, '.csv'), row.names = FALSE)
}



# Predictions in the bibs plots
for (i in 1:length(model.names)){
  # Get predictions
  modelType <- model.names[i]
  
  predictions.bibsplots <- read.csv(paste0(output_wd, "/", modelType, '/Predictions_bibsplots_', 
                                 modelType, '.csv'))
  # Make them spatial object
  predict.spatial <- st_as_sf(predictions.bibsplots, coords = c("x", "y"), crs = 25833)
  # predict.spatial <- st_as_sf(predictions.plots, coords = c("x", "y"), crs = 25833)
  
  
  # Make and save map
  default.breaks <- c(0, round(quantile(predictions.bibsplots$Abu.pred)[2:4], 0), round(max(predictions.bibsplots$Abu.pred), 0))
  default.title <- paste0(modelType, ". \nPredicted Abundance and Richness in BIBS plots")
  output <- paste0(output_wd, "/", modelType, '/Map_predictions_bibsplots_Abu_S_', 
                   modelType, '.pdf')
  
  see_save_pred_map(prediction.file = predict.spatial, bubbles.size = "S.pred", bubbles.col = "Abu.pred",  
                    breaks = default.breaks, main.title = default.title, 
                    print.map = TRUE, save.mymap = TRUE, output = output)
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BIRD TRANSECTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# READING THE GRID DATA FOR BIRD TRANSECTS
grid.tmp100 <- read.csv(file.path(Data_wd, "ForPrediction/predict_transects_env_cs_100m.csv"))
grid.tmp500 <- read.csv(file.path(Data_wd, "ForPrediction/predict_transects_env_cs_500m.csv"))

xy.transects <- read.csv(file.path(Data_wd, "ForPrediction/transects_coordinates_25833.csv"), sep=",")

## Plant abund cannot be used for predictions because no data
## Patch area not included because is problematic for birds

grid <- as.data.frame(cbind(imperv.100m = grid.tmp100$impervious_surface,
                            noise.100m = grid.tmp100$noise,
                            dist.water.100m = grid.tmp100$distance_water,
                            temp.100m = grid.tmp100$temp_day,
                            open.green.100m = grid.tmp100$open_green,
                            imperv.500m = grid.tmp500$impervious_surface,
                            open.green.500m = grid.tmp500$open_green,
                            temp.500m = grid.tmp500$temp_day))


### PREDICT MODELS IN A LOOP IN THE BIRD TRANSECTS

# Not using carabids1 this time because no data on plant abund for the prediction, 
# Not using bees1 because patch.area
model.names <- c(#"bees1", 
  "bees2",
  # "carabids1", 
  "carabids2",
  "grasshop1", "grasshop2", "grasshop3", "grasshop4")
#"spiders1", "spiders2"

for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata", sep = ""))
  load(file = filename)

  #DEFINING THE NEW STUDY DESIGN
  # nyNew = nrow(grid)
  dfPiNew <- xy.transects$routcode_2
  dfPiNew = as.data.frame(dfPiNew)
  colnames(dfPiNew) = m$levelNames
  dfPiAll=rbind(m$dfPi,dfPiNew)

  #DEFINING RANDOM EFFECTS THAT INCLUDE BOTH OLD (USED FOR MODEL FITTING) AND NEW (THE GRID DATA) UNITS
  rL1=m$rL[[1]]$clone()
  xyold = rL1$s
  xy = xy.transects[,2:3]
  rownames(xy) = dfPiNew[,1]
  colnames(xy) = colnames(xyold)
  xyall = rbind(xyold,xy)
  rL1$pi = dfPiAll[,1]
  rL1$s = xyall

  # GET PREDICTIONS
  predYR.all <- m$predict(dfPiNew = dfPiNew, XData = grid, rL = list(rL1), expected = TRUE, predictEtaMean = TRUE)  
  predYR <- apply(abind(predYR.all,along=3),c(1,2),mean)
  
  # GET THE SUMMARY IN ABUNDANCE OR RICHNESS
  Abu =rowSums(predYR) # In this case, abundance. Richness with pa data
  S = rowSums(predYR > 1) # richness taken as species presence when predicted abundance is higher than 1
  mapData=data.frame(xy,Abu,S)
  
  # SAVE PREDICTIONS
  predictions <- cbind(dfPiNew, xy, Abu, S)
  write.csv(predictions, file=paste0(output_wd, "/", modelType, '/Predictions_transects_', 
                                     modelType, '.csv'), row.names = FALSE)
}
  


# PLOT PREDICTED ABUNDANCE AND RICHNESS
# Abuplot <- ggplot(data = mapData, aes(x=x, y=y, color=Abu)) + 
#     geom_point(size=3)
# Abuplot + ggtitle(paste0(modelType, "Predicted species abundance")) + 
#     scale_color_gradient(low="lightblue1", high="navy")
# ggsave(paste0(output_wd, "/", modelType, '/Predicted_Abundance_', 
#                 modelType, '.pdf'))
# Splot <- ggplot(data = mapData, aes(x=x, y=y, color=S)) + 
#     geom_point(size=3)
# Splot + ggtitle(paste0(modelType, "Predicted species richness")) + 
#     scale_color_gradient(low="lightsalmon", high="orangered")
# ggsave(paste0(output_wd, "/", modelType, '/Predicted_Richness_', 
#                 modelType, '.pdf'))
#   
  



#####################
# Mapping

# Survey data
centroids <- read_sf(paste0(maps_wd, "/birds_points_25833.gpkg"))
plot(st_geometry(centroids), main = "Bird centroids")

# Environment
## berlin and land use
berlin <- read_sf(paste0(maps_wd, "/berlin_city_border_25833.gpkg"))
berlin.out <- read_sf(paste0(maps_wd,"/Berlin_border.gpkg"))

berlin.out <- berlin.out %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin.out <- st_transform(berlin.out, crs = 25833)

water <- read_sf(paste0(maps_wd, "/waterbodies_Berlin_25833.gpkg"))

# Predictions in the bird transects
for (i in 1:length(model.names)){
  # Get predictions
  modelType <- model.names[i]
  
  predictions <- read.csv(paste0(output_wd, "/", modelType, '/Predictions_transects_', 
                          modelType, '.csv'))
  # Make them spatial object
  predict.spatial <- st_as_sf(predictions, coords = c("x", "y"), crs = 25833)
  
  # Make and save map
  default.breaks <- c(0, round(quantile(predictions$Abu)[2:4], 0), round(max(predictions$Abu), 0))
  default.title <- paste0(modelType, ". \nPredicted Abundance and Richness in bird transects")
  output <- paste0(output_wd, "/", modelType, '/Map_predictions_transects_Abu_S_', 
         modelType, '.pdf')
  
  
  see_save_pred_map(prediction.file = predict.spatial, bubbles.size = "S", bubbles.col = "Abu",  
                    breaks = default.breaks, main.title = default.title, 
                    print.map = FALSE, save.mymap = TRUE, output = output)
}







################################################################################
################################################################################

# END

################################################################################
################################################################################





#####################
# MAPPING

# library (raster)
library(sp)
library(rgdal)
library (sf)
library(tmap)

maps_wd <- file.path(WorkDir, "2_Raw_data/ForPrediction/maps")



#####################
# Predictions

pred.spiders <- SpatialPointsDataFrame(prediction.spiders[,2:3], prediction.spiders)

#####################
# Survey data

plots <- read_sf(paste0(maps_wd, "/Dry_Grassland_Plots_Berlin_25833.gpkg"))
plot(st_geometry(plots), main = "BIBS plots")


pred.spiders <- st_as_sf(pred.spiders)
pred.spiders <- pred.spiders %>% st_set_crs(25833)



tmap_mode("plot")
tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  tm_shape(water) +
  tm_polygons("blue") +
  tm_shape(pred.spiders) +
  tm_bubbles("pred.rich.spiders", col = "pred.abu.spiders", 
             border.col = "black", border.alpha = .5, 
             style="fixed", breaks=c(seq(0, 200, by=50), Inf),
             palette="-RdYlBu", contrast=1, 
             title.size="Species richness", 
             title.col="Species abundance") +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(compass.type = "4star") +
  tm_compass(position = c("right", "top"))+
  tm_legend(frame = TRUE,
            bg.color="lightyellow")


