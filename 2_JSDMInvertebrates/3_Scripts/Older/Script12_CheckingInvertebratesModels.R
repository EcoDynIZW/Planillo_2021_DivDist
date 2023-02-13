## 23 10 2018

## JSDM invertebrates. ENVIRONMENTAL MODEL. species in 3 or more sites

# No phylogeny

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Model thin = 100
## Samples: 30000
## Chains: 3

rm(list = ls())
gc()

getMKLthreads()
setMKLthreads(4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
#install Rtools to be able to install BayesLogit unless you have it already
# library(httr)
# set_config(use_proxy(url="http://192.168.2.2", port=3128, username="user",password="password"))
# 
# install_url('https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz')
library(BayesLogit)
# install_github("gtikhonov/HMSC")
library(Hmsc)
library(stats)
library(abind)
library(corrplot)
library(ggplot2)
library(nnet)
library(coda)
library(ape)
library(phytools)
library(fields)
library(dplyr)
# MARKOV CHAINS
# devtools::install_github("njtierney/mmcc")
library(mmcc)
library(MCMCvis)
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



# predict_spiders <- file.path(WorkDir, "4_Processed_data/InvPredictions/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(WorkDir, "/3_Scripts/Functions_CheckModels.R"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop to go through all the models
# checking posteriors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load models
# gc()
model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4",
                 "spiders1", "spiders2")
# model.names <- "spiders2"


for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]

  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                      "_abu_3sites.Rdata", sep = ""))
  load(file = filename)

  #Evaluate model fit: mixing 
  #Extract the posteior as a coda object
  mpost <- m$convertToCodaObject()

  # Assess the convergence of the sampling
  # as we are evaluating the null model, it only has intercept in the fixed regression part
  assess.diagn(c.obj = mpost$Gamma, output = paste0(output_wd, "/", modelType), par = 'Gamma')
  
  
  # Examining posterior distributions numerically
  # environmental effects
  # tidy(mpost$Gamma)
  # glance(mpost$Gamma)
  # summary(mpost)
  summary(mpost$Gamma)
  ESS(mpost$Gamma)
  Gamma.results <- as.data.frame(MCMCsummary(mpost$Gamma))
  Beta.results <- as.data.frame(MCMCsummary(mpost$Beta))
  write.csv(Gamma.results, paste0(output_wd, "/", modelType, "/posterior_gamma_environment_", modelType, ".csv"))
  write.csv(Beta.results, paste0(output_wd, "/", modelType, "/posterior_beta_species_", modelType, ".csv"))

  # Markov chains
  MCMCtrace(mpost$Gamma, 
            pdf = TRUE, 
            open_pdf = FALSE,
            filename = paste0("MCMCtrace_", modelType),
            wd = paste0(output_wd, "/", modelType, "/"))

  par(mfrow=c(1,1))

  # Coef plot for gamma
  pdf(paste0(output_wd, "/", modelType, '/Gammas_coef_plot_', 
             modelType, '.pdf'))
  MCMCplot(mpost$Gamma, ref_ovl = TRUE)
  dev.off()

  # mpost$Beta[[1]]

  # Coef plot for betas for each species
  pdf(paste0(output_wd, "/", modelType, '/Betas_coef_plot_', 
             modelType, '.pdf'))
  MCMCplot(mpost$Beta, 
           ref_ovl = TRUE,
           rank = T,
           xlab = 'ESTIMATE',
           labels_sz = 0.3,
           med_sz = 1,
           thick_sz = 1,
           thin_sz = 1,
           ax_sz = 1,
           main_text_sz = 1)
  dev.off()

# Print a plot for each predictor
  plot.cov.betas(m) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT EFFECTS OF COVARIATES 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4",
                 "spiders1", "spiders2")
# model.names <- "spiders2"


for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste0("/model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata"))
  load(file = filename)

  postBeta = m$getPostEstimate(parName = "Beta")
  pdf(paste0(output_wd,"/", model.names[i], '/ParEstimates_Beta_model_', 
             model.names[i], '.pdf'))
  m$plotBeta(post = postBeta, param = 'Support', supportLevel = 0.95, plotTree = F, 
             spNamesNumbers = c(T,F))
  dev.off()
  
  postGamma = m$getPostEstimate(parName = "Gamma")
  pdf(paste0(output_wd, "/", model.names[i], '/ParEstimates_Gamma_model_', 
             model.names[i], '.pdf'))
  m$plotGamma(post = postGamma, param = "Support", supportLevel = 0.95)
  dev.off()
  
  # latent vars
  # postEta <- m$getPostEstimate(parName = 'Eta')
  # plot(postEta$mean[,1], postEta$mean[,2])
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop to go through all the models
# Evaluate model fit: explanatory power 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4", 
                 "spiders1", "spiders2")

for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata", sep = ""))
  load(file = filename)
  
  # Explanatory R2
  predY <- m$computePredictedValues()## returns the fitted values of models
  write.csv(predY, paste0(output_wd, "/", modelType, '/Predicted_values_sp_', 
                          modelType, '.csv'))
  # We compute the r2 manually
  R2.sp <- matrix(NA, m$ns, 1)
  for (i in 1:m$ns) {
    R2.sp[i, ] <- cor(predY[, i],m$Y[, i])^2
  }

  mean(R2.sp)
  pdf(paste0(output_wd, "/", modelType, '/Plot_R2_species_', 
             modelType, '.pdf'))
  plot(R2.sp,
       main = paste0(modelType, " - Explanatory R2 species \n R2 = ", round(mean(R2.sp), 2)),
       xlab = "Species index", 
       ylab = "R2 for each species")
  dev.off()

  R2.site <- matrix(NA, m$ny, 1)
  for (i in 1:m$ny) {
    R2.site[i, ] <- cor(predY[i, ], m$Y[i, ])^2
  }
  mean(R2.site)
  pdf(paste0(output_wd, "/", modelType, '/Plot_R2_sites_', 
             modelType, '.pdf'))
  plot(R2.site,
       main = paste0(modelType, " - Explanatory R2 sites \n R2 = ", round(mean(R2.site), 2)),
       xlab = "Site index", 
       ylab = "R2 for each site")
  dev.off()
  
  pdf(paste0(output_wd, "/", modelType, '/Plot_R2_avg_sp_abund_', 
             modelType, '.pdf'))
  plot(colSums(m$Y) / m$ny, R2.sp,
       main = paste0(modelType, " - Explanatory R2: average species abundance per site \n Mean = ", round(mean(R2.sp, na.rm = TRUE), 2)),
       ylim = c(0,1),
       xlab = "Avg Sp Abundance / Site",
       pch = 16)
  dev.off()

  #compare predicted and observed row sums (one point is a sampling unit)
  pdf(paste0(output_wd, "/", modelType, '/Plot_Obs_vs_Pred_sites_', 
             modelType, '.pdf'))
  plot(rowSums(m$Y), rowSums(predY),
       main = paste0(modelType, " - Observed vs Predicted Site abundance \n Mean R2 = ", round(mean(R2.site, na.rm = TRUE), 2)),
       xlim = c(min(rowSums(m$Y), rowSums(predY)), max(rowSums(m$Y), rowSums(predY))),
       ylim = c(min(rowSums(m$Y), rowSums(predY)), max(rowSums(m$Y), rowSums(predY))),
       ylab = "predicted abundance in sites",
       xlab = "observed abundance in sites",
       pch = 16)
  abline(0,1, col = "red", 
         lty = 4,
       lwd = 1)
  dev.off()

  
  # compare predicted and observed col sums (one data point is a species)
  pdf(paste0(output_wd, "/", modelType, '/Plot_Obs_vs_Pred_sp_', 
             modelType, '.pdf'))
  plot(colSums(predY)~colSums(m$Y),
       main = paste0(modelType, " - Observed vs Predicted Species abundance \n Mean R2 = ", round(mean(R2.sp, na.rm = TRUE), 2)),
       xlab = "observed total abundance",
       ylab = "predicted total abundance",
       pch = 16)
  abline(0,1, col = "red", 
         lty = 4,
         lwd = 1)
  dev.off()
}


  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOT EFFECTS OF COVARIATES (BEGINNING)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## PLOT SPECIES ASSOCIATIONS
# residual covariances among species

model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4")
#"spiders1", "spiders2"

# modelType <- "spiders2"

for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata", sep = ""))
  load(file = filename)
  
  # plot only those associations with 95% posterior support and order species 
  # into groups asssociated positively or negatively
  OmegaCor <- m$computeAssociations()
  supportLevel <- 0.95

  title <- as.character(m$XFormula)
  
  pdf(paste0(output_wd, "/", modelType, '/SpAssociations_', 
             modelType, '.pdf'))
  plotOrder <- corrMatOrder(OmegaCor[[1]]$mean, order = "alphabet")
  toPlot <- ((OmegaCor[[1]]$support > supportLevel) + 
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0)*OmegaCor[[1]]$mean
  corrplot(toPlot[plotOrder, plotOrder], 
           type = "lower", tl.col = "black", 
           tl.cex = 0.5, tl.srt = 45, 
           method = "square", col = colorRampPalette(c("red", "white", "green"))(200),
           # title = paste0(modelType, ". Species associations. Support level 0.95", "\n model = ", title[2], "+ Random: Site"),
           mar = c(0,0,3,0))
  mtext(paste0(modelType, ". Species associations. Support level 0.95", 
               "\n model = ", title[2], " + Random: Site"),
        at = 8.5, line = 1, cex = 1.2)
  dev.off()
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Postprocessing: parameter estimates 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


model.names <- c("bees1", "bees2",
                 "carabids1", "carabids2",
                 "grasshop1", "grasshop2", "grasshop3", "grasshop4")
#"spiders1", "spiders2"

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

bees1.groups <- c(3,1,2,2)
groupnames.patch <- c("Urban", "Green", "Patch size")
bees2.groups <- c(1,2)
groupnames.2 <- c("Urban", "Green")
carabids1.groups <- c(3,1,1,2,2,2)
groupnames.plants <- c("Urban", "Green", "Plant abund")
carabids2.groups <- c(1,2,2)
grasshop1.groups <- c(1,1,2,2)
grasshop2.groups <- c(1,2)
grasshop3.groups <- c(1,1,2,2)
grasshop4.groups <- c(1,2)
spiders1.groups <- c(3,1,1,2,2)
spiders2.groups <- c(1,1,2,2)

inv.groups <- list(bees1.groups, bees2.groups, carabids1.groups, carabids2.groups,
                   grasshop1.groups, grasshop2.groups, grasshop3.groups, grasshop4.groups,
                   spiders1.groups, spiders2.groups)
inv.groupsnames <- list(groupnames.patch, groupnames.2, groupnames.plants, groupnames.2,
                        groupnames.2, groupnames.2, groupnames.2, groupnames.2, 
                        groupnames.patch, groupnames.2)


## Variance partitioning

for (i in 1:length(model.names)){
  ## Define model to examine
  modelType = model.names[i]
  
  ## Read in the model
  filename = file.path(output_wd, paste("model_pois_", 
                                        modelType,
                                        "_abu_3sites.Rdata", sep = ""))
  load(file = filename)
  
  # Total variance explained by hte model can be partition into the contributions
  # of each fixed effect (or group) and random effect
  
  # Grouping the variables and giving names to the groups
  group <- inv.groups[[i]]
  groupnames <- inv.groupsnames[[i]]
  # compute var part
  VP <- m$computeVariancePartitioning(group = group, groupnames = groupnames)
  # plot var part
  pdf(paste0(output_wd, "/", modelType, '/Varpart_', 
             modelType, '.pdf'))
  m$plotVariancePartitioning(VP = VP)
  title(main = paste0("\n \n", modelType))
  dev.off()
}

#Order sp from more to less impact of urban variables
# VP$vals <- VP$vals[,order(VP$vals[1,], decreasing = TRUE)]
# m$plotVariancePartitioning(VP = VP)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m$covNames

#READING THE GRID DATA
grid.tmp100 <- read.csv(file.path(Data_wd, "ForPrediction/predict_transects_env_cs_100m.csv"))
grid.tmp500 <- read.csv(file.path(Data_wd, "ForPrediction/predict_transects_env_cs_500m.csv"))

xy.transects <- read.csv(file.path(Data_wd, "ForPrediction/transects_coordinates_25833.csv"), sep=",")


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

## Plant abund cannot be used for predictions because no data
grid <- as.data.frame(cbind(patch.area = grid.tmp100$patch_area,
                            imperv.100m = grid.tmp100$impervious_surface,
                            noise.100m = grid.tmp100$noise,
                            dist.water.100m = grid.tmp100$dist_wtr,
                            temp.100m = grid.tmp100$temp_day,
                            open.green.100m = grid.tmp100$open_green,
                            imperv.500m = grid.tmp500$impervious_surface,
                            open.green.500m = grid.tmp500$open_green,
                            temp.500m = grid.tmp500$temp_day))


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
  

#MODELS IN A LOOP
# Not using carabids1 this time because no data on plant abund for the prediction
model.names <- c("bees1", "bees2",
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
    
  # GET PREDICTIONS
  predYR.all <- m$predict(dfPiNew = dfPiNew, XData = grid, rL = list(rL1), expected = TRUE, predictEtaMean = TRUE)  
  predYR <- apply(abind(predYR.all,along=3),c(1,2),mean)
  

    # GET THE SUMMARY IN ABUNDANCE OR RICHNESS
  Abu =rowSums(predYR) # In this case, abundance. Richness with pa data
  S = rowSums(predYR > 1) # richness taken as species presence when predicted abundance is higher than 1
  mapData=data.frame(xy,Abu,S)
  
  # PLOT PREDICTED ABUNDANCE AND RICHNESS
  Abuplot <- ggplot(data = mapData, aes(x=x, y=y, color=Abu)) + 
    geom_point(size=3)
  Abuplot + ggtitle(paste0(modelType, "Predicted species abundance")) + 
    scale_color_gradient(low="lightblue1", high="navy")
  ggsave(paste0(output_wd, "/", modelType, '/Predicted_Abundance_', 
                     modelType, '.pdf'))
  Splot <- ggplot(data = mapData, aes(x=x, y=y, color=S)) + 
    geom_point(size=3)
  Splot + ggtitle(paste0(modelType, "Predicted species richness")) + 
    scale_color_gradient(low="lightsalmon", high="orangered")
  ggsave(paste0(output_wd, "/", modelType, '/Predicted_Richness_', 
                modelType, '.pdf'))
  
  # SAVE PREDICTIONS
  predictions <- cbind(dfPiNew, xy, Abu, S)
  write.csv(predictions, file=paste0(output_wd, "/", modelType, '/Predictions_transects_', 
                                            modelType, '.csv'), row.names = FALSE)
}



#####################
# MAPPING

# library (raster)
library (sf)
library(tmap)

maps_wd <- file.path(Data_wd, "ForPrediction/maps")



#####################
# Predictions

pred.grasshop <- st_as_sf(predictions, coords = c("x", "y"), crs = 25833)

#####################
# Survey data

plots <- read_sf(paste0(maps_wd, "/Dry_Grassland_Plots_Berlin_25833.gpkg"))
plot(st_geometry(plots), main = "BIBS plots")

centroids <- read_sf(paste0(maps_wd, "/birds_points_25833.gpkg"))
plot(st_geometry(centroids), main = "Bird centroids")


#####################
# ENVIRONMENTAL DATA

## berlin and land use

berlin <- read_sf(paste0(maps_wd, "/berlin_city_border_25833.gpkg"))
berlin.out <- read_sf(paste0(maps_wd,"/Berlin_border.gpkg"))

berlin.out <- berlin.out %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin.out <- st_transform(berlin.out, crs = 25833)

water <- read_sf(paste0(maps_wd, "/waterbodies_Berlin_25833.gpkg"))


##MAP

# pdf(paste0(output_carabids, "CarabidsPredBerlin_map.pdf"))

max(pred.grasshop$Abu)

# pdf(paste0(output_birds, "/Birds_Transects_map.pdf"))
tmap_mode("plot")
mymap <- tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  tm_shape(water) +
  tm_polygons("lightblue", border.col = NULL) +
  tm_shape(predictions) +
  tm_bubbles("S", col = "Abu", 
             border.col = "black", border.alpha = .8, 
             style="fixed", breaks = c(0, round(quantile(predictions$Abu)[2:4], 0),Inf),
             # breaks=c(-Inf, seq(20, 60, by=10), Inf),
             palette="OrRd", contrast=1, 
             title.size="Species richness", 
             title.col="Species abundance") +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(main.title = paste0(modelType, ". \nPredicted Abundance and Richness in bird transects"),
    compass.type = "4star") +
  tm_compass(position = c("right", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")
tmap_save(mymap, filename = paste0(output_wd, "/", modelType, '/Predictions_transects_Abu_S_', 
                      modelType, '.pdf'))

dev.off()



#####################################
#####################################
#####################
# Predictions

pred.birds <- SpatialPointsDataFrame(prediction.birds[,2:3], prediction.birds)

#####################
# Survey data


pred.birds <- st_as_sf(pred.birds)
pred.birds <- pred.birds %>% st_set_crs(25833)


tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  tm_shape(water) +
  tm_polygons("blue") +
  tm_shape(pred.birds) +
  tm_bubbles("pred.rich.birds", col = "pred.rich.birds", 
             border.col = "black", border.alpha = .5, 
             style="fixed", breaks=c(-Inf, seq(10, 30, by=5), Inf),
             palette="-RdYlBu", contrast=1) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(compass.type = "4star") +
  tm_compass(position = c("right", "top"))+
  tm_legend(frame = TRUE,
            bg.color="lightyellow")









##################################################################################################
#COMPUTE SPECIES RICHNESS (S), COMMUNITY WEIGHTED MEANS (predT), REGIONS OF COMMON PROFILE (RCP)
Abu =rowSums(predYR) # In this case, abundance. Richness with pa data
S = rowSums(predYR > 1) # richness taken as species presence when predicted abundance is higher than 1
predT = (predYR%*%m$Tr)/matrix(rep(Abu,m$nt),ncol=m$nt)
RCP = kmeans(predYR, 3)
RCP$cluster = as.factor(RCP$cluster)

#EXTRACT THE OCCURRENCE PROBABILITIES OF ONE EXAMPLE SPECIES
pred_species1 = predYR[,1]

##################################################################################################
#MAKE A DATAFRAME OF THE DATA TO BE PLOTTED
mapData=data.frame(xy,Abu,predT,pred_species1,RCP$cluster)

##################################################################################################
#PLOT PREDICTED OCCURRENCE OF A SINGLE SPECIES
sp <- ggplot(data = mapData, aes(x=x, y=y, color=pred_species1))+geom_point(size=3)
sp + ggtitle("Predicted species 1 occurrence") + scale_color_gradient(low="blue", high="red")


##################################################################################################
#PLOT PREDICTED ABUNDANCe
sp <- ggplot(data = mapData, aes(x=x, y=y, color=Abu))+geom_point(size=3)
sp + ggtitle("Predicted species abundance") + scale_color_gradient(low="blue", high="red")

##################################################################################################
#PLOT PREDICTED REGIONS OF COMMON PROFILE
sp <- ggplot(data = mapData, aes(x=x, y=y, color=RCP$cluster))+geom_point(size=3)
sp + ggtitle("Regions of common profile") 


## Save predictions in bird transects
pred.abu.spiders <- Abu
pred.rich.spiders <- S # richness taken as species presence when predicted abundance is higher than 1

#save predictions
prediction.spiders <- cbind(dfPiNew, xy, pred.abu.spiders, pred.rich.spiders)
# write.csv(prediction.spiders, file=paste0(predict_spiders, "/prediction_spiders_transects_env3sites_20180924.csv"), row.names = FALSE)


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


#####################
# ENVIRONMENTAL DATA

## berlin and land use

berlin <- read_sf(paste0(maps_wd, "/berlin_city_border_25833.gpkg"))
berlin.out <- read_sf(paste0(maps_wd,"/Berlin_border.gpkg"))

berlin.out <- berlin.out %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin.out <- st_transform(berlin.out, crs = 25833)

water <- read_sf(paste0(maps_wd, "/waterbodies_Berlin_25833.gpkg"))


# pdf(paste0(output_carabids, "CarabidsPredBerlin_map.pdf"))
?tm_bubbles

pdf(paste0(output_spiders, "/Spiders_plots_map.pdf"))
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
dev.off()

tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  tm_shape(water) +
  tm_polygons("blue") +
  tm_shape(pred.birds) +
  tm_bubbles("pred.rich.birds", col = "pred.rich.birds", 
             border.col = "black", border.alpha = .5, 
             style="fixed", breaks=c(-Inf, seq(10, 30, by=5), Inf),
             palette="-RdYlBu", contrast=1) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(compass.type = "4star") +
  tm_compass(position = c("right", "top"))+
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

