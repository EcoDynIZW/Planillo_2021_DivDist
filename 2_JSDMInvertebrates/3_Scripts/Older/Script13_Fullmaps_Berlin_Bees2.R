## 2 November 2018


## Prediction map for bees2 
## Used a grid of 500m

# Model: bees2
# without species only in 1 or 2 sites #

# Env predictors: imperv.500m + open.green.500m

rm(list = ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
library(BayesLogit)
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
library(sf)
library(tmap)
library(raster)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workspace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#WorkDir <- "O:/Cluster_Planillo_Aimara/JSDMs_TestingPredMaps"

WorkDir <- getwd()
Data_wd <- file.path(WorkDir, "2_Raw_data/ForPrediction")
Model_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_3sites")
Maps_wd <- file.path(WorkDir, "2_Raw_data/ForPrediction/maps")

Output_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_3sites/Maps_prediction")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the model 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename = file.path(Model_wd, "model_pois_bees2_abu_3sites.Rdata")
load(file = filename)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# See the formula of the model
m$XFormula


#READING THE GRID DATA (Covariate values extracted from a raster of 20m resolution, Berlin)
# Variable values are standardized!
grid.tmp1 <- read.csv(file.path(Data_wd, "/predict_grid_bibs_plots_cs_500m.csv"))

# Extract covariates in the model and coordinates
grid <- as.data.frame(cbind(x = grid.tmp1$x,
                            y = grid.tmp1$y,
                            imperv.500m = grid.tmp1$impervious_surface,
                            open.green.500m = grid.tmp1$open_green))

head(grid)

# Divide the data into groups of 100 rows so the prediction will use less memory 
# doing each group at a time
grid$split <- dismo::kfold( x = grid, k = floor(nrow(grid)/100)+1 )

# Transform into a list of dataframes based on the group assigned in the previous step
grid.list <- split(x = grid, f = grid$split )
head(grid.list[[1]])


# Make all the step for the prediction of the JSDM into a function
my.predict <- function(grid) {
  #DEFINING THE NEW STUDY DESIGN
  myNew <- nrow(grid)
  dfPiNew <- matrix(NA, myNew, 1)
  dfPiNew[,1] <- sprintf('new_Sample_%.3d', 1:myNew)
  dfPiNew = as.data.frame(dfPiNew)
  colnames(dfPiNew) = m$levelNames
  dfPiAll=rbind(m$dfPi,dfPiNew)
  #DEFINING RANDOM EFFECTS THAT INCLUDE BOTH OLD (USED FOR MODEL FITTING) AND NEW (THE GRID DATA) UNITS
  rL1=m$rL[[1]]$clone()
  xy = grid[,1:2]
  xyold = rL1$s
  rownames(xy) = dfPiNew[,1]
  colnames(xy) = colnames(xyold)
  xyall = rbind(xyold,xy)
  rL1$pi = dfPiAll[,1]
  rL1$s = xyall
  mystart <- proc.time()
  #predYR1=m$predict(dfPiNew = dfPiNew, XData = grid, rL = list(rL1), expected = FALSE, predictEtaMean = TRUE)
  predYR = apply(abind::abind(m$predict(dfPiNew = dfPiNew, XData = grid, rL = list(rL1), expected = FALSE, predictEtaMean = TRUE),along=3),c(1,2),mean)
  myend <- proc.time()
  a <- a+1
  print(a)
  print(myend -mystart)
  return(predYR)
}

fulltime <- proc.time()
a <- 0

# Test if the function works
# grid <- grid[1000:1005,]
# my.predict(grid)


# Apply the function to all elements in the list
res.list <- lapply(grid.list, my.predict)
# Paste all results in the same dataframe
res.df <- do.call("rbind", res.list)
# res.df <- res.df[order()]
# save results
write.csv(res.df , paste0(Output_wd, "/Predict_fullmap_bees2.csv"))

fulltime.end <- proc.time()
fulltime.end - fulltime

#Checking everything is correct
head(res.df)
nrow(res.df)
nrow(grid)

###################################################
###################################################

#y <- c( 116.25, 125.41, 145.24, 154.67, 170.56)
#x <- c( 10, 20, 30, 40, 50)
#mod <- lm(y ~ x)
#new <- data.frame(x = seq(0, 7000, 100))
#predict(lm(y ~ x), new, se.fit = TRUE)
#9752.662/60/60


###################################################
###################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~
## Ploting the map
#~~~~~~~~~~~~~~~~~~~~~~~~~

## Load data

bees.berlin <- read.csv(paste0(Output_wd, "/Predict_fullmap_bees2.csv"))
bees.berlin <- bees.berlin[order(bees.berlin$X),]
write.csv(bees.berlin , paste0(Output_wd, "/Predict_fullmap_bees2.csv"), row.names = F)
bees.berlin <- read.csv(paste0(Output_wd, "/Predict_fullmap_bees2.csv"))

head(bees.berlin)
str(bees.berlin)

# Extract abundance a richness 
bees.Abu <- rowSums(bees.berlin[,-1])# Abundance per site
bees.S <- rowSums(bees.berlin[,-1] > 1)# Richness per site


#####################
# Mapping

# Environment
## berlin and water
# berlin <- read_sf(paste0(Maps_wd, "/berlin_city_border_25833.gpkg"))
berlin.out <- read_sf(paste0(Maps_wd,"/Berlin_border.gpkg"))

berlin.out <- berlin.out %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin.out <- st_transform(berlin.out, crs = 25833)
plot(st_geometry(berlin.out))

water <- read_sf(paste0(Maps_wd, "/waterbodies_Berlin_25833.gpkg"))
plot(st_geometry(water))

# A raster to paste the values:
# Load example raster to have extension
focal.500m <- stack(paste0(Maps_wd, "/FocalMean500m_Stack.tif"))
names.focal <- read.csv(paste0(Maps_wd, "/names_focal_stack.csv"))
names(focal.500m) <- names.focal$x


# create a raster with same extension and resolution as one in the stack
raster1 <- raster(focal.500m$distance_water)
# Transform it to 500m resolution (as the original grid used for prediction)
raster.500m <- aggregate(raster1, fact = 25)
# Paste values of abundance and richness
raster.abu <- setValues(raster.500m, bees.Abu)
raster.S <- setValues(raster.500m, bees.S)


rasters.bees <- stack(raster.abu, raster.S)
names(rasters.bees) <- c("Bees.abund", "Bees.richness")


rasters.bees.mask <- mask(rasters.bees, berlin.out)


library(RColorBrewer)
library(colorRamps)

pdf(paste0(Output_wd, "/Berlin_map_bees_abund_res500m.pdf"))
tm_shape(rasters.bees.mask) +
  tm_raster( "Bees.abund", title = "Wild bees abundance",
             # palette = rev(heat.colors(6)),
             # breaks = c(-Inf, 15, 20, 25, 30, 35, Inf),
             # palette = colorRamps::matlab.like2(4), 
             # palette = "YlOrRd",
             palette = "OrRd",
             breaks = c(0, round(quantile(rasters.bees$Bees.abund)[2:4], 0), round(cellStats(rasters.bees$Bees.abund,max), 0))) +
  tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  # tm_text("spatial_al") 
  tm_shape(water) +
  tm_polygons("lightblue", border.col = NULL) +
  tm_layout(main.title = "Wild bees abundance \nres = 500m",
            compass.type = "4star",
            legend.title.size = 1.4,
            legend.text.size = 0.8) +
  tm_compass(position = c("left", "bottom"),
             size = 2) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_legend(position = c("right", "top"),
            frame = TRUE,
            bg.color="lightyellow")
dev.off()


pdf(paste0(Output_wd, "/Berlin_map_bees_richness_res500m.pdf"))
tm_shape(rasters.bees.mask) +
  tm_raster( "Bees.richness", title = "Wild bees richness",
             # palette = rev(heat.colors(6)),
             # breaks = c(-Inf, 15, 20, 25, 30, 35, Inf),
             # palette = colorRamps::matlab.like2(4), 
             # palette = "YlOrRd",
             palette = "Greens",
             breaks = c(0, round(quantile(rasters.bees$Bees.richness)[2:4], 0), round(cellStats(rasters.bees$Bees.richness,max), 0))) +
  tm_shape(berlin.out) +
  tm_borders("ivory4", lwd = 2.5) +
  # tm_text("spatial_al") 
  tm_shape(water) +
  tm_polygons("lightblue", border.col = NULL) +
  tm_layout(main.title = "Wild bees richness \nres = 500m",
            compass.type = "4star",
            legend.title.size = 1.4,
            legend.text.size = 0.8) +
  tm_compass(position = c("left", "bottom"),
             size = 2) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_legend(position = c("right", "top"),
            frame = TRUE,
            bg.color="lightyellow")
dev.off()



