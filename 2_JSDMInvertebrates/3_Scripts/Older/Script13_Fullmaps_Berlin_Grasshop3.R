## 5 November 2018


## Prediction map for grasshoppers3 
## Used a grid of 500m

# Model: grasshop3
# without species only in 1 or 2 sites #

# Env predictors: ~imperv.100m + noise.100m + temp.100m + open.green.100m


rm(list = ls())
setMKLthreads(1)
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

filename = file.path(Model_wd, "model_pois_grasshop3_abu_3sites.Rdata")
load(file = filename)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spatial prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# See the formula of the model
m$XFormula
# ~imperv.100m + noise.100m + temp.100m + open.green.100m


#READING THE GRID DATA (Covariate values extracted from a raster of 20m resolution, Berlin)
# Variable values are standardized!
grid.tmp1 <- read.csv(file.path(Data_wd, "/predict_grid_bibs_plots_cs_100m.csv"))

# Extract covariates in the model and coordinates
grid <- as.data.frame(cbind(x = grid.tmp1$x,
                            y = grid.tmp1$y,
                            imperv.100m = grid.tmp1$impervious_surface,
                            noise.100m = grid.tmp1$noise,
                            temp.100m = grid.tmp1$temp_day,
                            open.green.100m = grid.tmp1$open_green))

head(grid)

# Remove NAs
grid2 <- grid[complete.cases(grid), ]
head(grid2)
str(grid2)

# Divide the data into groups of 50 rows so the prediction will use less memory 
# doing each group at a time
grid2$split <- dismo::kfold( x = grid2, k = floor(nrow(grid2)/50)+1 )

# Transform into a list of dataframes based on the group assigned in the previous step
grid2.list <- split(x = grid2, f = grid2$split )
head(grid2.list[[1]])

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
res.list <- lapply(grid2.list, my.predict)
# Paste all results in the same dataframe
res.df <- do.call("rbind", res.list)
my.order <- as.numeric(row.names(res.df))
res.df <- cbind.data.frame(my.order, res.df)
head(res.df)
res.df <- res.df[order(res.df$my.order),]
# save results
write.csv(res.df , paste0(Output_wd, "/Predict_fullmap_grasshop3.csv"), row.names = F)

fulltime.end <- proc.time()
fulltime.end - fulltime
