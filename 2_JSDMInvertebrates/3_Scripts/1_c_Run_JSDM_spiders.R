
### Joint Species Distribution Models for Invertebrate Taxa in Berlin

## Script for Spiders

### Packages
library(Hmsc)


### Workspace
WorkDir <- getwd()

Env_wd <- file.path(WorkDir, "2_Raw_data/EnvValues/")
Species_wd <- file.path(WorkDir, "2_Raw_data/Species")
Output_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_Spiders_poisson/")


## Species data

# spiders data
spiders <- read.csv(paste0(Species_wd, "/spiders_42plots_3sites.csv"))
str(spiders)
head(spiders)
ncol(spiders)
# [1] 112 Species remaining
Y.spiders <- as.matrix(spiders)


## Environmental data
env100m.tmp <- read.csv(paste0(Env_wd, "/bibsplots_allenvir_100m.csv"))
str(env100m.tmp)
head(env100m.tmp)

Covariates <- cbind.data.frame(imperv.100m = env100m.tmp$impervious_surface,
                               noise.100m = env100m.tmp$noise,
                               open.green.100m = env100m.tmp$open_green,
                               temp.100m = env100m.tmp$temp_day)
summary(Covariates)


## spatial coordinates
xyData.tmp <- read.csv(paste0(Env_wd, '/env_42plots.csv'), sep = ';')
head(xyData.tmp)
xyData <- cbind(x = xyData.tmp$X, y = xyData.tmp$Y)


## Set up the model

### Regression formulas for environmental covariates
XFormula.spiders = ~ imperv.100m + noise.100m + temp.100m + open.green.100m


### Study design
studyDesign <- data.frame(site = xyData.tmp$Plot)
sRL <- xyData
rownames(sRL) <- xyData.tmp$Plot
rL <- HmscRandomLevel(sData = sRL)

## Fit models for abundance  data 

### Define MCMC parameters
thin <- 100
samples <- 10000
transient <- 10000
nChains <- 3
verbose <- 10000


### Poisson distribution ### 

# Set the components of each model
m <- Hmsc(Y = Y.spiders, XData = Covariates, XFormula = XFormula.spiders,
          studyDesign=studyDesign, ranLevels=list(site=rL),
          distr = "poisson")
# run Markov Chains
m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose,
                nParallel = 3)
#Save the model
filename = paste0(Output_wd, "/model_pois_spiders_abu_3sites.rds")
saveRDS(m,file=filename)



#########################################
# spiders model: model check
#########################################

### Packages
library(coda)
library(corrplot)
library(MCMCvis)
library(Hmsc)


#######################
## Model convergence 

# Evaluate convergence: Effective sample size and gelman-rubin diagnostic (potencial reduction factor)
mpost <- convertToCodaObject(m)
# Numerical output
ess.beta <- effectiveSize(mpost$Beta)
gd.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
ess.gamma <- effectiveSize(mpost$Gamma)
gd.gamma <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf 
ess.omega <- effectiveSize(mpost$Omega[[1]])
gd.omega <- gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf
convergence.names <- c("ess.beta", "ess.gamma", "ess.omega",
                       "gd.beta", "gd.gamma", "gd.omega")
convergence.list <- list(ess.beta, ess.gamma, ess.omega,
                         gd.beta, gd.gamma, gd.omega)
for (i in 1:length(convergence.names)){
  write.csv(convergence.list[[i]], paste0(Output_wd, "/", convergence.names[i], ".csv")) 
}

# Graphical output
png(paste0(Output_wd, "/spiders_model_convergence_hist.png"), width = 800, height = 1000,
    pointsize = 20)
par(mfrow=c(3,2))
hist(ess.beta, main = "ess(beta)_spiders")
hist(ess.gamma, main = "ess(gamma)_spiders")
hist(ess.omega, main = "ess(omega)_spiders")
hist(gd.beta, main = "psrf(beta)_spiders")
hist(gd.gamma, main = "psrf(gamma)_spiders")
hist(gd.omega, main = "psrf(omega)_spiders")
dev.off()

# Save plots of the chains
MCMCtrace(mpost$Beta, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "spiders_MCMCtrace_beta",
          wd = Output_wd)
MCMCtrace(mpost$Gamma, 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "spiders_MCMCtrace_gamma",
          wd = Output_wd)
MCMCtrace(mpost$Omega[[1]], 
          pdf = TRUE, 
          open_pdf = FALSE,
          filename = "spiders_MCMCtrace_omega",
          wd = Output_wd)

par(mfrow=c(1,1))



#######################
## Model fot 

# Explanatory R2. Get predictions for the observed values (in poisson scale)
preds <- computePredictedValues(m, expected = FALSE)
preds.values <- apply(abind(preds,along=3),c(1,2),median)

# R2 with the built in function
modelr2.explanatory <- evaluateModelFit(hM = m, predY = preds)
modelr2.explanatory

# R2 Mannually for the species
R2.sp <- matrix(NA, m$ns, 1)
for (i in 1:m$ns) {
  R2.sp[i, ] <- cor(preds.values[, i],m$Y[, i])^2
}

mean(R2.sp, na.rm=TRUE)
# [1] 0.722537

# R2 Mannually for the sites
R2.site <- matrix(NA, m$ny, 1)
for (i in 1:m$ny) {
  R2.site[i, ] <- cor(preds.values[i, ], m$Y[i, ])^2
}
mean(R2.site)
# [1] 0.9460594


############################################
## Estimates of Beta and Gamma parameters

# We are going to obtain the values numerically and visually.

## load function to plot the values
source(paste0(WorkDir, "/3_Scripts/Function_PlotBetas.R"))


# Beta values
Beta.results <- as.data.frame(MCMCsummary(mpost$Beta))
write.csv(Beta.results, paste0(Output_wd, "/spiders_beta_coeffients.csv"), 
          row.names = FALSE)
# Default beta plot in Hmsc package
postBeta <- getPostEstimate(m, parName = "Beta")

# Print a plot for each predictor
plot.cov.betas(m, modelType = "spiders") 

# Gamma values
Gamma.results <- as.data.frame(MCMCsummary(mpost$Gamma))
write.csv(Gamma.results, paste0(Output_wd, "/spiders_gamma_coeffients.csv"), 
          row.names = FALSE)
# Default gamma plot in Hmsc package
postGamma <- getPostEstimate(m, parName = "Gamma")

# Coef. gammma
MCMCplot(mpost$Gamma, ref_ovl = TRUE)

############################################
## Speceis co-occurrences 

OmegaCor <- computeAssociations(m)
supportLevel <- 0.95

# Default plot in Hmsc package
toPlot <- ((OmegaCor[[1]]$support > supportLevel)
           + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
png(paste0(Output_wd, "/spiders_default_omegaplot.png"))
corrplot(toPlot, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = paste0("random effect level: ", m$rLNames[1]), 
         tl.cex = 0.4, 
         mar = c(0,0,1,0))
dev.off()


############################################
## Variance partitioning

# The order of the variables, if they are continuous, is 
# 1. intercept(this can be in any group)
# 2. first variable
# 3. second variable
# ETC.

# The formula we used for running the model is: 
# XFormula.spiders 

head(m$X)

VP <- computeVariancePartitioning(m, group = c(1,1,1,2,2), groupnames = c("urban", "nature"))

# plot var part
png(paste0(Output_wd, "/spiders_default_VP.png"), 
    width = 800)
plotVariancePartitioning(m, VP = VP, las = 2, cex.names = 0.4)
title(main = "\n \nspiders")
dev.off()



#########################################
# Spiders model: predictions
#########################################

### Packages
library(abind)
library(ggplot2)
library(sf)
library(tmap)
library(Hmsc)


### Working directories
WorkDir <- getwd()

Pred_wd <- file.path(WorkDir, "2_Raw_data/ToPredict")
Output_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_Spiders_poisson/Predictions")


m$XFormula
# ~ imperv.100m + noise.100m + temp.100m + open.green.100m


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Predicting in Bird Transects #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### READING THE GRID DATA FOR BIRD TRANSECTS
grid.tmp100 <- read.csv(file.path(Pred_wd, "/transects_allenvir_100m.csv"))

grid <- as.data.frame(cbind(imperv.100m = grid.tmp100$impervious_surface,
                            noise.100m = grid.tmp100$noise,
                            temp.100m = grid.tmp100$temp_day,
                            open.green.100m = grid.tmp100$open_green))
xy.transects <- read.csv(file.path(Pred_wd, "/transects_coordinates_25833.csv"), sep=",")

### Defining the new study design
StudyDesignNew <- as.data.frame(xy.transects$routcode_2)
colnames(StudyDesignNew) <- "site"
StudyDesignAll <- rbind(m$studyDesign,StudyDesignNew)

### Defining random effects that include both old (model fitting) and 
# new (grid data) units  
rL1 <- m$rL[[1]]
xyold <- rL1$s
xy <- xy.transects[,2:3]
rownames(xy) <- StudyDesignNew[,1]
colnames(xy) <- colnames(xyold)
xyall <- rbind(xyold,xy)
rL1$pi <- StudyDesignAll[,1]
rL1$s <- xyall

# Predictions in poisson scale
predYR2 <- predict(m, studyDesign=StudyDesignNew, XData=grid, ranLevels=list(site=rL1),
                   expected = FALSE, predictEtaMean=TRUE)
# We use the median instead of the mean to reflect better the results and respect the poisson distribution
predY2 <- apply(abind(predYR2,along=3),c(1,2), median)
rownames(predY2) <- xy.transects$routcode_2
m$Y
write.csv(predY2, paste0(Output_wd, "/spiders_predictions_pois_birdtransects.csv"))


# COMPUTE ABUNDANCE (Abu), SPECIES RICHNESS (S), COMMUNITY WEIGHTED MEANS (predT),
# REGIONS OF COMMON PROFILE (RCP)

Abu <- rowSums(predY2)
S <- rowSums(predY2 > 0) # species presence with more than 0.5 abundance
predT <- (predY2%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
RCP <- kmeans(predY2, 4) # four regions
RCP$cluster <- as.factor(RCP$cluster)
# EXTRACT THE OCCURRENCE PROBABILITIES OF ONE EXAMPLE SPECIES
pred_sp10 <- predY2[,10]
# MAKE A DATAFRAME OF THE DATA TO BE PLOTTED
mapData <- data.frame(xy,Abu, S,predT,pred_sp10,RCP$cluster)

# We are now ready to map the predicted values onto the coordinates 
sp10 <- ggplot(data = mapData, aes(x=x, y=y, color=pred_sp10))+geom_point(size=1)
sp10 + ggtitle("Predicted Species 10 occurrence") +
  xlab("East coordinate") + ylab("North coordinate") +
  scale_color_gradient(low="blue", high="red", name ="Abundance") # with p/a this is Occurrence probability

sp.RCP <- ggplot(data = mapData, aes(x=x, y=y, color=RCP$cluster))+geom_point(size=1)
sp.RCP + ggtitle("Common profile zones") +
  xlab("East coordinate") + ylab("North coordinate") +
  scale_color_gradient(low="blue", high="red", name ="RCP group") # with p/a this is Occurrence probability

sp.abu <- ggplot(data = mapData, aes(x=x, y=y, color=Abu))+geom_point(size=1)
sp.abu + ggtitle("Predicted species abundance") +
  xlab("East coordinate") + ylab("North coordinate") +
  scale_color_gradient(low="blue", high="red", name ="Abundance")

sp.richness <- ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=1)
sp.richness + ggtitle("Predicted species richness") +
  xlab("East coordinate") + ylab("North coordinate") +
  scale_color_gradient(low="blue", high="red", name ="Species richness")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Predicting **carabid** community in Berlin #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# We are going to use **100 m resolution**
# The goal of the prediction is to obtain a value that will be used to predict bird abundance in the city. 

# For the process we follow these steps:
# 1. Load the JSDM model that we'll use for prediction: *model_pois_carabid_abu_3sites.rds*
# 2. A dataframe that contains the environmental variables of a 100m resolution grid over the city
# with spatial coordinates
# 3. Divide the new environmental data into subsets to **parallelize** the predictions, or do it in "small" groups.
# 4. Run the predictions and extract total abundance and species richness for each grid cell. 


### Packages
library(Hmsc)
library(dplyr) # rearrange data
library(dismo) # subset for parallelizing


### Workspace
WorkDir <- getwd()
Env_wd <- file.path(WorkDir, "2_Raw_data/EnvValues/")
Preds_wd <- file.path(WorkDir, "2_Raw_data/ToPredict/")
Output_wd <- file.path(WorkDir, "5_Results/InvertebrateModels_Spiders_poisson/Predictions/")


### Load environmental grid
grid.topredict.100mres <- read.csv(paste0(Preds_wd, "/grid_topredict_100mres.csv"))
head(grid.topredict.100mres)

# The variables we are going need to have the same name as in the model

m$XFormula # ~ imperv.100m + dist.water.100m + temp.100m
colnames(grid.topredict.100mres)
grid.100mres.spiders <- grid.topredict.100mres %>%
  dplyr::select(x = x_100m, y = y_100m, 
                imperv.100m = impervious_surface_100m,
                noise.100m = noise_100m,
                temp.100m = temp_day_100m,
                open.green.100m = open_green_100m)
head(grid.100mres.spiders)


# Remove NAs to improve code efficiency
grid2.100mres <- grid.100mres.spiders[complete.cases(grid.100mres.spiders), ]
head(grid2.100mres)
str(grid2.100mres)


### Subset for parallelize
nrow(grid2.100mres)
# [1] 90574

# Assign a group to each observation. k is the number of groups. +1 is used to include all the observation
# not included in the exact division
grid2.100mres$split <- dismo::kfold( x = grid2.100mres, k = floor(nrow(grid2.100mres)/500)+1 )

# Transform into a list of dataframes based on the group assigned in the previous step
grid2.100mres.list <- split(x = grid2.100mres, f = grid2.100mres$split )

print(paste0("Number of folds = ", length(grid2.100mres.list), " (Groups of 500 rows)"))
length(grid2.100mres.list)
head(grid2.100mres.list[[1]])


### Function to predict 
# To simplify our procces we make all the step for the prediction of the JSDM into a function

my.predict.fun <- function(grid) {
  #DEFINING THE NEW STUDY DESIGN
  myNew <- nrow(grid)
  StudyDesignNew <- matrix(NA, myNew, 1)
  StudyDesignNew[,1] <- sprintf('new_Sample_%.3d', 1:myNew)
  StudyDesignNew <- as.data.frame(StudyDesignNew)
  colnames(StudyDesignNew) = colnames(m$studyDesign)
  StudyDesignAll <- rbind(m$studyDesign,StudyDesignNew)
  #DEFINING RANDOM EFFECTS THAT INCLUDE BOTH OLD (USED FOR MODEL FITTING) AND NEW (THE GRID DATA) UNITS
  rL1 <- m$rL[[1]]
  xyold <- rL1$s
  x <- grid$x
  y <- grid$y
  xy <- cbind(x, y) 
  rownames(xy) <- StudyDesignNew[,1]
  colnames(xy) <- colnames(xyold)
  xyall <- rbind(xyold,xy)
  rL1$pi <- StudyDesignAll[,1]
  rL1$s <- xyall
  # Predict poisson (integers)
  predYR <- predict(m, studyDesign=StudyDesignNew, XData=grid, ranLevels=list(site=rL1),
                   expected = FALSE, predictEtaMean=TRUE)
  # We use the median instead of the mean to reflect better the results and respect the poisson   distribution in the final results
  predY <- apply(abind(predYR,along=3),c(1,2), median)
  return(predY)
}


## Run predictions With parallelization
#Make the code use several cores at the same time
cl <- makeCluster(detectCores()-2) # Using all but two cores 
clusterExport(cl, c("m")) # Copy objects to all cores
clusterEvalQ(cl, library("Hmsc")) # copy libraries to all cores
res.list <- parLapply(cl, grid2.100mres.list, my.predict.fun) # Apply parallel functions
stopCluster(cl) # stop parallelization

# Paste all results in the same dataframe
res.df <- do.call("rbind", res.list)

# Reorder back to original order
my.order <- as.numeric(row.names(res.df))
res.df <- cbind.data.frame(my.order, res.df)
res.df <- res.df[order(res.df$my.order),]
head(res.df)
tail(res.df)

## Summary data
# Abundance per site. We remove the first columns because it indicates the number of the line
spiders.pred.abu <- rowSums(res.df[,-1], na.rm=TRUE)
# Richness per site. We decide to use an abundance higher than 0.5 as presence
spiders.pred.rich <- rowSums(res.df[,-1] > 0.5, na.rm=TRUE)
# Add coordinates and put all the data together
spiders.pred <- cbind(x = grid2.100mres$x, y = grid2.100mres$y, res.df, 
                      spiders.pred.abu = spiders.pred.abu, 
                      spiders.pred.rich = spiders.pred.rich)
# Save data
write.csv(spiders.pred , paste0(Output_wd, "/Predict_Berlin_100m_spiders_3sites.csv"), 
          row.names = F)
