#############################
## SPATIAL SCALE SELECTION ##
#############################

## All spatial information uses epgs 25833

# We select the spatial scale and the most relevant environmental layers by
# Multivariate GLM models for community response to environment 
# using mvabund package (see Warton et al. 2012)

# We are going to do this for the arthropod data obtained in the arthropod sampling plots
# and for the bird data obtained in the monitoring transects

# First, we extracted the environmental variables for the plots and the transects, at the 
# five spatial scales (resolutions) that we are going to compare: 100m, 500m, 1km, 2km, 5km (script 1).
# Then, we test for correlations among environmental variables for the plots and 
# for the transects at all the resolutions (script 2). 
# Finally, we run multiple multivariate GLM at all the resolutions and compare the 
# results via AICc to select the best spatial scale and the best variables.







#############################
## SPATIAL SCALE SELECTION ##
#############################

#####################
# LIBRARIES
library(dplyr)
library(stringr)
library(mvabund)

#####################
# WORKING SPACE
work_wd <- getwd()
sp_wd <- paste0(work_wd, "/2_Raw_data/Species/")
env_wd <- paste0(work_wd, "/4_Processed_data/EnvValues/")
procsp_wd <- paste0(work_wd, "/4_Processed_data/SpValues/")
output_wd <- paste0(work_wd, "/5_Results/ModelSelection3Sites/")



#####################
### Functions 
# Function for AICc
source(file = paste0(work_wd, "/3_Scripts/Functions_MultivariateModelSelection.R"))

# Model formulas
formulas.42plots.3sites <- readRDS(paste0(work_wd, "/4_Processed_data/formulas_42plots_3sites.rds"))
# Family: "negative binomial"


#####################
# INVERTEBRATE SPECIES ANALYSES

### PREDICTORS
# plant data
selected.plots <- read.csv(paste0(sp_wd, "Invertebrates_shared_plots.csv"))
selected.plots <- selected.plots %>% rename(ID_plot = Plot_ID)

plants.tmp <- read.csv(paste0(sp_wd, "plants_sitexsp.csv"), sep = ";")
plants.tmp <- plants.tmp %>% rename(ID_plot = X)
plants <- merge(selected.plots, plants.tmp, by="ID_plot")
plants$abund <- rowSums(plants[2:length(plants)])
plants$richness <- rowSums(plants[,2:length(plants)] != 0)

# ENVIRONMENTAL DATA
env.100m.cs <- read.csv(paste0(env_wd, "bibsplots_env_centered_100m.csv"), row.names = 1)
env.500m.cs <- read.csv(paste0(env_wd, "bibsplots_env_centered_500m.csv"), row.names = 1)
env.1km.cs <- read.csv(paste0(env_wd, "bibsplots_env_centered_1km.csv"), row.names = 1)
env.2km.cs <- read.csv(paste0(env_wd, "bibsplots_env_centered_2km.csv"), row.names = 1)
env.5km.cs <- read.csv(paste0(env_wd, "bibsplots_env_centered_5km.csv"), row.names = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CARABIDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Species data
carabids <- read.csv(paste0(sp_wd, "carabids_sitexsp_namesphylo_42plots.csv"))
carabids.abund <- as.matrix(carabids[,2:length(carabids)]) #remove sites 
ncol(carabids.abund)
# [1] 104 Species in original data

# Remove single and doubletons
carabids.3sites <- carabids.abund[,colSums(carabids.abund>0)>2] # Removing species present only in 1 or 2 sites
ncol(carabids.3sites)
# [1] 73 Species remaining
rowSums(carabids.3sites)
write.csv(carabids.3sites, paste0(procsp_wd, "carabids_42plots_3sites.csv"))

# Abundance matrix
carabids.mvabund.3sites <-  mvabund(carabids.3sites)
plot(carabids.mvabund.3sites)

myresponse <- carabids.mvabund.3sites

### MODELS 
# Get names of the models
names.forms <- names(formulas.42plots.3sites)

splitnames <- strsplit(names.forms, ".", fixed = TRUE)
model.type <- vector()
for (i in 1:length(splitnames)){
  model.type[i] <- paste0(splitnames[[i]][2], ".", splitnames[[i]][3]) 
}
model.type

# Multivariate models with mvabund for each formula
models.carabids.3sites <- list()
for (i in 1:length(model.type)){
  models.carabids.3sites[[i]] <- manyglm(formulas.42plots.3sites[[i]],
                                         family = "negative.binomial")
  names(models.carabids.3sites)[i] <- model.type[i]
}

models.carabids.3sites[[2]]
saveRDS(models.carabids.3sites, paste0(output_wd, "carabids_mvmodels_3sites.rds"))

# This is to check models outputs. Not running because computational time
# summary(models.carabids.3sites[[2]])
# anova(models.carabids.3sites[[2]])

# AICc values of the models
aic.full <- vector()
aicc.full <- vector()
aicc.avg <- vector()

for (i in 1:length(models.carabids.3sites)){
  aic.full[i] <- models.carabids.3sites[[i]]$AIC # AIC value as sum given by the model
  aicc.full[i] <- get.AICc(mymodel = models.carabids.3sites[[i]], myresponse = carabids.mvabund.3sites) # AICc value of the sum of univariate models
  aicc.avg[i] <- get.avgAICc(mymodel = models.carabids.3sites[[i]], myresponse = carabids.mvabund.3sites) # AICc value averaged by number of univariate models
}

AIC.carabids.3sites <- cbind.data.frame (model.type, aic.full, aicc.full, aicc.avg)
AIC.carabids.3sites

get.avgAICc(mymodel = models.carabids.3sites[[2]], myresponse = carabids.mvabund.3sites)
get.AICc(mymodel = models.carabids.3sites[[2]], myresponse = carabids.mvabund.3sites)/ncol(carabids.mvabund.3sites)
AIC.carabids.3sites$aic.full[2]
models.carabids.3sites[[2]]$AIC

write.csv(AIC.carabids.3sites, paste0(output_wd, "carabids_AIC_3sites.csv"))


## Look at significance of seleted models
# Only shown for one model due to computational time  
anovam4.carabids <- anova(models.carabids.3sites[[4]], p.uni = "adjusted")
anovam4.carabids$table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRASSHOPPERS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
# Species data
grasshop <- read.csv(paste0(sp_wd, "grasshoppers_sitexsp_namesphylo_42plots.csv"))
grasshop.abund <- as.matrix(grasshop[,2:length(grasshop)]) #remove sites 
ncol(grasshop.abund)
# [1] 20 Species in original data

# Remove single and doubletons
grasshop.3sites <- grasshop.abund[,colSums(grasshop.abund>0)>2] # Removing species present only in 1 or 2 sites
ncol(grasshop.3sites)
# [1] 18 Species remaining
write.csv(grasshop.3sites, paste0(procsp_wd, "grasshoppers_42plots_3sites.csv"))
rowSums(grasshop.3sites)

# Abundance matrix
grasshop.mvabund.3sites <-  mvabund(grasshop.3sites)
plot(grasshop.mvabund.3sites)

myresponse <- grasshop.mvabund.3sites

### MODELS 
# Get names of the models
names.forms <- names(formulas.42plots.3sites)

splitnames <- strsplit(names.forms, ".", fixed = TRUE)
model.type <- vector()
for (i in 1:length(splitnames)){
  model.type[i] <- paste0(splitnames[[i]][2], ".", splitnames[[i]][3]) 
}
model.type

# Multivariate models with mvabund for each formula
models.grasshop.3sites <- list()
for (i in 1:length(model.type)){
  models.grasshop.3sites[[i]] <- manyglm(formulas.42plots.3sites[[i]],
                                         family = "negative.binomial")
  names(models.grasshop.3sites)[i] <- model.type[i]
}

models.grasshop.3sites[[2]]
saveRDS(models.grasshop.3sites, paste0(output_wd, "grasshoppers_mvmodels_3sites.rds"))

# This is to check models outputs. Not running because computational time
# summary(models.carabids.3sites[[2]])
# anova(models.carabids.3sites[[2]])

# AICc values of the models
aic.full <- vector()
aicc.full <- vector()
aicc.avg <- vector()

for (i in 1:length(models.grasshop.3sites)){
  aic.full[i] <- models.grasshop.3sites[[i]]$AIC # AIC value as sum given by the model
  aicc.full[i] <- get.AICc(mymodel = models.grasshop.3sites[[i]], myresponse = grasshop.mvabund.3sites) # AICc value of the sum of univariate models
  aicc.avg[i] <- get.avgAICc(mymodel = models.grasshop.3sites[[i]], myresponse = grasshop.mvabund.3sites) # AICc value averaged by number of univariate models
}

AIC.grasshop.3sites <- cbind.data.frame (model.type, aic.full, aicc.full, aicc.avg)
AIC.grasshop.3sites

get.avgAICc(mymodel = models.grasshop.3sites[[2]], myresponse = grasshop.mvabund.3sites)
get.AICc(mymodel = models.grasshop.3sites[[2]], myresponse = grasshop.mvabund.3sites)/ncol(grasshop.mvabund.3sites)
AIC.grasshop.3sites$aic.full[2]
models.grasshop.3sites[[2]]$AIC

write.csv(AIC.grasshop.3sites, paste0(output_wd, "grasshoppers_AIC_3sites.csv"))


## Look at significance of seleted models
# Only shown for one model due to computational time  
anovam4.grasshop <- anova(models.grasshop.3sites[[4]], p.uni = "adjusted")
anovam4.grasshop$table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPIDERS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Species data
spiders <- read.csv(paste0(sp_wd, "spiders_sitexsp_namesphylo_42plots.csv"))
spiders.abund <- as.matrix(spiders[,2:length(spiders)]) #remove sites 
ncol(spiders.abund)
# [1] 182 Species in original data

# Remove single and doubletons
spiders.3sites <- spiders.abund[,colSums(spiders.abund>0)>2] # Removing species present only in 1 or 2 sites
ncol(spiders.3sites)
# [1] 112 Species remaining
write.csv(spiders.3sites, paste0(procsp_wd, "spiders_42plots_3sites.csv"))
rowSums(spiders.3sites)

# Abundance matrix
spiders.mvabund.3sites <-  mvabund(spiders.3sites)
plot(spiders.mvabund.3sites)

myresponse <- spiders.mvabund.3sites

### MODELS 
# Get names of the models
names.forms <- names(formulas.42plots.3sites)

splitnames <- strsplit(names.forms, ".", fixed = TRUE)
model.type <- vector()
for (i in 1:length(splitnames)){
  model.type[i] <- paste0(splitnames[[i]][2], ".", splitnames[[i]][3]) 
}
model.type

# Multivariate models with mvabund for each formula
models.spiders.3sites <- list()
for (i in 1:length(model.type)){
  models.spiders.3sites[[i]] <- manyglm(formulas.42plots.3sites[[i]],
                                         family = "negative.binomial")
  names(models.spiders.3sites)[i] <- model.type[i]
}

models.spiders.3sites[[2]]
saveRDS(models.spiders.3sites, paste0(output_wd, "spiders_mvmodels_3sites.rds"))

# AICc values of the models
aic.full <- vector()
aicc.full <- vector()
aicc.avg <- vector()

for (i in 1:length(models.spiders.3sites)){
  aic.full[i] <- models.spiders.3sites[[i]]$AIC # AIC value as sum given by the model
  aicc.full[i] <- get.AICc(mymodel = models.spiders.3sites[[i]], myresponse = spiders.mvabund.3sites) # AICc value of the sum of univariate models
  aicc.avg[i] <- get.avgAICc(mymodel = models.spiders.3sites[[i]], myresponse = spiders.mvabund.3sites) # AICc value averaged by number of univariate models
}

AIC.spiders.3sites <- cbind.data.frame (model.type, aic.full, aicc.full, aicc.avg)
AIC.spiders.3sites

get.avgAICc(mymodel = models.spiders.3sites[[2]], myresponse = spiders.mvabund.3sites)
get.AICc(mymodel = models.spiders.3sites[[2]], myresponse = spiders.mvabund.3sites)/ncol(spiders.mvabund.3sites)
AIC.spiders.3sites$aic.full[2]
models.spiders.3sites[[2]]$AIC

write.csv(AIC.spiders.3sites, paste0(output_wd, "spiders_AIC_3sites.csv"))


## Look at significance of seleted models
# Only shown for one model due to computational time  
anovam4.spiders <- anova(models.spiders.3sites[[4]], p.uni = "adjusted")
anovam4.spiders$table




#####################
# BIRD SPECIES ANALYSES

### PREDICTORS
trans.100m.cs <- read.csv(paste0(env_wd, "transects_env_centered_100m.csv"), row.names = 1)
trans.500m.cs <- read.csv(paste0(env_wd, "transects_env_centered_500m.csv"), row.names = 1)
trans.1km.cs <- read.csv(paste0(env_wd, "transects_env_centered_1km.csv"), row.names = 1)
trans.2km.cs <- read.csv(paste0(env_wd, "transects_env_centered_2km.csv"), row.names = 1)
trans.5km.cs <- read.csv(paste0(env_wd, "transects_env_centered_5km.csv"), row.names = 1)

#####################
# SPECIES DATA

# Species data
birds.tmp <- read.csv(paste0(sp_wd, "Birds_abund_noexo_norare2_noaquatic_bykm.csv"))
birds.tmp$site <- tolower(birds.tmp$site) 

## Put environment and birds in the same order
cbind(rownames(trans.100m.cs), birds.tmp$site)

birds <- birds.tmp[match(rownames(trans.100m.cs), birds.tmp$site),]

birds.abund <- as.matrix(birds[,1:ncol(birds)-1]) #remove sites 

# Abundance matrix
birds.mvabund <-  mvabund(birds.abund)
plot(birds.mvabund)

myresponse <- birds.mvabund


## MODELS
# Model formulas
formulas.birds <- readRDS(paste0(work_wd, "/4_Processed_data/formulas_birs.rds"))
# Family: "negative binomial"

# Get names of the models
names.forms <- names(formulas.birds)

splitnames <- strsplit(names.forms, ".", fixed = TRUE)
model.type <- vector()
for (i in 1:length(splitnames)){
  model.type[i] <- paste0(splitnames[[i]][3], ".", splitnames[[i]][4]) 
}
model.type

## Multivariate models with mvabund for each formula
models.birds <- list()
for (i in 1:length(model.type)){
  models.birds[[i]] <- manylm(formulas.birds[[i]])
  names(models.birds)[i] <- model.type[i]
}

models.birds[[5]]
names(models.birds)[5]
saveRDS(models.birds, paste0(output_wd, "birds_mvmodels_3sites.rds"))


## AICc value
aic.full <- vector()
aicc.full <- vector()
aicc.avg <- vector()

for (i in 1:length(models.birds)){
  aic.full[i] <- sum(AIC(models.birds[[i]])) # total AIC value
  aicc.full[i] <- get.AICc.gaus(mymodel = models.birds[[i]], myresponse = birds.mvabund) # AICc value of the sum of univariate models
  aicc.avg[i] <- get.avgAICc.gaus(mymodel = models.birds[[i]], myresponse = birds.mvabund) # AICc value averaged by number of univariate models
}

AIC.birds <- cbind.data.frame (model.type, aic.full, aicc.full, aicc.avg)
AIC.birds

get.avgAICc(mymodel = models.birds[[2]], myresponse = birds.mvabund)
get.AICc(mymodel = models.birds[[2]], myresponse = birds.mvabund)/ncol(birds.mvabund)
AIC.birds$aic.full[2]
sum(AIC(models.birds[[2]]))

write.csv(AIC.birds, paste0(output_wd, "birds_AIC_3sites.csv"))


## Look at significance of seleted models
anovam4.birds <- anova(models.birds[[4]])
anovam4.birds
anovam4.birds$table


