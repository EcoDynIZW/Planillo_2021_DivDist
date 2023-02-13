####################################################
### TEST FOR CORRELATIONS AMONG ENVIRONMENTAL VALUES 
####################################################

### Packages
library(raster)
library(sf)
library(dplyr)
library(corrplot)
extract = raster::extract # set the function extract of raster package as priority

### WORK SPACE
work_wd <- getwd()
data_wd <- paste0(work_wd, "/2_Raw_data/Environment/")
sp_wd <- paste0(work_wd, "/2_Raw_data/Species/")
fmeans_wd <- paste0(work_wd, "/4_Processed_data/")
env_wd <- paste0(work_wd, "/4_Processed_data/EnvValues/")
output_wd <- paste0(work_wd, "/5_Results/")


#############################################################
## TESTING FOR CORRELATIONS IN ENVIRONMENTAL VALUES: BIBS PLOTS 

## load plot values
plots.100m.all <- read.csv(paste0(env_wd, "bibsplots_allenvir_100m.csv"), row.names = 1)
plots.500m.all <- read.csv(paste0(env_wd, "bibsplots_allenvir_500m.csv"), row.names = 1)
plots.1km.all <- read.csv(paste0(env_wd, "bibsplots_allenvir_1km.csv"), row.names = 1)
plots.2km.all <- read.csv(paste0(env_wd, "bibsplots_allenvir_2km.csv"), row.names = 1)
plots.5km.all <- read.csv(paste0(env_wd, "bibsplots_allenvir_5km.csv"), row.names = 1)


### DATA: 100 x 100m
# Center and Scale 
plots.100m.all

centered.var.100m <- plots.100m.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.100m, paste0(env_wd, "bibsplots_env_centered_100m.csv"))

# Correlations 
M <- cor(centered.var.100m)

pdf(file = paste0(output_wd,"correlations/correlations_100m.pdf"), width = 5,height = 5) # this line is to save the plot
corrplot(M, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(M, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()

### DATA: 500 x 500m
# Center and Scale 
centered.var.500m <- plots.500m.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.500m, paste0(env_wd, "bibsplots_env_centered_500m.csv"))

# Correlations
M <- cor(centered.var.500m) 

pdf(file = paste0(output_wd,"correlations/correlations_500m.pdf"), width = 5,height = 5)
corrplot(M, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(M, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()

### DATA: 1 x 1km
# Center and Scale 
centered.var.1km <- plots.1km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.1km, paste0(env_wd, "bibsplots_env_centered_1km.csv"))

# Correlations
M <- cor(centered.var.1km) 

pdf(file = paste0(output_wd,"correlations/correlations_1km.pdf"), width = 5,height = 5)
corrplot(M, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(M, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()

### DATA: 2 x 2km
# Center and Scale 
centered.var.2km <- plots.2km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.2km, paste0(env_wd, "bibsplots_env_centered_2km.csv"))

# Correlations
M <- cor(centered.var.2km) 

pdf(file = paste0(output_wd,"correlations/correlations_2km.pdf"), width = 5,height = 5)
corrplot(M, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(M, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()


### DATA: 5 x 5km
# Center and Scale 
centered.var.5km <- plots.5km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.5km, paste0(env_wd, "bibsplots_env_centered_5km.csv"))

# Correlations
M <- cor(centered.var.5km) 

pdf(file = paste0(output_wd,"correlations/correlations_5km.pdf"), width = 5,height = 5)
corrplot(M, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(M, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()



#############################################################
## TESTING FOR CORRELATIONS IN ENVIRONMENTAL VALUES: BIRD TRANSECTS

transects.100m.all <- read.csv(paste0(env_wd, "transects_allenvir_100m.csv"), row.names = 1)
transects.500m.all <- read.csv(paste0(env_wd, "transects_allenvir_500m.csv"), row.names = 1)
transects.1km.all <- read.csv(paste0(env_wd, "transects_allenvir_1km.csv"), row.names = 1)
transects.2km.all <- read.csv(paste0(env_wd, "transects_allenvir_2km.csv"), row.names = 1)
transects.5km.all <- read.csv(paste0(env_wd, "transects_allenvir_5km.csv"), row.names = 1)

### DATA: 100 x 100m
# Center and Scale 
centered.var.100m.trans <- transects.100m.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.100m.trans, paste0(env_wd, "transects_env_centered_100m.csv"))

# Correlations 
Mt <- cor(centered.var.100m.trans)

pdf(file = paste0(output_wd,"correlations/correlations_birds_100m.pdf"), width = 5,height = 5)
corrplot(Mt, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(Mt, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()


### DATA: 500 x 500m
# Center and Scale 
centered.var.500m.trans <- transects.500m.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.500m.trans, paste0(env_wd, "transects_env_centered_500m.csv"))

# Correlations 
Mt <- cor(centered.var.500m.trans)

pdf(file = paste0(output_wd,"correlations/correlations_birds_500m.pdf"), width = 5,height = 5)
corrplot(Mt, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(Mt, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()


### DATA: 1 x 1km
# Center and Scale 
centered.var.1km.trans <- transects.1km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.1km.trans, paste0(env_wd, "transects_env_centered_1km.csv"))

# Correlations 
Mt <- cor(centered.var.1km.trans)

pdf(file = paste0(output_wd,"correlations/correlations_birds_1km.pdf"), width = 5,height = 5)
corrplot(Mt, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(Mt, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()

### DATA: 2 x 2km
# Center and Scale 
centered.var.2km.trans <- transects.2km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.2km.trans, paste0(env_wd, "transects_env_centered_2km.csv"))

# Correlations 
Mt <- cor(centered.var.2km.trans)

pdf(file = paste0(output_wd,"correlations/correlations_birds_2km.pdf"), width = 5,height = 5)
corrplot(Mt, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(Mt, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()


### DATA: 5 x 5km
# Center and Scale 
centered.var.5km.trans <- transects.5km.all %>% 
        transmute_all(.funs = funs(scale(.))) 
write.csv(centered.var.5km.trans, paste0(env_wd, "transects_env_centered_5km.csv"))

# Correlations 
Mt <- cor(centered.var.5km.trans)

pdf(file = paste0(output_wd,"correlations/correlations_birds_5km.pdf"), width = 5,height = 5)
corrplot(Mt, method = "circle", type = "upper", 
         cl.cex = 0.6, tl.col = "black", tl.pos = "lt",
         tl.cex = 0.7)
corrplot(Mt, method = "number", type = "lower", 
         cl.cex = 0.8, diag = FALSE, tl.pos = "n", cl.pos = "n",
         add = TRUE,
         number.cex = 0.6, number.digits = 2)
dev.off()





#############################################################
## TESTING FOR CORRELATIONS ACROSS SCALES: BIBS PLOTS 

### Packages
library(corrplot)


### Workspace
work_wd <- getwd()
sp_wd <- paste0(work_wd, "/2_Raw_data/Species/")
env_wd <- paste0(work_wd, "/4_Processed_data/EnvValues/")
output_wd <- paste0(work_wd, "/5_Results/correlations/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPERVIOUS SURFACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
imperv <- cbind.data.frame(imperv.100m = plots.100m.all$impervious_surface,
                           imperv.500m = plots.500m.all$impervious_surface,
                           imperv.1km = plots.1km.all$impervious_surface,
                           imperv.2km = plots.2km.all$impervious_surface,
                           imperv.5km = plots.5km.all$impervious_surface)

apply(imperv, 2, FUN = range)
apply(imperv, 2, FUN = mean)

cor(imperv)
#             imperv.100m imperv.500m imperv.1km imperv.2km imperv.5km
# imperv.100m   1.0000000   0.7669408  0.6064266  0.5700458  0.4691011
# imperv.500m   0.7669408   1.0000000  0.9528966  0.8993105  0.7300875
# imperv.1km    0.6064266   0.9528966  1.0000000  0.9705439  0.8182940
# imperv.2km    0.5700458   0.8993105  0.9705439  1.0000000  0.8865559
# imperv.5km    0.4691011   0.7300875  0.8182940  0.8865559  1.0000000

# corrplot
pdf(paste0(output_wd, "corrplot_imperv_scales.pdf"))
corrplot(cor(imperv),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of impervious surface across scales",
         mar = c(4,1,3,1))
corrplot(cor(imperv),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_imperv_scales.pdf"))
par(mfrow = c(3,2))
hist(imperv$imperv.100m, xlim = c(0,100), ylim = c(0, nrow(imperv)), 
     main = "Impervious surfave. Histogram of values in 100m scale")
hist(imperv$imperv.500m, xlim = c(0,100), ylim = c(0, nrow(imperv)), 
     main = "Impervious surfave. Histogram of values in 500m scale")
hist(imperv$imperv.1km, xlim = c(0,100), ylim = c(0, nrow(imperv)), 
     main = "Impervious surfave. Histogram of values in 1km scale")
hist(imperv$imperv.2km, xlim = c(0,100), ylim = c(0, nrow(imperv)), 
     main = "Impervious surfave. Histogram of values in 2km scale")
hist(imperv$imperv.5km, xlim = c(0,100), ylim = c(0, nrow(imperv)), 
     main = "Impervious surfave. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# POPULATION 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop <- cbind.data.frame(pop.100m = plots.100m.all$pop_100m,
                        pop.500m = plots.500m.all$pop_500m,
                           pop.1km = plots.1km.all$pop_1km,
                           pop.2km = plots.2km.all$pop_2km,
                           pop.5km = plots.5km.all$pop_5km)

cor(pop)
# corrplot
pdf(paste0(output_wd, "corrplot_population_scales.pdf"))
corrplot(cor(pop),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of population across scales",
         mar = c(4,1,3,1))
corrplot(cor(pop),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_population_scales.pdf"))
par(mfrow = c(3,2))
hist(pop$pop.100m, xlim = c(0,max(pop[,1])), ylim = c(0, nrow(pop)), 
     main = "Population. Histogram of values in 100m scale")
hist(pop$pop.500m, xlim = c(0,max(pop[,2])), ylim = c(0, nrow(pop)), 
     main = "Population. Histogram of values in 500m scale")
hist(pop$pop.1km, xlim = c(0,max(pop[,3])), ylim = c(0, nrow(pop)), 
     main = "Population. Histogram of values in 1km scale")
hist(pop$pop.2km, xlim = c(0,max(pop[,4])), ylim = c(0, nrow(pop)), 
     main = "Population. Histogram of values in 2km scale")
hist(pop$pop.5km, xlim = c(0,max(pop[,5])), ylim = c(0, nrow(pop)), 
     main = "Population. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DISTANCE TO WATER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dist.wtr <- cbind.data.frame(dist.100m = plots.100m.all$distance_water,
                        dist.500m = plots.500m.all$distance_water,
                        dist.1km = plots.1km.all$distance_water,
                        dist.2km = plots.2km.all$distance_water,
                        dist.5km = plots.5km.all$distance_water)
cor(dist.wtr)

# corrplot
pdf(paste0(output_wd, "corrplot_distwtr_scales.pdf"))
corrplot(cor(dist.wtr),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of distance to water across scales",
         mar = c(4,1,3,1))
corrplot(cor(dist.wtr),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_distwtr_scales.pdf"))
par(mfrow = c(3,2))
hist(dist.wtr$dist.100m, xlim = c(0, 2000), ylim = c(0, nrow(dist.wtr)), 
     main = "Distance to water. Histogram of values in 100m scale")
hist(dist.wtr$dist.500m, xlim = c(0,2000), ylim = c(0, nrow(dist.wtr)), 
     main = "Distance to water. Histogram of values in 500m scale")
hist(dist.wtr$dist.1km, xlim = c(0,2000), ylim = c(0, nrow(dist.wtr)), 
     main = "Distance to water. Histogram of values in 1km scale")
hist(dist.wtr$dist.2km, xlim = c(0,2000), ylim = c(0, nrow(dist.wtr)), 
     main = "Distance to water. Histogram of values in 2km scale")
hist(dist.wtr$dist.5km, xlim = c(0,2000), ylim = c(0, nrow(dist.wtr)), 
     main = "Distance to water. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LIGHT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
light <- cbind.data.frame(light.100m = plots.100m.all$light,
                          light.500m = plots.500m.all$light,
                          light.1km = plots.1km.all$light,
                          light.2km = plots.2km.all$light,
                          light.5km = plots.5km.all$light)

cor(light)

# corrplot
pdf(paste0(output_wd, "corrplot_light_scales.pdf"))
corrplot(cor(light),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of light across scales",
         mar = c(4,1,3,1))
corrplot(cor(light),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_light_scales.pdf"))
par(mfrow = c(3,2))
hist(light$light.100m, xlim = c(0,255), ylim = c(0, nrow(light)), 
     main = "Light. Histogram of values in 100m scale")
hist(light$light.500m, xlim = c(0,255), ylim = c(0, nrow(light)), 
     main = "Light. Histogram of values in 500m scale")
hist(light$light.1km, xlim = c(0,255), ylim = c(0, nrow(light)), 
     main = "Light. Histogram of values in 1km scale")
hist(light$light.2km, xlim = c(0,255), ylim = c(0, nrow(light)), 
     main = "Light. Histogram of values in 2km scale")
hist(light$light.5km, xlim = c(0,255), ylim = c(0, nrow(light)), 
     main = "Light. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOISE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
noise <- cbind.data.frame(noise.100m = plots.100m.all$noise,
                          noise.500m = plots.500m.all$noise,
                          noise.1km = plots.1km.all$noise,
                          noise.2km = plots.2km.all$noise,
                          noise.5km = plots.5km.all$noise)
cor(noise)
# corrplot
pdf(paste0(output_wd, "corrplot_noise_scales.pdf"))
corrplot(cor(noise),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of noise across scales",
         mar = c(4,1,3,1))
corrplot(cor(noise),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_noise_scales.pdf"))
par(mfrow = c(3,2))
hist(noise$noise.100m, xlim = c(0,100), ylim = c(0, nrow(noise)), 
     main = "Noise. Histogram of values in 100m scale")
hist(noise$noise.500m, xlim = c(0,100), ylim = c(0, nrow(noise)), 
     main = "Noise. Histogram of values in 500m scale")
hist(noise$noise.1km, xlim = c(0,100), ylim = c(0, nrow(noise)), 
     main = "Noise. Histogram of values in 1km scale")
hist(noise$noise.2km, xlim = c(0,100), ylim = c(0, nrow(noise)), 
     main = "Noise. Histogram of values in 2km scale")
hist(noise$noise.5km, xlim = c(0,100), ylim = c(0, nrow(noise)), 
     main = "Noise. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DAY TEMPERATURE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
temp <- cbind.data.frame(temp.100m = plots.100m.all$temp_day,
                          temp.500m = plots.500m.all$temp_day,
                          temp.1km = plots.1km.all$temp_day,
                          temp.2km = plots.2km.all$temp_day,
                          temp.5km = plots.5km.all$temp_day)

cor(temp)

# corrplot
pdf(paste0(output_wd, "corrplot_daytemp_scales.pdf"))
corrplot(cor(temp),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of day temperature across scales",
         mar = c(4,1,3,1))
corrplot(cor(temp),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_daytemp_scales.pdf"))
par(mfrow = c(3,2))
hist(temp$temp.100m, xlim = c(10,40), ylim = c(0, nrow(temp)), 
     main = "Temp. Histogram of values in 100m scale")
hist(temp$temp.500m, xlim = c(10,40), ylim = c(0, nrow(temp)), 
     main = "Temp. Histogram of values in 500m scale")
hist(temp$temp.1km, xlim = c(10,40), ylim = c(0, nrow(temp)), 
     main = "Temp. Histogram of values in 1km scale")
hist(temp$temp.2km, xlim = c(10,40), ylim = c(0, nrow(temp)), 
     main = "Temp. Histogram of values in 2km scale")
hist(temp$temp.5km, xlim = c(10,40), ylim = c(0, nrow(temp)), 
     main = "Temp. Histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TREE COVER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tree <- cbind.data.frame(tree.100m = plots.100m.all$tree_cover,
                         tree.500m = plots.500m.all$tree_cover,
                         tree.1km = plots.1km.all$tree_cover,
                         tree.2km = plots.2km.all$tree_cover,
                         tree.5km = plots.5km.all$tree_cover)

cor(tree)

# corrplot
pdf(paste0(output_wd, "corrplot_treecover_scales.pdf"))
corrplot(cor(tree),
         type = "upper",
         method = "number",
         tl.pos = "td",
         title = "corrplot of tree cover across scales",
         mar = c(4,1,3,1))
corrplot(cor(tree),
         type = "lower",
         method = "circle",
         tl.pos = "l",
         add = TRUE)
dev.off()

# Histograms of values
pdf(paste0(output_wd, "histograms_treecover_scales.pdf"))
par(mfrow = c(3,2))
hist(tree$tree.100m, xlim = c(0,100), ylim = c(0, nrow(tree)), 
     main = "Tree cover. histogram of values in 100m scale")
hist(tree$tree.500m, xlim = c(0,100), ylim = c(0, nrow(tree)), 
     main = "Tree cover. histogram of values in 500m scale")
hist(tree$tree.1km, xlim = c(0,100), ylim = c(0, nrow(tree)), 
     main = "Tree cover. histogram of values in 1km scale")
hist(tree$tree.2km, xlim = c(0,100), ylim = c(0, nrow(tree)), 
     main = "Tree cover. histogram of values in 2km scale")
hist(tree$tree.5km, xlim = c(0,100), ylim = c(0, nrow(tree)), 
     main = "Tree cover. histogram of values in 5km scale")
dev.off()
par(mfrow = c(1,1))

