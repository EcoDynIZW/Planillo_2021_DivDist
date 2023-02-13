##################################
## PREPARE ENVIRONMENTAL LAYERS ##
##################################

## All spatial information uses epgs 25833

### Packages
library(raster)
library(fasterize)
library(sf)
library(dplyr)

### Workspace
work_wd <- getwd()
data_wd <- paste0(work_wd, "/2_Raw_data/Environment/")
output_wd <- paste0(work_wd, "/4_Processed_data/")


### Load layers vector layers to rasterize them
# Vector: waterbody polygons
water.poly <- read_sf(paste0(data_wd, "waterbodies_Berlin_25833.gpkg"))
plot(water.poly)

# Vector: temperature polygons
temp.poly <- read_sf(paste0(data_wd, "temperature_Berlin_2015_25833.gpkg"))
plot(temp.poly["T2M14HMEA"])

# Raster: imperviousness, 20m res as template raster
imperv <- raster(paste0(data_wd, "impervious_surface_20m_25833.tif"))
plot(imperv)


#############################
### Rasterize polygons

# create raster for template
my.raster <- raster(imperv)

## Transform water polygon layer into raster with presence/absence (1-0)
water.raster <- rasterize(water.poly, my.raster, field = 1, background = 0)
plot(water.raster)
# Save new raster
writeRaster(water.raster, filename = paste0(data_wd, "waterbodies_20m_25833.tif"))

## Transform temperature polygon layer into raster with temp at different times
temp.day.raster <- fasterize(temp.poly, my.raster, field = "T2M14HMEA") # day temperature measured at 14:00
plot(temp.day.raster)

temp.night.raster <- fasterize(temp.poly, my.raster, field = "T2M04HMEA") # night temp measured at  04:00
plot(temp.night.raster)

temp.cooldown.raster <- fasterize(temp.poly, my.raster, field = "ABKUEHLMEA") # cooldown index day-night
plot(temp.cooldown.raster)

# save rasters
writeRaster(temp.day.raster, filename = paste0(data_wd, "temp_day_20m_25833.tif"), overwrite = TRUE)
writeRaster(temp.night.raster, filename = paste0(data_wd, "temp_night_20m_25833.tif"), overwrite = TRUE)
writeRaster(temp.cooldown.raster, filename = paste0(data_wd, "temp_cooldown_20m_25833.tif"), overwrite = TRUE)


##################################
## Plot environmental values maps  
##################################






#######################
## Compute focal means  
#######################


## Obtained mean values at different scales using a moving window approach without 
# transforming the resolution of the rasters

### Load data: all rasters as a stack
all.rasters <- list.files(data_wd, pattern = "tif$")
# we remove population density and water bodies because we will use them in a different way
my.rasters <- all.rasters[!all.rasters %in% 
                            c("population_density_20m_25833.tif", "waterbodies_20m_25833.tif")]
multi.raster <- stack(paste0(data_wd, my.rasters))


# population density as independent raster to get the sum instead of mean
pop.density <- raster(paste0(data_wd, "population_density_20m_25833.tif"))


##### Focal mean: POPULATION DENSITY #####

# Each pixel represents people in 20x20 m. To compute the population in different areas, we should sum the pixels. 
# But focal function does the mean, so we are going to cheat it by recalculating the value of 
# each pixel as population density in the area of interest

## 100 m
# Proportion new pixel/original pixel for multiplying the layer
(100*100)/(20*20)
# [1] 25
pop.dens.100m <- pop.density * 25

#Relation new/original pixel for the matrix in focal means
100/20
# [1] 5

focal.pop.100m <- focal(pop.dens.100m, matrix(1/(5*5), nrow =  5, ncol = 5),
                        na.rm = TRUE , pad=TRUE, padValue=0)
plot(focal.pop.100m)

writeRaster(focal.pop.100m, filename = paste0(output_wd,  "Focal_pop_100m.tif"))

## 500 m
(500*500)/(20*20)
# [1] 625
pop.dens.500m <- pop.density * 625

500/20
# [1] 25

focal.pop.500m <- focal(pop.dens.500m, matrix(1/(25*25), nrow =  25, ncol = 25),
                        na.rm = TRUE , pad=TRUE, padValue=0)
plot(focal.pop.500m)
writeRaster(focal.pop.500m, filename = paste0(output_wd,  "Focal_pop_500m.tif"))

## 1 KM
(1000*1000)/(20*20)
# [1] 2500
pop.dens.1km <- pop.density * 2500

1000/20
# [1] 50

focal.pop.1km <- focal(pop.dens.1km, matrix(1/(51*51), nrow =  51, ncol = 51),
                       na.rm = TRUE , pad=TRUE, padValue=0)

plot(focal.pop.1km)
writeRaster(focal.pop.1km, filename = paste0(output_wd,  "Focal_pop_1km.tif"))

### 2km
(2000*2000)/(20*20)
# [1] 10000
pop.dens.2km <- pop.density * 10000

2000/20
# [1] 100

focal.pop.2km <- focal(pop.dens.2km, matrix(1/(101*101), nrow =  101, ncol = 101),
                       na.rm = TRUE , pad=TRUE, padValue=0)
plot(focal.pop.2km)
writeRaster(focal.pop.2km, filename = paste0(output_wd,  "Focal_pop_2km.tif"))

### 5km
(5000*5000)/(20*20)
# [1] 62500
5000/20
# [1] 250

pop.dens.5km <- pop.density * 62500
focal.pop.5km <- focal(pop.dens.5km, matrix(1/(251*251), nrow =  251, ncol = 251),
                       na.rm = TRUE , pad=TRUE, padValue=0)
plot(focal.pop.5km)
writeRaster(focal.pop.5km, filename = paste0(output_wd,  "Focal_pop_5km.tif"))



##### Focal mean: RASTER STACK #####

## function to run focal mean in stack
focal_mean_for_stack <- function(x, y, z) {
  lapply(1:nlayers(x), function(i){
    scratch <-focal(x[[i]], w=matrix(1/(y*z), nrow=y, ncol=z), 
                    na.rm = TRUE , pad=TRUE, padValue=0); 
    return(scratch) })  
}

###100 M

## specify the dimension of the matrix to cumpute the focal mean 
# 5 x 5 pixels : 100 x 100 m grid
number_rows <- 5 
number_cols <- 5

multi.raster %>% 
  focal_mean_for_stack(number_rows,number_cols) %>% 
  stack() -> my.result.100m

## rename the raster layer
names(my.result.100m) <- names(multi.raster)

## plot
plot(my.result.100m)

## export results in tif format but without names
writeRaster(my.result.100m, filename=paste0(output_wd, "FocalMean100m_Stack.tif"), 
            options="INTERLEAVE=BAND", overwrite=TRUE)
names <- c("distance_water", "impervious_surface", "light", "noise", "open_green", 
           "temp_cooldown", "temp_day", "temp_night",  "tree_cover")
write.csv(names, paste0(output_wd, "names_focal_stack.csv"), row.names = F)


### FOCAL MEAN 500 M

# 500/20 = 25 x 25 pixels : 500 x 500 m grid
number_rows <- 25 
number_cols <- 25

multi.raster %>% focal_mean_for_stack(number_rows,number_cols) %>% stack() -> my.result.500m
names(my.result.500m) <- names(multi.raster)

plot(my.result.500m)

writeRaster(my.result.500m, filename=paste0(output_wd, "FocalMean500m_Stack.tif"), 
            options="INTERLEAVE=BAND", overwrite=TRUE) # same names as 100m results


### FOCAL MEAN 1km

number_rows <- 51 
number_cols <- 51

multi.raster %>% focal_mean_for_stack(number_rows,number_cols) %>% stack() -> my.result.1km

names(my.result.1km) <- names(multi.raster)

plot(my.result.1km)
writeRaster(my.result.1km, filename=paste0(output_wd, "FocalMean1km_Stack.tif"), 
            options="INTERLEAVE=BAND", overwrite=TRUE) # same names as 100m results


### FOCAL MEAN 2km

number_rows <- 101 
number_cols <- 101

multi.raster %>% focal_mean_for_stack(number_rows,number_cols) %>% stack() -> my.result.2km

names(my.result.2km) <- names(multi.raster)

plot(my.result.2km)

writeRaster(my.result.2km, filename=paste0(output_wd, "FocalMean2km_Stack.tif"), 
            options="INTERLEAVE=BAND", overwrite=TRUE) # same names as 100m results


### FOCAL MEAN 5km

## Too slow if done as a stack. Rasters one by one

# 5000/20 = 250 x 250 pixels : 5000 x 5000 m grid

y <- 251 
z <- 251

names <- c("distance_water", "impervious_surface", "light", "noise", "open_green", 
           "temp_cooldown", "temp_day", "temp_night",  "tree_cover")
multi.raster

dist.wtr <- multi.raster[[1]]
plot(dist.wtr)

imperv <- multi.raster[[2]]
light <- multi.raster[[3]]
noise <- multi.raster[[4]]
open.green <- multi.raster[[5]]
temp_cool <- multi.raster[[6]]
temp_day <- multi.raster[[7]]
temp_night <- multi.raster[[8]]
tree.cover <- multi.raster[[9]]


focal5km.distwtr <- focal(dist.wtr,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                          na.rm = TRUE , pad=TRUE, padValue=0)
plot(focal5km.distwtr)

focal5km.imperv <- focal(imperv,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                         na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.light <- focal(light,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                        na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.noise <- focal(noise,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                        na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.opengreen <- focal(open.green,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                            na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.temp_cool <- focal(temp_cool,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                            na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.temp_day <- focal(temp_day,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                           na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.temp_night <- focal(temp_night,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                             na.rm = TRUE , pad=TRUE, padValue=0)

focal5km.treecover <- focal(tree.cover,  w=matrix(1/(y*z), nrow =  y, ncol = z),
                            na.rm = TRUE , pad=TRUE, padValue=0)

# make them a stack again
my.results.5km <- stack(focal5km.distwtr, focal5km.imperv, focal5km.light, focal5km.noise, focal5km.opengreen, 
                        focal5km.temp_cool, focal5km.temp_day, focal5km.temp_night, focal5km.treecover)

writeRaster(my.results.5km, filename=paste0(output_wd, "FocalMean5km_Stack.tif"), 
            options="INTERLEAVE=BAND", overwrite=TRUE) 




####################################################
### EXTRACT ENVIRONMENTAL VALUES FOR TRANSECTS AND PLOTS
####################################################


### Packages
library(raster)
library(sf)
library(dplyr)
extract = raster::extract # set the function extract of raster package as priority

### WORK SPACE
work_wd <- getwd()
data_wd <- paste0(work_wd, "/2_Raw_data/Environment/")
sp_wd <- paste0(work_wd, "/2_Raw_data/Species/")
fmeans_wd <- paste0(work_wd, "/4_Processed_data/")
env_wd <- paste0(work_wd, "/4_Processed_data/EnvValues/")
output_wd <- paste0(work_wd, "/5_Results/")



### Load data

all.plots <- read_sf(paste0(sp_wd, "Dry_Grassland_Plots_Berlin_25833.gpkg")) # all BIBS plots
selected.plots <- read.csv(paste0(sp_wd, "Invertebrates_shared_plots.csv")) # Plots used in the arthropod analyses

plots <- all.plots %>% 
  dplyr::filter(ID_plot %in% selected.plots$Plot_ID)
plot(st_geometry(plots), main = "Arthropod plots")


transects <- read_sf(paste0(sp_wd, "birds_transects_final_25833.gpkg")) # bird transects
plot(st_geometry(transects), main = "Bird transects")

focal.100m <- stack(paste0(fmeans_wd, "FocalMean100m_Stack.tif")) # Focal means at 100m
names.focal <- read.csv(paste0(fmeans_wd, "names_focal_stack.csv"))
names(focal.100m) <- names.focal$x
plot(focal.100m)

plot(focal.100m$temp_night, 
     main = "temperature at night in summer")

focal.500m <- stack(paste0(fmeans_wd, "FocalMean500m_Stack.tif")) # Focal means at 500m
names(focal.500m) <- names.focal$x

focal.1km <- stack(paste0(fmeans_wd, "FocalMean1km_Stack.tif")) # Focal means at 1km
names(focal.1km) <- names.focal$x

focal.2km <- stack(paste0(fmeans_wd, "FocalMean2km_Stack.tif")) # Focal means at 2km
names(focal.2km) <- names.focal$x

focal.5km <- stack(paste0(fmeans_wd, "FocalMean5km_Stack.tif")) # Focal means at 5km
names(focal.5km) <- names.focal$x

pop.100m <- raster(paste0(fmeans_wd, "Focal_pop_100m.tif")) # population density raster at 100m 
pop.500m <- raster(paste0(fmeans_wd, "Focal_pop_500m.tif")) # population density raster at 500m 
pop.1km <- raster(paste0(fmeans_wd, "Focal_pop_1km.tif")) # population density raster at 1km 
pop.2km <- raster(paste0(fmeans_wd, "Focal_pop_2km.tif")) # population density raster at 2km 
pop.5km <- raster(paste0(fmeans_wd, "Focal_pop_5km.tif")) # population density raster at 5km 




###############################################################
## Extract values of environmental layers for BIBS PLOTS 
###############################################################

# Raster stack values
plots.100m <- extract(focal.100m, plots, fun = mean, na.rm = TRUE)
row.names(plots.100m) <- plots$ID_plot

plots.500m <- extract(focal.500m, plots, fun = mean, na.rm = TRUE)
row.names(plots.500m) <- plots$ID_plot

plots.1km <- extract(focal.1km, plots, fun = mean, na.rm = TRUE)
row.names(plots.1km) <- plots$ID_plot

plots.2km <- extract(focal.2km, plots, fun = mean, na.rm = TRUE)
row.names(plots.2km) <- plots$ID_plot

plots.5km <- extract(focal.5km, plots, fun = mean, na.rm = TRUE)
row.names(plots.5km) <- plots$ID_plot


# Human population density values 
plots.pop100m <- extract(pop.100m, plots)
plots.pop100m <- as.data.frame(cbind(ID_plot = plots$ID_plot, pop.100m = plots.pop100m))

plots.pop500m <- extract(pop.500m, plots)
plots.pop500m <- as.data.frame(cbind(ID_plot = plots$ID_plot, pop.500m = plots.pop500m))

plots.pop1km <- extract(pop.1km, plots)
plots.pop1km <- as.data.frame(cbind(ID_plot = plots$ID_plot, pop.1km = plots.pop1km))

plots.pop2km <- extract(pop.2km, plots)
plots.pop2km <- as.data.frame(cbind(ID_plot = plots$ID_plot, pop.2km = plots.pop2km))

plots.pop5km <- extract(pop.5km, plots)
plots.pop5km <- as.data.frame(cbind(ID_plot = plots$ID_plot, pop.5km = plots.pop5km))

# save extracted values
cbind.data.frame(rownames(plots.100m), plots.pop100m$ID_plot)
as.data.frame(cbind(plots.100m, as.numeric(plots.pop100m$pop.100m)))

plots.100m.all <- as.data.frame(cbind(pop_100m = as.numeric(plots.pop100m$pop.100m), plots.100m))
plots.500m.all <- as.data.frame(cbind(pop_500m = as.numeric(plots.pop500m$pop.500m), plots.500m))
plots.1km.all <- as.data.frame(cbind(pop_1km = as.numeric(plots.pop1km$pop.1km), plots.1km))
plots.2km.all <- as.data.frame(cbind(pop_2km = as.numeric(plots.pop2km$pop.2km), plots.2km))
plots.5km.all <- as.data.frame(cbind(pop_5km = as.numeric(plots.pop5km$pop.5km), plots.5km))

write.csv(plots.100m.all, paste0(env_wd, "bibsplots_allenvir_100m.csv"))
write.csv(plots.500m.all, paste0(env_wd, "bibsplots_allenvir_500m.csv"))
write.csv(plots.1km.all, paste0(env_wd, "bibsplots_allenvir_1km.csv"))
write.csv(plots.2km.all, paste0(env_wd, "bibsplots_allenvir_2km.csv"))
write.csv(plots.5km.all, paste0(env_wd, "bibsplots_allenvir_5km.csv"))


###############################################################
## Extract values of environmental layers for BIRD TRANSECTS
###############################################################



# Raster stack values
transects.100m <- extract(focal.100m, transects, fun = mean, na.rm = TRUE)
row.names(transects.100m) <- transects$routcode_2

transects.500m <- extract(focal.500m, transects, fun = mean, na.rm = TRUE)
row.names(transects.500m) <- transects$routcode_2

transects.1km <- extract(focal.1km, transects, fun = mean, na.rm = TRUE)
row.names(transects.1km) <- transects$routcode_2

transects.2km <- extract(focal.2km, transects, fun = mean, na.rm = TRUE)
row.names(transects.2km) <- transects$routcode_2

transects.5km <- extract(focal.5km, transects, fun = mean, na.rm = TRUE)
row.names(transects.5km) <- transects$routcode_2


# Human population density values 
transects.pop100m <- extract(pop.100m, transects, fun = mean, na.rm = TRUE)
transects.pop100m <- as.data.frame(cbind(routcode_2 = transects$routcode_2, pop.100m = transects.pop100m))

transects.pop500m <- extract(pop.500m, transects, fun = mean, na.rm = TRUE)
transects.pop500m <- as.data.frame(cbind(routcode_2 = transects$routcode_2, pop.500m = transects.pop500m))

transects.pop1km <- extract(pop.1km, transects, fun = mean, na.rm = TRUE)
transects.pop1km <- as.data.frame(cbind(routcode_2 = transects$routcode_2, pop.1km = transects.pop1km))

transects.pop2km <- extract(pop.2km, transects, fun = mean, na.rm = TRUE)
transects.pop2km <- as.data.frame(cbind(routcode_2 = transects$routcode_2, pop.2km = transects.pop2km))

transects.pop5km <- extract(pop.5km, transects, fun = mean, na.rm = TRUE)
transects.pop5km <- as.data.frame(cbind(routcode_2 = transects$routcode_2, pop.5km = transects.pop5km))


# save extracted values
cbind.data.frame(rownames(transects.100m), transects.pop100m$routcode_2)
transects.100m.all <- as.data.frame(cbind(pop_100m = as.numeric(transects.pop100m$pop.100m), 
                                          transects.100m))
transects.500m.all <- as.data.frame(cbind(pop_500m = as.numeric(transects.pop500m$pop.500m), 
                                          transects.500m))
transects.1km.all <- as.data.frame(cbind(pop_1km = as.numeric(transects.pop1km$pop.1km), 
                                         transects.1km))
transects.2km.all <- as.data.frame(cbind(pop_2km = as.numeric(transects.pop2km$pop.2km), 
                                         transects.2km))
transects.5km.all <- as.data.frame(cbind(pop_5km = as.numeric(transects.pop5km$pop.5km), 
                                         transects.5km))

write.csv(transects.100m.all, paste0(env_wd, "transects_allenvir_100m.csv"))
write.csv(transects.500m.all, paste0(env_wd, "transects_allenvir_500m.csv"))
write.csv(transects.1km.all, paste0(env_wd, "transects_allenvir_1km.csv"))
write.csv(transects.2km.all, paste0(env_wd, "transects_allenvir_2km.csv"))
write.csv(transects.5km.all, paste0(env_wd, "transects_allenvir_5km.csv"))




