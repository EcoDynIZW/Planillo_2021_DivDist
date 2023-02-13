####################################################
# Creating pretty maps of the environmental values #
# Jul 19, A. Planillo
####################################################


### Packages
library(tmap)
library(sf)
library(raster)
library(RColorBrewer)
library(viridis)
library(dplyr)


### Workspace
work_wd <- getwd()
data_wd <- paste0(work_wd, "/2_Raw_data/Environment/")
procdata_wd <- paste0(work_wd, "/4_Processed_data/")
output_wd <- file.path(work_wd, "5_Results/Mymaps_envValues")



### Data
# Berlin border
berlin <- read_sf(paste0(Rawmaps_wd, "/berlin_border.gpkg"))
plot(st_geometry(berlin))


berlin <- berlin %>% st_set_crs(NA) %>% st_set_crs(3035)
berlin <- st_transform(berlin, crs = 25833)


# environment at 100 m resolution
focal.100m <- stack(paste0(procdata_wd, "/FocalMean100m_Stack.tif"))
names.focal <- read.csv(paste0(procdata_wd, "/names_focal_stack.csv"))
names(focal.100m) <- names.focal$x

focal.100m


# Human population density
pop.100m <- raster(paste0(procdata_wd, "/Focal_pop_100m.tif"))
plot(pop.100m)


### Maps

##### Natural environment ####
tmap_mode("plot")


# temp of day during summer

temp <- subset(focal.100m, "temp_day", drop = FALSE)
temp_map <- tm_shape(temp) +
  tm_raster("temp_day", title = "Day temperature (ºC)",
            palette = rev(heat.colors(6)),
            breaks = c(-Inf, 15, 20, 25, 30, 35),
            legend.show = TRUE) +
  tm_shape(berlin) +
  tm_borders("grey20", lwd = 2.5) +
  # tm_text("spatial_al") 
  tm_scale_bar(position=c("left", "bottom")) +
  # tm_layout(frame = FALSE) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            # legend.outside.position = "right") +
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

tmap_save(temp_map, paste0(output_wd, "/tempday_100m.png"), 
          dpi = 600)


# Tree cover
# use the mask function to create a raster for berlin
tree.cover <- mask(focal.100m$tree_cover, berlin)

tree_map <- tm_shape(tree.cover) +
  tm_raster("tree_cover", title = "Tree cover (%)",
            palette = "Greens")+
  # palette = rev(terrain.colors(5)),
  # breaks = c(-Inf, 25, 50, 75, 100)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  # tm_text("spatial_al") 
  # tm_shape(water) +
  # tm_polygons("blue") +
  tm_scale_bar(position=c("left", "bottom")) +
  # tm_layout(frame = FALSE) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            # legend.outside.position = "right") +
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

tree_map
tmap_save(tree_map, paste0(output_wd, "/tree_100m.png"), 
          dpi = 600)


# Open green area
# use the mask function to create a raster for berlin
open.green <- mask(focal.100m$open_green, berlin)

# Define colors for plotting
my.palette.ogreen <- c("khaki", "darkgreen", "greenyellow", "yellowgreen")

ogreen_map <- tm_shape(open.green) +
  tm_raster("open_green", title = "Green area (%)",
            palette = "YlOrBr") +
  # palette = rev(terrain.colors(5)),
  # breaks = c(-Inf, 25, 50, 75, 100)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  # tm_text("spatial_al") 
  # tm_shape(water) +
  # tm_polygons("blue") +
  tm_scale_bar(position=c("left", "bottom")) +
  # tm_layout(frame = FALSE) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            # legend.outside.position = "right") +
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

tmap_save(ogreen_map, paste0(output_wd, "/open_green_100m.png"), 
          dpi = 600)


# distance to water

# use the mask function to create a raster for berlin
dist.water <- mask(focal.100m$distance_water, berlin)

dwater_map <- tm_shape(dist.water) +
  tm_raster("distance_water", title = "Distance to water (m)",
            palette = "-Blues",
            breaks = c(-Inf, 100, 500, 1000, Inf)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

tmap_save(dwater_map, paste0(output_wd, "/dist_water_100m.png"), 
          dpi = 600)



##### Anthropogenic environment ####

# Impervious surface

# use the mask function to create a raster for berlin
imperv <- mask(focal.100m$impervious_surface, berlin)

imperv_map <- tm_shape(imperv) +
  tm_raster("impervious_surface", title = "Impervious surface (%)",
            palette = "Greys") +
  # breaks = c(-Inf, 100, 500, 1000, Inf)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

tmap_save(imperv_map, paste0(output_wd, "/impervious_surface_100m.png"), 
          dpi = 600)


# Noise

noise <- mask(focal.100m$noise, berlin)


noise_map <- tm_shape(noise) +
  tm_raster("noise", title = "Noise (dBa)",
            palette = "inferno",
            breaks = c(-Inf, 40, 60, 70, 80, 90)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

noise_map
tmap_save(noise_map, paste0(output_wd, "/noise_100m.png"), 
          dpi = 600)



# Human population

pop_map <- tm_shape(pop.100m) +
  tm_raster("Focal_pop_100m", title = "Number of inhabitants",
            palette = "YlOrRd",
            breaks = c(-Inf, 50, 100, 200, 500, Inf)) +
  tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

pop_map
tmap_save(pop_map, paste0(output_wd, "/population_100m.png"), 
          dpi = 600)





##################################
## Maps of the observations ####

species_wd <- file.path(WorkDir, "2_Raw_data/Species")

list.files(species_wd)
bibsplots <- read_sf(paste0(species_wd, "/Dry_Grassland_Plots_Berlin_25833.gpkg"))
transects2 <- read_sf(paste0(species_wd, "/birds_transects_final_25833.gpkg"))

plot(st_geometry(berlin))
plot(st_geometry(bibsplots), col = "blue", add = TRUE)
plot(st_geometry(transects), add = TRUE)
plot(st_geometry(transects2), add = TRUE)


# Maps
bibsplots_map <- tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_shape(bibsplots) +
  tm_symbols(col = "darkgreen", alpha = 1, 
             border.lwd = 1, border.col = "black",  
             size = 0.5) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")
bibsplots_map

transects_lines_map <- tm_shape(berlin) +
  tm_borders("gray20", lwd = 2.5) +
  tm_shape(transects2) +
  tm_lines(col = "darkred", alpha = 1, 
           # border.lwd = 1, border.col = "black",  
           lwd = 2) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_layout(frame = TRUE, compass.type = "4star", legend.outside = FALSE,
            inner.margins = c(0.05,0.03,0.05,0.03),
            legend.position = c("right", "top"))+
  tm_compass(position = c("left", "top")) +
  tm_legend(frame = TRUE,
            bg.color="lightyellow")

transects_lines_map
