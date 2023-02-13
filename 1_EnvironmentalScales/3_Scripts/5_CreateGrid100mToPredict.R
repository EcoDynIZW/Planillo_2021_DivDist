## Script that extracts the environmental variables at a 100m grid for future predictions

### Packages
library(raster)
library(dplyr)

### Workspace
work_wd <- getwd()
procdata_wd <- paste0(work_wd, "/4_Processed_data/")
env_wd <- paste0(work_wd, "/4_Processed_data/EnvValues/")


### Extract environmental values in a 100m grid for future predictions 

# we are going to extract values from 100m and 500m focal means, at 100m scale 
## load data
focal.100m <- stack(paste0(procdata_wd, "FocalMean100m_Stack.tif"))
names.focal <- read.csv(paste0(procdata_wd, "names_focal_stack.csv"))
names(focal.100m) <- names.focal$x

pop.100m <- raster(paste0(procdata_wd, "/Focal_pop_100m.tif"))

focal.500m <- stack(paste0(procdata_wd, "FocalMean500m_Stack.tif"))
names(focal.500m) <- names.focal$x

pop.500m <- raster(paste0(procdata_wd, "/Focal_pop_500m.tif"))


# extract raster as template
r <- raster(focal.100m$temp_day)
res(r) # 20 x 20m

# make it 100 m resolution 
r1 <- aggregate(r, fact = 5)
res(r1)

# Now get the coordinates of each pixel
p <- as(r1, 'SpatialPixels')
p1 <- as.data.frame(p)

## extract the variables that were at 100m scale
test.100m <- extract(focal.100m, p1, na.rm = TRUE)
pop.pred100 <- extract(pop.100m, p1, na.rm = TRUE)
grid.100m <- as.data.frame(cbind(p1, pop.pred100, test.100m))
head(grid.100m)

range(grid.100m$tree_cover)

## extract the variables that were at 500m scale
test.500m <- extract(focal.500m, p1, na.rm = TRUE)
pop.pred500 <- extract(pop.500m, p1, na.rm = TRUE)
grid.500m <- as.data.frame(cbind(p1, pop.pred500, test.500m))
head(grid.500m)


grid.100m.names <- grid.100m %>% 
  rename_all(.funs = ~ paste0(.x, "_100m"))
grid.500m.names <- grid.500m %>% 
  rename_all(.funs = ~ paste0(.x, "_500m"))

nrow(grid.100m.names)
nrow(grid.500m.names)

grid.final <- cbind.data.frame(grid.100m.names, grid.500m.names)

write.csv(grid.final, paste0(env_wd, "grid_topredict_100mres.csv"), row.names = FALSE)
