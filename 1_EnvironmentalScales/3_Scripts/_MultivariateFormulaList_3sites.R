#####################
#### LIST OF MODEL FORMULAS FOR SCALE AND VARIABLE SELECTION
#####################


### BIBS PLOTS
## Formula list for BIBS plots model selection : 42 plots

## We remove singletons or doubletons species (species in only one or two sites)

## the formulas used for bibs plots multivariate model selection:
## Keeping same variables, only changing scale

# FORMULAS BIBS PLOTS

# Null model
form.null <- mvformula(myresponse ~ 1)

# plant abundance and richness
form.plant <- mvformula(myresponse ~ plants$richness + plants$abund)

# Scale: 100m
form.100m.full <- mvformula(myresponse ~ env.100m.cs$impervious_surface + env.100m.cs$noise + 
                              env.100m.cs$distance_water + env.100m.cs$temp_day + env.100m.cs$open_green)
form.100m.human <- mvformula(myresponse ~ env.100m.cs$impervious_surface + env.100m.cs$noise)
form.100m.nat <- mvformula(myresponse ~ env.100m.cs$distance_water + env.100m.cs$temp_day + env.100m.cs$open_green)

# Scale: 500m
form.500m.full <- mvformula(myresponse ~ env.500m.cs$impervious_surface + env.500m.cs$noise + 
                              env.500m.cs$distance_water + env.500m.cs$temp_day + env.500m.cs$open_green)
form.500m.human <- mvformula(myresponse ~ env.500m.cs$impervious_surface + env.500m.cs$noise)
form.500m.nat <- mvformula(myresponse ~ env.500m.cs$distance_water + env.500m.cs$temp_day + env.500m.cs$open_green)

# Scale: 1km
form.1km.full <- mvformula(myresponse ~ env.1km.cs$impervious_surface + env.1km.cs$noise + 
                              env.1km.cs$distance_water + env.1km.cs$temp_day + env.1km.cs$open_green)
form.1km.human <- mvformula(myresponse ~ env.1km.cs$impervious_surface + env.1km.cs$noise)
form.1km.nat <- mvformula(myresponse ~ env.1km.cs$distance_water + env.1km.cs$temp_day + env.1km.cs$open_green)

# Scale: 2km
form.2km.full <- mvformula(myresponse ~ env.2km.cs$impervious_surface + env.2km.cs$noise + 
                              env.2km.cs$distance_water + env.2km.cs$temp_day + env.2km.cs$open_green)
form.2km.human <- mvformula(myresponse ~ env.2km.cs$impervious_surface + env.2km.cs$noise)
form.2km.nat <- mvformula(myresponse ~ env.2km.cs$distance_water + env.2km.cs$temp_day + env.2km.cs$open_green)

# Scale: 5km
form.5km.full <- mvformula(myresponse ~ env.5km.cs$impervious_surface + env.5km.cs$noise + 
                              env.5km.cs$distance_water + env.5km.cs$temp_day + env.5km.cs$open_green)
form.5km.human <- mvformula(myresponse ~ env.5km.cs$impervious_surface + env.5km.cs$noise)
form.5km.nat <- mvformula(myresponse ~ env.5km.cs$distance_water + env.5km.cs$temp_day + env.5km.cs$open_green)

formulas.42plots.3sites <- as.list(c(form.null, form.plant, form.100m.full, form.100m.human, form.100m.nat,
                              form.500m.full, form.500m.human, form.500m.nat, form.1km.full, form.1km.human, form.1km.nat,
                              form.2km.full, form.2km.human, form.2km.nat, form.5km.full, form.5km.human, form.5km.nat))
names.form42plots <- c("form.null", "form.plant", "form.100m.full", "form.100m.human", "form.100m.nat",
                       "form.500m.full", "form.500m.human", "form.500m.nat", 'form.1km.full', "form.1km.human", "form.1km.nat",
                       "form.2km.full", "form.2km.human", "form.2km.nat", 'form.5km.full', "form.5km.human", "form.5km.nat")
for (i in 1:length(formulas.42plots.3sites)){
  names(formulas.42plots.3sites)[i] <- names.form42plots[i]
}
formulas.42plots.3sites

saveRDS(formulas.42plots.3sites, file=paste0("./4_Processed_data/formulas_42plots_3sites.rds"))




#####################
#### BIRDS TRANSECTS
#####################


# FORMULAS Birds TRANSECTS

# Null model
form.bird.null <- mvformula(myresponse ~ 1)

# Scale: 100m
form.bird.100m.full <- mvformula(myresponse ~ trans.100m.cs$tree_cover + trans.100m.cs$open_green + trans.100m.cs$distance_water +
                              trans.100m.cs$noise + trans.100m.cs$pop_100m)
form.bird.100m.human <- mvformula(myresponse ~ trans.100m.cs$noise + trans.100m.cs$pop_100m)
form.bird.100m.nat <- mvformula(myresponse ~ trans.100m.cs$tree_cover + trans.100m.cs$open_green + trans.100m.cs$distance_water)

# Scale: 500m
form.bird.500m.full <- mvformula(myresponse ~ trans.500m.cs$tree_cover + trans.500m.cs$open_green + trans.500m.cs$distance_water +
                              trans.500m.cs$noise + trans.500m.cs$pop_500m)
form.bird.500m.human <- mvformula(myresponse ~ trans.500m.cs$noise + trans.500m.cs$pop_500m)
form.bird.500m.nat <- mvformula(myresponse ~ trans.500m.cs$tree_cover + trans.500m.cs$open_green + trans.500m.cs$distance_water)

# Scale: 1km
form.bird.1km.full <- mvformula(myresponse ~ trans.1km.cs$tree_cover + trans.1km.cs$open_green + trans.1km.cs$distance_water +
                              trans.1km.cs$noise + trans.1km.cs$pop_1km)
form.bird.1km.human <- mvformula(myresponse ~ trans.1km.cs$noise + trans.1km.cs$pop_1km)
form.bird.1km.nat <- mvformula(myresponse ~ trans.1km.cs$tree_cover + trans.1km.cs$open_green + trans.1km.cs$distance_water)

# Scale: 2km
form.bird.2km.full <- mvformula(myresponse ~ trans.2km.cs$tree_cover + trans.2km.cs$open_green + trans.2km.cs$distance_water +
                             trans.2km.cs$noise + trans.2km.cs$pop_2km)
form.bird.2km.human <- mvformula(myresponse ~ trans.2km.cs$noise + trans.2km.cs$pop_2km)
form.bird.2km.nat <- mvformula(myresponse ~ trans.2km.cs$tree_cover + trans.2km.cs$open_green + trans.2km.cs$distance_water)

# Scale: 5km
form.bird.5km.full <- mvformula(myresponse ~ trans.5km.cs$tree_cover + trans.5km.cs$distance_water +
                             trans.5km.cs$noise + trans.5km.cs$pop_5km)
form.bird.5km.human <- mvformula(myresponse ~ trans.5km.cs$noise + trans.5km.cs$pop_5km)
form.bird.5km.nat <- mvformula(myresponse ~ trans.5km.cs$tree_cover + trans.5km.cs$distance_water)

formulas.birds <- as.list(c(form.bird.null, form.bird.100m.full, form.bird.100m.human, form.bird.100m.nat,
                              form.bird.500m.full, form.bird.500m.human, form.bird.500m.nat, form.bird.1km.full, form.bird.1km.human, form.bird.1km.nat,
                              form.bird.2km.full, form.bird.2km.human, form.bird.2km.nat, form.bird.5km.full, form.bird.5km.human, form.bird.5km.nat))
names.formbirds <- c("form.bird.null", "form.bird.100m.full", "form.bird.100m.human", "form.bird.100m.nat",
                       "form.bird.500m.full", "form.bird.500m.human", "form.bird.500m.nat", 'form.bird.1km.full', "form.bird.1km.human", "form.bird.1km.nat",
                       "form.bird.2km.full", "form.bird.2km.human", "form.bird.2km.nat", 'form.bird.5km.full', "form.bird.5km.human", "form.bird.5km.nat")
for (i in 1:length(formulas.birds)){
  names(formulas.birds)[i] <- names.formbirds[i]
}
formulas.birds

saveRDS(formulas.birds, file=paste0("./4_Processed_data/formulas_birs.rds"))
