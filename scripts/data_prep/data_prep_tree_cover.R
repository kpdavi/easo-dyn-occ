# Purpose: Compile tree cover data to include as covariates in the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Notes ----
## Anonymity of survey locations ----
# The locations of the monitoring sites are not public, so the longitude and latitude coordinates in the sites.csv were randomly generated. Thus, the values in the output file created in this script will not align with the data values used in the manuscript, but the code will work, and the data with the real covariate values are available on Zenodo.
# Output file for this script: processed/covs_tree.csv. That file isn't saved in this script because the covariate values aren't accurate (given that in this script they are generated for fake survey locations; see comment immediately above).

## Data downloads ----
# Tree cover rasters from the Rangeland Analysis Platform (https://rangelands.app/rap/?biomass_t=herbaceous&ll=39.0000,-98.0000&z=5) were downloaded on April 18, 2022. They were cropped to the extent of Fort Collins, CO in ArcGIS.


# Required packages ----
# install.packages("pacman")
pacman::p_load(
  tidyverse,
  sf,
  terra,
  landscapemetrics
  )


# Global variables ----
root_dir <- "/Users/kristindavis/Library/CloudStorage/OneDrive-NewMexicoStateUniversity/PhD/Data/EASO/publication_data/"
# root_dir <- "/Users/kristindavis/Library/CloudStorage/OneDrive-NewMexicoStateUniversity/PhD/Data/EASO/publication_data/"
buf_size <- 250  # for 250-m radii buffer used to calculate tree cover variables
year_begin <- 2012 # for calculating a 'year' column after calculating the cohesion metric


# Import data ----
## Transects ----
sites <- read.csv(paste0(root_dir, "raw/sites.csv"), header = TRUE)
sites$ID <- as.numeric(as.factor(sites$Transect))
sites_sp <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)

## Tree cover ----
tree_rasts <- terra::rast(paste0(root_dir, "raw/RAP_tree_cover_2013-2021.tif"))
plot(tree_rasts)

## Get unique values of tree cover, to check values are reasonable
tree_rast_values <- terra::unique(tree_rasts)
tree_rast_values_unique_by_year <- lapply(tree_rast_values, unique)
print(tree_rast_values_unique_by_year)


# Prepare data ----
# Create polygons using a 250-m-radii buffer ----
sites_poly <- sites_sp |>
  st_cast("POINT") |>
  st_buffer(dist = buf_size)

# Ensure sites and rasters are in the same coordinate reference system
sites_poly <- st_transform(sites_poly, crs = crs(tree_rasts))


# Calculate tree cover variables ----
## Average proportion of tree cover ----
tree_avg_extract <- terra::extract(
  x = tree_rasts, 
  y = sites_poly, 
  fun = mean, 
  method = "bilinear",
  ID = TRUE,
  raw = TRUE,
  na.rm = TRUE
  )

# Format output
tree_avg_extract_long <- tree_avg_extract |>
  pivot_longer(
    cols = matches("\\d{4}"), 
    names_to = "names", 
    values_to = "tree.prop"
    ) |>
  mutate(
    year = as.numeric(str_extract(names, "\\d{4}")),
    tree.prop = round(tree.prop, 3)
    ) |>
  dplyr::select(ID, year, tree.prop)

tree_avg <- sites |>
  dplyr::select(Transect, ID) |>
  left_join(tree_avg_extract_long, by = "ID") |>
  dplyr::select(-ID)


# Cohesion of tree cover ----
## Reclassify raster ----
### Binary layer where tree cover > 10% is 1 and tree cover <= 10% is 0
reclass_values <- c(0, 10, 0, 10, 100, 1)
reclass_mat <- matrix(reclass_values, ncol = 3, byrow = TRUE)
tree_rasts_reclass <- classify(tree_rasts, reclass_mat)
plot(tree_rasts_reclass)


## Check format of reclassified raster ----
### Rasters need to be in a metric coordinate reference system (crs) if units of results are based on cell sizes and/or distances
check_landscape(tree_rasts_reclass) 
crs(tree_rasts_reclass)  # not a metric crs 

### Reproject data into a metric crs
crs_utmzone_13 <- "26913" # EPSG code for NAD83 UTM Zone 13, where Fort Collins is located

sites_sp_utm13 <- st_transform(sites_sp, crs = paste0("epsg:", crs_utmzone_13))
sites_poly_utm13 <- st_transform(sites_poly, crs = paste0("epsg:", crs_utmzone_13))

tree_rasts_reclass_utm13 <- project(tree_rasts_reclass, y = paste0("epsg:", crs_utmzone_13), method = "near")
check_landscape(tree_rasts_reclass_utm13)  # passes check for metric crs


## Calculate cohesion ----
tree_coh_extract <- sample_lsm(
  tree_rasts_reclass_utm13, 
  y = sites_poly_utm13, 
  plot_id = sites_poly_utm13$Transect, 
  shape = "circle", 
  size = 250, 
  what = "lsm_c_cohesion"
  )

# Format output
tree_coh_form <- tree_coh_extract |>
  mutate(
    year = layer + year_begin,
    value = round(value, 3)
    ) |>
  rename(
    Transect = plot_id,
    tree.coh = value
    ) |>
  dplyr::select(Transect, year, tree.coh, class)

tree_coh <- sites |>
  dplyr::select(Transect) |>
  left_join(tree_coh_form, by = "Transect") |>
  complete(Transect, nesting(year, class), fill = list(tree.coh = 0)) |>
  filter(class == 1) |>  # 1 is the class for tree cover
  dplyr::select(-class)

# Join tree cover variables ----
covs_tree = tree_avg |>
  left_join(tree_coh, by = c("Transect", "year"))

# end script