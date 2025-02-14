# Purpose: Compile climate data to include as covariates in the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Notes ----
## Anonymity of survey locations ----
# The locations of the monitoring sites are not public, so the longitude and latitude coordinates in the sites.csv were randomly generated. However, they were generated using the minimum and maximum longitude / latitude values present in the real data, so the climate variables calculated in this script are those that were used as covariates in the model.
# Output file for this script: processed/covs_climate.csv. The data to run the analyses for this manuscript is in Zenodo.

## Data downloads ----
# DAYMET data were downloaded and compiled in April 2022.


# Required packages and functions ----
# install.packages("pacman")
pacman::p_load(
  tidyverse,
  sf,
  terra,
  ncdf4,
  FedData
)

# Function for calculating climate variables
source("/Users/kristindavis/Library/CloudStorage/OneDrive-NewMexicoStateUniversity/PhD/Code/R code/EASO/Data compilation/utils_func_calculate_climate_variables.R")


# Global variables ----
root_dir <- "/Users/kristindavis/Library/CloudStorage/OneDrive-NewMexicoStateUniversity/PhD/Data/EASO/publication_data/"
clim_vars <- c("prcp", "tmin")  # names of the focal climate variables in DAYMET
year_vec <- c(2013:2022)  # focal years


# Import data ----
## Transects ----
sites <- read.csv(paste0(root_dir, "raw/sites.csv"), header = TRUE)
sites_sp <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)

## Climate ----
### Create bounding box ----
## Used as the extent within which focal climate variables are summarized
### Identify min and max values for longitude / latitude for determining the bounding box edges
sites_coords_df <- as.data.frame(st_coordinates(sites_sp))

### Create raster of the bounding box
bbox_fc <- rast()
ext(bbox_fc) <- c(min(sites_coords_df$X), max(sites_coords_df$X), min(sites_coords_df$Y), max(sites_coords_df$Y))

### Download DAYMET data ----
## Note: it takes some time for the FedData::get_daymet() function to run
fc_daymet_select <- FedData::get_daymet(
  # Supply the focal area
  template = bbox_fc,
  #--- label ---#
  label = "FC_select",
  #--- variables to download ---#
  elements = clim_vars,
  #--- years ---#
  years = year_vec
)

# Calculate climate variables ----
## Specify date range for breeding season and winter ----
### Breeding season: March 1-June 15 (Birds of the World species account) 
### Winter: November 1-February 28/29 preceding the breeding season
### Picked a year at random in order to correctly specify the date range
breeding_start <- "2022-03-01"
breeding_end <- "2022-06-15"
breeding_dates_full <- make.names(format(seq(ymd(breeding_start), ymd(breeding_end), by = 'days'), "%m-%d")) # make.names() changes the hyphens in the vector to periods, which matches the format of the layer names in the DAYMET data
breeding_dates <- gsub("X", ".", paste(breeding_dates_full, collapse = '|'))

winter_start <- "2019-11-01"
winter_end <- "2020-02-29"
winter_dates_full <- make.names(format(seq(ymd(winter_start), ymd(winter_end), by = 'days'), "%m-%d"))
winter_dates <- gsub("X", ".", paste(winter_dates_full, collapse = '|'))

## Apply calculation function ----
## Breeding season average minimum temperature
tmin_breeding_avg <- daymet.vars.calc(fc_daymet_select, "tmin", breeding_dates, mean, "tmin.breeding.avg")
tmin_breeding_avg <- transform(tmin_breeding_avg, tmin.breeding.avg.lag1 = c(tmin.breeding.avg[-1],NA)) # 1-year-lagged variable (i.e., t - 1)

## Winter average minimum temperature
tmin_winter_avg <- daymet.vars.calc(fc_daymet_select, "tmin", winter_dates, mean, "tmin.winter.avg")
tmin_winter_avg <- transform(tmin_winter_avg, tmin.winter.avg.lag1 = c(tmin.winter.avg[-1],NA)) # 1-year-lagged variable (i.e., t - 1)

## Winter cumulative precipitation
prcp_winter_sum <- daymet.vars.calc(fc_daymet_select, "prcp", winter_dates, sum, "winter.prcp.sum")


# Join climate variables ----
clim_calcs_comb <- plyr::join_all(list(tmin_breeding_avg, tmin_winter_avg, prcp_winter_sum), by = "year", type = "left")

covs_climate <- clim_calcs_comb |>
  mutate(year = as.numeric(gsub("X", "" , year))) |>
  filter(year != max(year)) 

# end script