# Purpose: Compile detection covariate data for the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
pacman::p_load(
  here,
  tidyverse,
  suncalc
)


# Import data ----
## Surveys ----
surveys <- read.csv(here("data", "raw", "surveys.csv"), header = TRUE)


# Generate detection covariates ----
## Ordinal date ----
surveys$date.ord <- yday(surveys$date)

## Moon phase ----
# Create a combined date / time column (in POSIXct format)
date_time <- surveys |>
  dplyr::select(Transect, year, date, start.time) |>
  mutate(start.time.form = as.POSIXct(paste(date, start.time), format = "%Y-%m-%d %H:%M"))

# Calculate moon phase
moon_frac <- getMoonIllumination(date = date_time$start.time.form, keep = "fraction")

# Reformat moon phase data
moon_frac_form <- moon_frac |>
  rename(moon.phase = fraction,
         date.start.time = date) |>
  mutate(moon.phase = round(moon.phase, 3),
         date = ymd(as.Date(date.start.time, tz = "America/Denver")),
         start.time = format(as.POSIXct(date.start.time), format = "%H:%M")) |>
  dplyr::select(date, start.time, moon.phase)

moon_phase <- date_time |>
  dplyr::select(Transect, year) |>
  bind_cols(moon_frac_form)


# Join detection covariates to survey data ----
surveys_covs_det <- surveys |>
  mutate(date = ymd(date)) |>
  left_join(moon_phase, by = c("Transect", "year", "date", "start.time"))

# write.csv(surveys_covs_det, here("data", "processed", "surveys_covs_det.csv"), row.names = FALSE)

# end script
