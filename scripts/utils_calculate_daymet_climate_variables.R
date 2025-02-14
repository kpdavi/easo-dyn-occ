# Purpose: Function for calculating climate variables from DAYMET data; used for calculating climate covariates for the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

daymet.vars.calc = function(daymet_data_list, daymet_char, dates_vec, fun_name, daymet_var_foc) {
  
  # Select focal climate variable and subset focal dates
  daymet_stack <- daymet_data_list[[daymet_char]]
  dates_foc <- terra::subset(daymet_stack, subset = grep(paste(dates_vec, collapse = '|'), names(daymet_stack)))

  # Create index for year
  index_yr <- as.integer(substr(names(dates_foc), 1, 4))
  
  # Calculate focal metrics (e.g., average annual daily temperature: average daily temperature then averaged by year)
  daymet_stat_daily <- tapp(dates_foc, index_yr, fun_name, na.rm = TRUE)
  daymet_stats <- data.frame(daymet_stat_mean = round(global(daymet_stat_daily, "mean"), 3))
  names(daymet_stats) = daymet_var_foc
  
  # Format metrics dataframe
  daymet_stats_format <- tibble::rownames_to_column(daymet_stats, "year")
  
  return(daymet_stats_format)
}

# end script
