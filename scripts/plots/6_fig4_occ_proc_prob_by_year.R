# Purpose: Create figure 4 from the dynamic occupancy model analysis for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
pacman::p_load(
  here,
  tidyverse,
  runjags
)

# Source plotting function ----
source(here("scripts", "utils_plot_occ_prob_by_year.R"))


# Import data and model output ----
load(here("data", "processed", "model_data.RData"))
mod_out_occ <- readRDS(here("output", "mod_out_occ.rds"))
mod_out_occ_df <- as.data.frame(t(as.matrix(mod_out_occ$mcmc)))

# Identify start and end years (for filtering and plotting)
yr_first <- as.numeric(min(easo_obs$year))
yr_last <- as.numeric(max(easo_obs$year))
yr_first_min_1 <- min(easo_obs$year) - 1
yr_last_min_1 <- max(easo_obs$year) - 1


# Format occupancy process parameters for plotting ----
## Calculate posterior means and 80 and 95% credible intervals ----
### Occupancy processes ----
occ_params_df <- mod_out_occ_df |>
  rownames_to_column() |>
  rename(mon_param = rowname) |>

  # Filter out detection parameters
  filter(!str_detect(mon_param, "p\\[")) |>

  # Extract site, year, and parameter info from the model output for monitored parameters
  mutate(
    site_id = if_else(
      str_detect(mon_param, "\\[\\d+,\\d+\\]"),
      as.numeric(str_extract(mon_param, "(?<=\\[)\\d+")),
      NA_real_
    ),
    year_id = if_else(
      str_detect(mon_param, "\\[\\d+,\\d+\\]"),
      as.numeric(str_extract(mon_param, "(?<=,)\\d+")),
      as.numeric(str_extract(mon_param, "\\d+"))
    ),
    param = str_extract(mon_param, "^[^\\[]*")) |>
  relocate(site_id, year_id, param, .before = V1) |>

  dplyr::select(param, site_id, year_id, everything(), -mon_param) |>

  # Convert dataframe from wide to long format
  pivot_longer(cols = contains("V"), names_to = "n_iter", values_to = "samp") |>

  # Calculate means and credible intervals by year (across sites)
  group_by(param, year_id) |>
  summarise(
    mean = mean(samp),
    lci.95 = quantile(samp, 0.025),
    uci.95 = quantile(samp, 0.975),
    lci.80 = quantile(samp, 0.10),
    uci.80 = quantile(samp, 0.90)
  ) |>

  # Calculate year column in YYYY format
  mutate(
    year = year_id + yr_first_min_1
  )

### Detection ----
p_psi_samples <- as.matrix(mod_out_occ$mcmc)
p_samples_raw <- as.data.frame(t(p_psi_samples[, grep("p[", colnames(p_psi_samples), fixed = T)]))

# Calculate survey-level mean detection
p_samples <- p_samples_raw |>

  # Calculate survey-level means
  mutate(p_surv_mean = rowMeans(across(everything()))) |>

  # Extract site, year, and parameter info from the model output for monitored parameters
  rownames_to_column() |>
  rename(mon_param = rowname) |>

  mutate(
    survey_id = if_else(
      str_detect(mon_param, "\\[\\d+,\\d+\\]"),
      as.numeric(str_extract(mon_param, "(?<=,)\\d+")),
      as.numeric(str_extract(mon_param, "\\d+"))
    ),
    param = str_extract(mon_param, "^[^\\[]*")) |>
  relocate(survey_id, param, .before = V1) |>
  dplyr::select(-mon_param)

# Calculate means and credible intervals by year (across surveys)
p_site_year <- easo_obs |>
  mutate(survey_id = row_number()) |>
  dplyr::select(Transect, year, survey_id) |>
  left_join(p_samples, by = "survey_id") |>
  dplyr::select(-survey_id) |>
  pivot_longer(cols = starts_with("V"), names_to = "p_samples_name", values_to = "p_samples") |>
  dplyr::select(-p_samples_name) |>
  group_by(param, year) |>
  summarise(
    mean = mean(p_samples),
    lci.95 = quantile(p_samples, 0.025),
    uci.95 = quantile(p_samples, 0.975),
    lci.80 = quantile(p_samples, 0.10),
    uci.80 = quantile(p_samples, 0.90)
  ) |>
  ungroup()


# Plot annual probabilities -----
(plot_phi <- plot.prob.annual.func(occ_params_df, "phi", "a) Persistence", "Probability", yr_first:yr_last_min_1))
(plot_col <- plot.prob.annual.func(occ_params_df, "gamma", "b) Colonization", "", yr_first:yr_last_min_1))
(plot_psi <- plot.prob.annual.func(occ_params_df, "phi", "c) Occupancy", "", yr_first:yr_last_min_1))
(plot_p <- plot.prob.annual.func(p_site_year, "p", "d) Detection", "", yr_first:yr_last))

(plots_occ_det = cowplot::plot_grid(
  plot_phi, plot_col, plot_psi, plot_p,
  nrow = 1,
  ncol = 4,
  align = "h"))

## Save plot
# ggsave(filename = here("output", "plots", "fig4_occ_proc_prob.png"),
#        plot = plots_occ_det,
#        width = 24,
#        height = 5,
#        units = "in",
#        dpi = 600)

# end script
