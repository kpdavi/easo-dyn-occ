# Purpose: Create figure 3 from the dynamic occupancy model analysis for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
pacman::p_load(
  here,
  tidyverse,
  runjags,
  cowplot
)

# Create / source functions ----
# Calculate the proportion of the posterior in the same direction as the mean
proportion_same_direction <- function(samples) {
  mean_posterior <- mean(samples)
  direction <- sign(mean_posterior)
  mean(sign(samples) == direction)
}

# Predict data from model
source(here("scripts", "utils", "utils_predict_data_from_model.R"))

# Plot predicted data
source(here("scripts", "utils", "utils_plot_predicted_data.R"))


# Import data and model output ----
## Data
surveys_cov_det <- read.csv(here("data", "processed", "surveys_covs_det.csv"), header = TRUE)
cov_tree <- read.csv(here("data", "processed", "covs_tree.csv"), header = TRUE)
cov_climate <- read.csv(here("data", "processed", "covs_climate.csv"), header = TRUE)
cov_lidar <- read.csv(here("data", "processed", "covs_lidar.csv"), header = TRUE)
coeff_names <- read.csv(here("data", "processed", "coeff_names.csv"), header = TRUE)

## Model output
mod_out <- readRDS(here("output", "mod_out_coeffs.rds"))


# Identify occupancy process covariates with probability of direction (pd) > 0.8 ----
mod_out_mat <- as.matrix(mod_out$mcmc)
mod_out_mat_colnames <- colnames(mod_out_mat)
covs_foc <- grep("^[PGH]", mod_out_mat_colnames)
covs_foc <- mod_out_mat_colnames[grep("^[PGH]", mod_out_mat_colnames)]

posterior_samples_foc_list <- lapply(covs_foc, function(col) mod_out_mat[, col, drop = FALSE])
names(posterior_samples_foc_list) <- covs_foc

## Create a data frame with the probability of direction
prop_dir_mn <- sapply(posterior_samples_foc_list, proportion_same_direction)
prop_dir_mn_df <- data.frame(param = names(prop_dir_mn),
                             pd = round(prop_dir_mn, 3))
prop_dir_mn_filt <- prop_dir_mn_df |>
  filter(pd > 0.8) |>
  left_join(coeff_names, by = "param")


# Predict and plot new data from the model ----
## Persistence ----
### Cohesion of tree cover ----
predict.cov.data(mod_out, cov_tree, "tree.coh", "H_occ", 2, "Cohesion of tree cover", "Persistence", "pred_phi_tree.coh", "pred_phi_tree.coh_sum")
(plot_phi_tree.coh <- plot.cov.preds(pred_phi_tree.coh_sum, "Cohesion of tree cover (%)\n", "Persistence", "Phi"))

### Breeding season minimum temperature ----
predict.cov.data(mod_out, cov_climate, "avg.breeding.tmin", "H_occ", 6, "Breeding season minimum temperature", "Persistence", "pred_phi_avg.breeding.tmin", "pred_phi_avg.breeding.tmin_sum")
(plot_phi_avg.breeding.tmin <- plot.cov.preds(pred_phi_avg.breeding.tmin_sum, "Breeding season\nminimum temperature (°C)", "Persistence", "Phi"))

### Winter minimum temperature ----
predict.cov.data(mod_out, cov_climate, "avg.winter.tmin", "H_occ", 4, "Winter minimum temperature", "Persistence", "pred_phi_avg.winter.tmin", "pred_phi_avg.winter.tmin_sum")
(plot_phi_avg.winter.tmin <- plot.cov.preds(pred_phi_avg.winter.tmin_sum, "Winter minimum\ntemperature (°C)", "Persistence", "Phi"))

### Winter minimum temperature lag 1 year ----
predict.cov.data(mod_out, cov_climate, "lag.avg.winter.tmin", "H_occ", 5, "Winter minimum temperature lag 1 year", "Persistence", "pred_phi_lag.avg.winter.tmin", "pred_phi_lag.avg.winter.tmin_sum")
(plot_phi_lag.avg.winter.tmin <- plot.cov.preds(pred_phi_lag.avg.winter.tmin_sum, "Winter minimum\ntemperature (°C) - lag 1 year", "Persistence", "Phi"))


## Colonization ----
### Cohesion of tree cover ----
predict.cov.data(mod_out, cov_tree, "tree.coh", "G_occ", 2, "Cohesion of tree cover", "Colonization", "pred_gamma_tree.coh", "pred_gamma_tree.coh_sum")
(plot_gamma_tree.coh <- plot.cov.preds(pred_gamma_tree.coh_sum, "Cohesion of tree cover (%)\n", "Colonization", "gamma"))

### Breeding season minimum temperature ----
predict.cov.data(mod_out, cov_climate, "avg.breeding.tmin", "G_occ", 6, "Breeding season minimum temperature", "Colonization", "pred_gamma_avg.breeding.tmin", "pred_gamma_avg.breeding.tmin_sum")
(plot_gamma_avg.breeding.tmin <- plot.cov.preds(pred_gamma_avg.breeding.tmin_sum, "Breeding season\nminimum temperature (°C)", "Colonization", "gamma"))

### Winter precipitation ----
predict.cov.data(mod_out, cov_climate, "sum.winter.prcp", "G_occ", 3, "Winter precipitation", "Colonization", "pred_gamma_sum.winter.prcp", "pred_gamma_sum.winter.prcp_sum")
(plot_gamma_sum.winter.prcp <- plot.cov.preds(pred_gamma_sum.winter.prcp_sum, "Winter precipitation (mm)\n", "Colonization", "gamma"))


## Initial occupancy ----
### Proportion tree cover ----
predict.cov.data(mod_out, cov_tree, "tree.prop", "P_occ", 2, "Proportion of tree cover", "Initial occupancy", "pred_psi1_tree.prop", "pred_psi1_tree.prop_sum")
(plot_psi1_tree.prop <- plot.cov.preds(pred_psi1_tree.prop_sum, "Tree cover (%)", "Initial occupancy", "psi", add_subscript = TRUE))


# Create combined plot ----
empty_grob <- grid::nullGrob()

## Plot
plot_preds_comb <- cowplot::plot_grid(

  # Persistence
  plot_phi_tree.coh, plot_phi_avg.breeding.tmin, plot_phi_avg.winter.tmin, plot_phi_lag.avg.winter.tmin,

  # Colonization
  plot_gamma_tree.coh, plot_gamma_avg.breeding.tmin, plot_gamma_sum.winter.prcp, empty_grob,

  # Initial occupancy
  plot_psi1_tree.prop, empty_grob, empty_grob, empty_grob,

  nrow = 3,
  ncol = 4,
  align = "h",
  axis = "l")

## Save plot
# ggsave(filename = here("output", "plots", "fig3_coeff_param_preds.png"),
#        plot = plot_preds_comb,
#        width = 18,
#        height = 11,
#        units = "in",
#        dpi = 600)

# end script
