# Purpose: Summarize output / calculate model diagnostics from the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
pacman::p_load(
  here,
  runjags,
  tidyverse,
  pROC
)

# Functions ----
## Calculate the proportion of the posterior in the same direction as the mean
proportion_same_direction <- function(samples) {
  mean_posterior <- mean(samples)
  direction <- sign(mean_posterior)
  mean(sign(samples) == direction)
}

# Import data / model outputs ----
surveys <- read.csv(here("data", "raw", "surveys.csv"))
coeff_names <- read.csv(here("data", "processed", "coeff_names.csv"), header = TRUE)
mod_out_occ_jags <- readRDS(here("output", "mod_out_occ.rds"))
mod_out_auc_jags <- readRDS(here("output", "mod_out_auc.rds"))
mod_out_coeffs_jags <- readRDS(here("output", "mod_out_coeffs.rds"))

## Format names of coefficients
det_covs <- data.frame(
  proc = rep("Detection", 3),
  param = paste0("A_det[", 1:3, "]"),
  names_real = c("Intercept", "Ordinal date", "Moon phase")
)

coeff_names <- coeff_names |>
  dplyr::select(param, proc, names_real) |>
  bind_rows(det_covs) |>
  rename(
    Parameter = param,
    Covariate = names_real,
    Process = proc
  )

# Model diagnostics ----
## Bayesian p value for log likelihood ----
mod_out_coeff_jags_sum <- add.summary(mod_out_coeffs_jags, confidence = c(0.95))
mod_out_coeff <- data.frame(mod_out_coeff_jags_sum$summaries)

(bpv_loglik <- mod_out_coeff |>
    rownames_to_column() |>
    rename(param = rowname) |>
    filter(param == "log.bpv") |>
    pull(Mean))

## Calculate area under the receiver operator curve (AUC) ----
## Calculate AUC for the prevalence process (i.e., at the level of the survey). Compares the observed data to data generated from the model.
pocc_samples <- as.matrix(mod_out_auc_jags$mcmc)
pocc_mn <- data.frame(pocc = colnames(pocc_samples),
                      mn = apply(pocc_samples, 2, mean))
y_obs <- surveys$presence
(auc_prev_pocc <- round(as.numeric(roc(y_obs, pocc_mn$mn)$auc), 3))

# end script
