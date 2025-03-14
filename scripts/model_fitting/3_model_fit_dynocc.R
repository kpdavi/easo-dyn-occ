# Purpose: Fit and summarize output / calculate model diagnostics from the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Note: JAGS must be installed on one's computer to fit the model. Download JAGS at https://sourceforge.net/projects/mcmc-jags/

# Required packages ----
pacman::p_load(
  here,
  runjags,
  tidyverse,
  pROC
)

set.seed = 4242

# Import data ----
load(here("data", "processed", "model_data.RData"))


# Prepare data for the model ----
# Data list
data_list <- list(
  y = easo_obs$presence,
  site = easo_obs$site_id,
  year = easo_obs$year_id,
  nsites = max(easo_obs$site_id),
  nyears = max(easo_obs$year_id),
  nsurveys = as.numeric(nrow(easo_obs)),
  n_det_covs = cov_det_num,
  n_init_occ_covs = cov_psi1_num,
  n_col_covs = cov_gamma_num,
  n_pers_covs = cov_phi_num,
  X_det = X_det,
  X_psi1 = X_psi1,
  X_phi = X_phi,
  X_gam = X_gam,
  position.lidar = position_lidar
)

# Generate initial values
my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(1,
                 ncol = data_list$nyears,
                 nrow = data_list$nsites),
      P_occ = rnorm(cov_psi1_num, 0, 0.1),
      G_occ = rnorm(cov_gamma_num, 0, 0.1),
      H_occ = rnorm(cov_phi_num, 0, 0.1),
      A_det = rnorm(cov_det_num, 0, 0.1),
      mu.lidar = rnorm(1, 0, 0.1),
      sd.lidar = rgamma(1, 0.1, 0.1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# Specify parameters monitored ----
## Model output files can be large if all parameters of interest are monitored in a single model fit, so the code below fits the model with different groups of monitored parameters

## Coefficients, parameters for missing LiDAR data, and Bayesian p-value for log likelihood (scalars)
monitor_params <- c("P_occ", "G_occ", "H_occ", "A_det", "mu.lidar", "sd.lidar", "log.bpv")

## Occupancy process parameters
# monitor_params <- c("phi", "gamma", "psi", "p")

## Parameters for calculating area under the receiver operating characteristic curve (AUC)
monitor_params <- c("pocc.sim")


# Fit the model ----
mod_out <- run.jags(
  model = here("scripts", "model_fitting", "model_code_dynocc.txt"),
  monitor = monitor_params,
  data = data_list,
  inits = my_inits,
  n.chains = 4,
  adapt = 1000,
  burnin = 1000,
  sample = 10000,
  thin = 5,
  method = "parallel"
)

# Save output ----
# Save output based on whether coefficients / vectors, parameters for calculating AUC, or occupancy process parameters are being monitored
if (length(monitor_params) > 4) {
  # Coefficients / scalars
  saveRDS(mod_out,
          file = here("output", "mod_out_coeffs.rds"))
} else if (length(monitor_params) == 4) {
  # Occupancy process parameters
  saveRDS(mod_out,
          file = here("output", "mod_out_occ.rds"))
} else {
  # Parameters for calculating AUC
  saveRDS(mod_out,
          file = here("output", "mod_out_auc.rds"))
}


# Model diagnostics ----
## Bayesian p value for log likelihood
mod_out_occ <- readRDS(here("output", "mod_out_occ.rds"))
mod_out_occ_jags_sum <- add.summary(mod_out_occ, confidence = c(0.95))
mod_out_occ <- data.frame(mod_out_occ_jags_sum$summaries)

(bpv_loglik <- mod_out_occ |>
  rownames_to_column() |>
  rename(param = rowname) |>
  filter(param == "log.bpv") |>
  pull(Mean))

## Calculate area under the receiver operator curve (AUC) ----
## Calculate AUC for the prevalence process (i.e., at the level of the survey)
mod_out_auc <- readRDS(here("output", "mod_out_auc.rds"))
pocc_samples <- as.matrix(mod_out_auc$mcmc)
pocc_mn <- data.frame(pocc = colnames(pocc_samples),
                      mn = apply(pocc_samples, 2, mean))
y_obs <- easo_obs$presence
(auc.prev.pocc <- round(as.numeric(roc(y_obs, pocc_mn$mn)$auc), 3)) # compares the observed data to data generated from the model

# end script
