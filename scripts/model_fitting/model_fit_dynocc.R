# Purpose: Fit the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Note: JAGS must be installed on one's computer to fit the model. Download JAGS at https://sourceforge.net/projects/mcmc-jags/

# Required packages ----
# install.packages("pacman")
pacman::p_load(
  runjags
)

# library(car)
# library(MCMCvis)
# library(pROC)

set.seed = 4242

# Import data ----
load("data/processed/model_data.RData")

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
  position_lidar = position_lidar
)

# Initial values
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


# Fit the model ----
monitor_params <- c("P_occ", "G_occ", "H_occ", "A_det", "mu.lidar", "sd.lidar", "log.bpv")
# monitor_params <- c("pocc.sim", "psi")
# monitor_params <- c("phi", "gamma", "psi", "p")

my_mod <- run.jags(
  model = 'scripts/model_fitting/model_code_dynocc.txt',
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

if (length(monitor_params) > 4) {
  saveRDS(my_mod, 
          file = paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, ".rds"))
} else if (length(monitor_params) == 4) {
  saveRDS(my_mod, 
          file = paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_occ_process_params_saved.rds"))
} else {
  saveRDS(my_mod, 
          file = paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_auc_params_saved.rds"))
}
        
# View and save results ----
plot(my_mod, vars = c("P_occ", 
                      "G_occ", 
                      "H_occ", 
                      "A_det", 
                      "mu.lidar", 
                      "sd.lidar"),
     file = paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, ".pdf"))

# Summary
sims_obj <- add.summary(my_mod, confidence = c(0.95))
sim_sum <- data.frame(sims_obj$summaries)

# Explore the psrf and effective sample size
sim_sum %>% filter(psrf > 1.1)
sim_sum %>% filter(SSeff < 1000)

# Save output
sim_sum <- round(sim_sum, 3)
write.csv(sim_sum, paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, ".sim.sum.csv"))

sim_samples <- as.matrix(my_mod$mcmc)
write.csv(sim_samples, paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, ".sim.samples.csv"))

sink(paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, "_model_results.txt"))
print(my_mod)
sink()

sink(paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, "_model_summary.txt"))
summary(my_mod)
sink()


# Model checking ----
## AUC ----
# AUC for prevalance process
## Prevalence is at the level of a survey
my_mod <- readRDS(paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_auc_params_saved.rds"))
sim_samples <- as.matrix(my_mod$mcmc)
pocc_samples <- sim_samples[, grepl("pocc.sim", colnames(sim_samples))]
pocc_mn <- data.frame(pocc = colnames(pocc_samples),
                      mn = apply(pocc_samples, 2, mean))
y_obs <- easo_obs$presence
(auc.prev.pocc <- round(as.numeric(roc(y_obs,pocc_mn$mn)$auc), 3)) # compares the observed data to detection probability times occupancy probability

sink(paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_auc_prev.txt"))
print(auc.prev.pocc)
sink()

# AUC for occupancy process
## Occupancy is at the level of site and year
psi_samples <- t(sim_samples[, grepl("psi", colnames(sim_samples))])
psi_mn <- apply(psi_samples, 1 , mean)
psi_mn_vec <- numeric(nrow(easo_surv_sites))
for (i in 1:nrow(easo_surv_sites)) {
  psi_mn_vec[i] <- psi_mn[easo_surv_sites$site_surv_id[i]]
}
y_site_yr_obs <- easo_surv_sites$site_yr_pres
(auc.occ <- round(suppressMessages(as.numeric(roc(y_site_yr_obs, psi_mn_vec)$auc)), 3)) # <- this tells us how well the model is figuring out where EASO is on the landscape. Gives a 0 - 1 rate of success.

sink(paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_auc_occ.txt"))
print(auc.occ)
sink()


# Explore model output ----
my_mod <- readRDS(file = paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_out1_", date, ".rds"))
mymodel_mat = as.matrix(my_mod$mcmc)

## Tables ----
col_names <- colnames(mymodel_mat)
foc_params <- grep("^[PGHA]", col_names)
foc_params <- col_names[grep("^[PGHA]", col_names)]

posterior_samples_list <- lapply(foc_params, function(col) mymodel_mat[, col, drop = FALSE])
names(posterior_samples_list) <- foc_params

# Apply the function to each focal parameter and create a data frame
prop_dir_mn <- sapply(posterior_samples_list, proportion_same_direction)
prop_dir_mn_df <- data.frame(Parameter = names(prop_dir_mn), 
                             ProportionSameDirection = round(prop_dir_mn, 3))
write.csv(prop_dir_mn_df, paste0(root, model_structure,"/", subf, "/easo_", subf_ch, "_prop_dir_mn_", date, ".csv"), row.names = FALSE)


## Plots ----
### Caterpillar plot ----
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_cat_plot_", date, ".png"), 
# png(paste0(root, model_structure,"/", subf, "/plots/easo_covs_cat_plot_", date, ".png"), 
    width = 11, 
    height = 7, 
    units = "in", 
    res = 600)
(MCMCplot(mymodel_mat, 
         params = c("P_occ", 
                    "G_occ", 
                    "H_occ", 
                    "A_det"), 
         ci = c(80, 95),
         main = "EASO dynamic occupancy model",
         labels = c(expression(psi[1] * ": intercept"),
                    expression(psi[1] * ": proportion tree cover (2013)"),
                    expression(psi[1] * ": SD in veg height (2013)"),
                    expression(gamma * ": intercept"),
                    expression(gamma * ": cohesion"),
                    expression(gamma * ": winter precip (mm)"),
                    expression(gamma * ": average winter minimum temperature"),
                    expression(gamma * ": average winter minimum temperature, 1 yr lag"),
                    expression(gamma * ": average breeding season minimum temperature"),
                    expression(phi * ": intercept"),
                    expression(phi * ": cohesion"),
                    expression(phi * ": winter precip (mm)"),
                    expression(phi * ": average winter minimum temperature"),
                    expression(phi * ": average winter minimum temperature, 1 yr lag"),
                    expression(phi * ": average breeding season minimum temperature"),
                    expression(phi * ": average breeding season minimum temperature, 1 yr lag"),
                    expression(p * ": intercept"),
                    expression(p * ": moon phase"),
                    expression(p * ": ordinal date")),
         ref_ovl = TRUE))
dev.off()

### Caterpillar plot - faceted ----
mod_param_ests = dplyr::select(as.data.frame(mymodel_mat), matches("^G_|^P_|^H_"))
mod_coef_vect <- apply(mod_param_ests, 2, mean)
mod_lci_80 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.1)))
mod_uci_80 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.9)))
mod_lci_95 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.025)))
mod_uci_95 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.975)))

### Replacement vector to change column names to covariates
mod_param_replacement_vector <- c(
  `P_occ[1]` = "Initial occupancy-Intercept",
  `P_occ[2]` = "Initial occupancy-Proportion of tree cover",
  `P_occ[3]` = "Initial occupancy-Standard deviation in vegetation height (m)",
  `G_occ[1]` = "Colonization-Intercept",
  `G_occ[2]` = "Colonization-Cohesion of tree cover (%)",
  `G_occ[3]` = "Colonization-Cumulative winter precipitation (mm)",
  `G_occ[4]` = "Colonization-Minimum winter temperature (\u00B0C)",
  `G_occ[5]` = "Colonization-Minimum winter temperature lag 1 year (\u00B0C)",
  `G_occ[6]` = "Colonization-Minimum breeding season temperature (\u00B0C)",
  `H_occ[1]` = "Persistence-Intercept",
  `H_occ[2]` = "Persistence-Cohesion of tree cover (%)",
  `H_occ[3]` = "Persistence-Cumulative winter precipitation (mm)",
  `H_occ[4]` = "Persistence-Minimum winter temperature (\u00B0C)",
  `H_occ[5]` = "Persistence-Minimum winter temperature lag 1 year (\u00B0C)",
  `H_occ[6]` = "Persistence-Minimum breeding season temperature (\u00B0C)",
  `H_occ[7]` = "Persistence-Minimum breeding season temperature lag 1 year (\u00B0C)"
)
mod_param_ests_renamed <- mod_param_ests |>
  rename_with(~ mod_param_replacement_vector)

### Format means and estimates data for plotting
mod_cat_dat <- data.frame(coeff.names = colnames(mod_param_ests_renamed),
                          mn = mod_coef_vect, 
                          lci.80 = mod_lci_80,
                          uci.80 = mod_uci_80,
                          lci.95 = mod_lci_95,
                          uci.95 = mod_uci_95)
mod_cat_dat <- mod_cat_dat |>
  rownames_to_column() |>
  rename(model_param = rowname) |>
  separate(coeff.names, into = c("proc_name", "coeff"), sep = "-") |>
  mutate(proc_name = fct_relevel(proc_name, c("Persistence", "Colonization", "Initial occupancy"))) |>
  group_by(proc_name) |>
  mutate(coeff = factor(coeff, levels = c("Intercept", setdiff(unique(coeff), "Intercept"))),
         coeff = fct_rev(coeff))

### Caterpillar plot (manual)
cat_plot_labels <- as.vector(c("a) Persistence", "b) Colonization", "c) Initial occupancy"))
names(cat_plot_labels) <- c("Persistence", "Colonization", "Initial occupancy")

(man_cat_plot_occ <- ggplot(mod_cat_dat, aes(x = mn, y = coeff)) + 
    geom_point(size = 3) +
    geom_segment(aes(x = lci.95, xend = uci.95, y = coeff, yend = coeff), color = "gray50", linewidth = 0.75) +
    geom_segment(aes(x = lci.80, xend = uci.80, y = coeff, yend = coeff), color = "black", linewidth = 1.5) +
    geom_vline(xintercept = 0, colour = "blue", linetype = 2) + 
    facet_wrap(vars(proc_name), labeller = labeller(proc_name = cat_plot_labels)) +
    labs(x = "Parameter estimate", 
         y = "") +
    theme_bw() +
    theme(axis.text = element_text(size = 20, color = "black"),
          axis.title = element_text(size = 22, color = "black"),
          strip.text = element_text(size = 20, face = "bold", color = "black", hjust = 0),
          strip.background = element_blank())
)
ggsave(filename = paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_caterpillar_plot_faceted_", date, ".png"),
       plot = man_cat_plot_occ,
       width = 23, 
       height = 6.5, 
       units = "in",
       dpi = 600)


### Caterpillar plot - detection ----
mod_param_ests_det = dplyr::select(as.data.frame(mymodel_mat), matches("^A_"))
mod_coef_vect_det <- apply(mod_param_ests_det, 2, mean)
mod_lci_80_det <- apply(mod_param_ests_det, 2, function(x) quantile(x, probs = c(0.1)))
mod_uci_80_det <- apply(mod_param_ests_det, 2, function(x) quantile(x, probs = c(0.9)))
mod_lci_95_det <- apply(mod_param_ests_det, 2, function(x) quantile(x, probs = c(0.025)))
mod_uci_95_det <- apply(mod_param_ests_det, 2, function(x) quantile(x, probs = c(0.975)))

### Replacement vector to change column names to covariates
mod_param_replacement_vector_det <- c(
  `A_det[1]` = "Intercept",
  `A_det[2]` = "Moon phase (%)",
  `A_det[3]` = "Ordinal date"
)
mod_param_ests_renamed_det <- mod_param_ests_det |>
  rename_with(~ mod_param_replacement_vector_det)

### Format means and estimates data for plotting
mod_cat_dat_det <- data.frame(coeff.names = colnames(mod_param_ests_renamed_det),
                          mn = mod_coef_vect_det, 
                          lci.80 = mod_lci_80_det,
                          uci.80 = mod_uci_80_det,
                          lci.95 = mod_lci_95_det,
                          uci.95 = mod_uci_95_det)
mod_cat_dat_det <- mod_cat_dat_det |>
  mutate(coeff.names = factor(coeff.names, levels = c("Intercept", "Moon phase (%)", "Ordinal date")),
         coeff.names = fct_rev(coeff.names))

(man_cat_plot_det <- ggplot(mod_cat_dat_det, aes(x = mn, y = coeff.names)) + 
    geom_point(size = 3) +
    geom_segment(aes(x = lci.95, xend = uci.95, y = coeff.names, yend = coeff.names), color = "gray50", linewidth = 0.75) +
    geom_segment(aes(x = lci.80, xend = uci.80, y = coeff.names, yend = coeff.names), color = "black", linewidth = 1.5) +
    geom_vline(xintercept = 0, colour = "blue", linetype = 2) + 
    labs(x = "Parameter estimate", 
         y = "") +
    theme_bw() +
    theme(axis.text = element_text(size = 20, color = "black"),
          axis.title = element_text(size = 22, color = "black"))
)
ggsave(filename = paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_caterpillar_plot_detection_", date, ".png"),
       plot = man_cat_plot_det,
       width = 8, 
       height = 6, 
       units = "in",
       dpi = 600)

### Annual estimates of occupancy and detection processes ----
### Calculate mean and CRIs
my_mod_occ_params <- readRDS(paste0(root, model_structure,"/", subf, "/model diagnostics/easo_", subf_ch, "_out1_", date, "_occ_process_params_saved.rds"))
mymodel_occ_params_df <- as.data.frame(t(as.matrix(my_mod_occ_params$mcmc)))

#### Occupancy process
occ_params_df <- mymodel_occ_params_df |>
  rownames_to_column() |>
  rename(mon_param = rowname) |>
  mutate(
    site = if_else(
      str_detect(mon_param, "\\[\\d+,\\d+\\]"),
      as.numeric(str_extract(mon_param, "(?<=\\[)\\d+")),
      NA_real_
    ),
    year = if_else(
      str_detect(mon_param, "\\[\\d+,\\d+\\]"),
      as.numeric(str_extract(mon_param, "(?<=,)\\d+")),
      as.numeric(str_extract(mon_param, "\\d+"))
    ),
    param = str_extract(mon_param, "^[^\\[]*")) |>
  filter(param != "p") |>
  dplyr::select(param, site, year, everything(), -mon_param) |>
  pivot_longer(cols = contains("V"), names_to = "n_iter", values_to = "samp") |>
  group_by(param, year) |>
  summarise(
    mean = mean(samp),
    lci = quantile(samp, 0.025),
    uci = quantile(samp, 0.975)
  ) |>
  rename(year_num = year) |>
  mutate(
    year = year_num + 2012
  )

#### Detection
p_psi_samples <- as.matrix(my_mod_occ_params$mcmc)
p_samples_raw <- as.data.frame(t(p_psi_samples[, grep("p[", colnames(p_psi_samples), fixed = T)]))
p_samples <- p_samples_raw |>
  mutate(p_surv_mean = rowMeans(across(starts_with("p"))),
  )
p_site_year = easo_df_thres |>
  dplyr::select(Transect, year) |>
  bind_cols(p_samples) |>
  pivot_longer(cols = starts_with("V"), names_to = "p_samples_name", values_to = "p_samples") |>
  dplyr::select(-p_samples_name) |>
  group_by(year) |>
  summarise(
    mean = mean(p_samples),
    lci = quantile(p_samples, 0.025),
    uci = quantile(p_samples, 0.975)
  )
p_site_year <- as.data.frame(p_site_year)

## Plot
(col_plot <- occ.plot.func(occ_params_df, "gamma", "colonization", "γ", 2013:2020))
(phi_plot <- occ.plot.func(occ_params_df, "phi", "persistence", "ϕ", 2013:2020))
(psi_plot <- occ.plot.func(occ_params_df, "psi", "occupancy", "ψ", 2013:2021))

(p_plot <- ggplot(p_site_year, aes(year, mean)) +
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) + 
    labs(x = "Year",
         y = "Probability of detection (p)") +
    scale_x_continuous(breaks = 2013:2021) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
)

comb_occ_det_plots = cowplot::plot_grid(
  psi_plot, phi_plot, col_plot, p_plot,
  labels = c("a) Occupancy", "b) Persistence", "c) Colonization", "d) Detection"),
  nrow = 1, 
  ncol = 4,
  label_size = 18,
  align = "h")
ggsave(filename = paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_occ_process_det_plots_", date, ".png"),
       plot = comb_occ_det_plots,
       width = 24, 
       height = 5, 
       units = "in",
       dpi = 600)

### Posteriors plotted with priors ----
mymodel_mat = as.matrix(my_mod$mcmc)
gen_prior <- rlogis(nrow(mymodel_mat), 0, 1)

#### Detection
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_det_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("A_det[1]", "A_det[2]", "A_det[3]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

### Initial occupancy
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_init_occ_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("P_occ[1]", "P_occ[2]", "P_occ[3]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

# png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_init_occ_cov4_post_prior_", date, ".png"), 
#     width = 9, 
#     height = 7, 
#     units = "in", 
#     res = 600)
# MCMCtrace(mymodel_mat, 
#           params = c("P_occ[4]"),
#           ISB = FALSE,
#           exact = TRUE,
#           priors = gen_prior,
#           pdf = FALSE,
#           post_zm = FALSE)
# dev.off()

### Colonization
#### Covariates 1 - 3
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_col_covs1-3_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("G_occ[1]", "G_occ[2]", "G_occ[3]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

#### Covariates 4 - 6
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_col_covs4-6_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("G_occ[4]", "G_occ[5]", "G_occ[6]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

### Persistence
#### Covariates 1 - 3
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pers_covs1-3_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("H_occ[1]", "H_occ[2]", "H_occ[3]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

#### Covariates 4 - 6
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pers_covs4-6_post_prior_", date, ".png"), 
    width = 9, 
    height = 7, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("H_occ[4]", "H_occ[5]", "H_occ[6]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()

#### Covariate 7
png(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pers_covs7_post_prior_", date, ".png"), 
    width = 5, 
    height = 3, 
    units = "in", 
    res = 600)
MCMCtrace(mymodel_mat, 
          params = c("H_occ[7]"),
          ISB = FALSE,
          exact = TRUE,
          priors = gen_prior,
          pdf = FALSE,
          post_zm = FALSE)
dev.off()



### Predictions ----
#### Initial occupancy ----
##### Proportion tree cover ----
# Generate new data
sd.perc.cont <- seq(min(easo_df_thres_unsc$perc.cont, na.rm = TRUE), max(easo_df_thres_unsc$perc.cont, na.rm = TRUE), length.out = 50)
sd.perc.cont_mean <- mean(easo_df_thres_unsc$perc.cont, na.rm = TRUE); sd.perc.cont_sd <- sd(easo_df_thres_unsc$perc.cont, na.rm = TRUE)
sd.perc.cont_s <- (sd.perc.cont - sd.perc.cont_mean) / sd.perc.cont_sd

# Separate out focal mcmc chains
P_int = combine.mcmc(my_mod, vars = "P_occ[1]")
P_perc_cont = combine.mcmc(my_mod, vars = "P_occ[2]")

# Predict new data
psi1_pred <- matrix(NA, nrow = length(sd.perc.cont_s), ncol = length(P_perc_cont))
for (i in seq_along(sd.perc.cont_s)) {
  psi1_pred[i, ] <- 1/(1 + exp(-1*(P_int + P_perc_cont*sd.perc.cont_s[i])))
}

newdat.sd.perc.cont <- data.frame(psi1 = apply(psi1_pred, 1, mean), 
                                  psi1_lci_95 = apply(psi1_pred, 1, function(x) quantile(x, 0.025)),
                                  psi1_uci_95 = apply(psi1_pred, 1, function(x) quantile(x, 0.975)),
                                  psi1_lci_80 = apply(psi1_pred, 1, function(x) quantile(x, 0.1)),
                                  psi1_uci_80 = apply(psi1_pred, 1, function(x) quantile(x, 0.9)),
                                  sd.perc.cont = sd.perc.cont)

# Plot predictions
(sd.perc.cont.plot <- ggplot(newdat.sd.perc.cont, aes(x = sd.perc.cont, y = psi1)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    # scale_x_continuous(breaks = seq(0.5, 1.5, .25)) +
    labs(x = "Tree cover in 2013 (%)", 
         y = bquote(Probability ~of ~initial ~occupancy ~(psi[1]))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_psi1_perc_cont_", date, ".png"),
       sd.perc.cont.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)


#### Colonization ----
##### Cohesion of tree cover ----
# Generate new data
sd.coh <- seq(min(easo_df_thres_unsc$cohesion, na.rm = TRUE), max(easo_df_thres_unsc$cohesion, na.rm = TRUE), length.out = 50)
sd.coh_mean <- mean(easo_df_thres_unsc$cohesion, na.rm = TRUE); sd.coh_sd <- sd(easo_df_thres_unsc$cohesion, na.rm = TRUE)
sd.coh_s <- (sd.coh - sd.coh_mean) / sd.coh_sd

# Separate out focal mcmc chains
G_int = combine.mcmc(my_mod, vars = "G_occ[1]")
G_coh = combine.mcmc(my_mod, vars = "G_occ[2]")

# Predict new data
gamma_pred <- matrix(NA, nrow = length(sd.coh_s), ncol = length(G_coh))
for (i in seq_along(sd.coh_s)) {
  gamma_pred[i, ] <- 1/(1 + exp(-1*(G_int + G_coh*sd.coh_s[i])))
}

newdat.sd.coh <- data.frame(gamma = apply(gamma_pred, 1, mean), 
                                  psi1_lci_95 = apply(gamma_pred, 1, function(x) quantile(x, 0.025)),
                                  psi1_uci_95 = apply(gamma_pred, 1, function(x) quantile(x, 0.975)),
                                  psi1_lci_80 = apply(gamma_pred, 1, function(x) quantile(x, 0.1)),
                                  psi1_uci_80 = apply(gamma_pred, 1, function(x) quantile(x, 0.9)),
                                  sd.coh = sd.coh)

# Plot predictions
(sd.coh.plot <- ggplot(newdat.sd.coh, aes(x = sd.coh, y = gamma)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Cohesion of tree cover (%)", 
         y = bquote(Probability ~of ~colonization ~(gamma))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_gamma_coh_", date, ".png"),
       sd.coh.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)


##### Cumulative winter precipitation ----
# Generate new data
sd.wp <- seq(min(easo_df_thres_unsc$sum.winter.prcp, na.rm = TRUE), max(easo_df_thres_unsc$sum.winter.prcp, na.rm = TRUE), length.out = 50)
sd.wp_mean <- mean(easo_df_thres_unsc$sum.winter.prcp, na.rm = TRUE); sd.wp_sd <- sd(easo_df_thres_unsc$sum.winter.prcp, na.rm = TRUE)
sd.wp_s <- (sd.wp - sd.wp_mean) / sd.wp_sd

# Separate out focal mcmc chains
G_int = combine.mcmc(my_mod, vars = "G_occ[1]")
G_wp = combine.mcmc(my_mod, vars = "G_occ[3]")

# Predict new data
gamma_wp_pred <- matrix(NA, nrow = length(sd.wp_s), ncol = length(G_wp))
for (i in seq_along(sd.wp_s)) {
  gamma_wp_pred[i, ] <- 1/(1 + exp(-1*(G_int + G_wp*sd.wp_s[i])))
}

newdat.sd.wp <- data.frame(gamma = apply(gamma_wp_pred, 1, mean), 
                            psi1_lci_95 = apply(gamma_wp_pred, 1, function(x) quantile(x, 0.025)),
                            psi1_uci_95 = apply(gamma_wp_pred, 1, function(x) quantile(x, 0.975)),
                            psi1_lci_80 = apply(gamma_wp_pred, 1, function(x) quantile(x, 0.1)),
                            psi1_uci_80 = apply(gamma_wp_pred, 1, function(x) quantile(x, 0.9)),
                            sd.wp = sd.wp)

# Plot predictions
(sd.wp.plot <- ggplot(newdat.sd.wp, aes(x = sd.wp, y = gamma)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Cumulative winter precipitation (mm)", 
         y = bquote(Probability ~of ~colonization ~(gamma))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_gamma_wp_", date, ".png"),
       sd.wp.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)

##### Average winter minimum temperature - lag 1 year ----
# Generate new data
sd.wtl <- seq(min(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE), max(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE), length.out = 50)
sd.wtl_mean <- mean(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE); sd.wtl_sd <- sd(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE)
sd.wtl_s <- (sd.wtl - sd.wtl_mean) / sd.wtl_sd

# Separate out focal mcmc chains
G_int = combine.mcmc(my_mod, vars = "G_occ[1]")
G_wtl = combine.mcmc(my_mod, vars = "G_occ[5]")

# Predict new data
gamma_wtl_pred <- matrix(NA, nrow = length(sd.wtl_s), ncol = length(G_wtl))
for (i in seq_along(sd.wtl_s)) {
  gamma_wtl_pred[i, ] <- 1/(1 + exp(-1*(G_int + G_wtl*sd.wtl_s[i])))
}

newdat.sd.wtl <- data.frame(gamma = apply(gamma_wtl_pred, 1, mean), 
                            psi1_lci_95 = apply(gamma_wtl_pred, 1, function(x) quantile(x, 0.025)),
                            psi1_uci_95 = apply(gamma_wtl_pred, 1, function(x) quantile(x, 0.975)),
                            psi1_lci_80 = apply(gamma_wtl_pred, 1, function(x) quantile(x, 0.1)),
                            psi1_uci_80 = apply(gamma_wtl_pred, 1, function(x) quantile(x, 0.9)),
                            sd.wtl = sd.wtl)

# Plot predictions
(sd.gamma.wtl.plot <- ggplot(newdat.sd.wtl, aes(x = sd.wtl, y = gamma)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Average winter minimum temperature (°C) - lag 1 year", 
         y = bquote(Probability ~of ~colonization ~(gamma))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_gamma_wtl_", date, ".png"),
       sd.gamma.wtl.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)

#### Persistence ----
##### Cohesion of tree cover ----
# Generate new data
sd.coh <- seq(min(easo_df_thres_unsc$cohesion, na.rm = TRUE), max(easo_df_thres_unsc$cohesion, na.rm = TRUE), length.out = 50)
sd.coh_mean <- mean(easo_df_thres_unsc$cohesion, na.rm = TRUE); sd.coh_sd <- sd(easo_df_thres_unsc$cohesion, na.rm = TRUE)
sd.coh_s <- (sd.coh - sd.coh_mean) / sd.coh_sd

# Separate out focal mcmc chains
H_int = combine.mcmc(my_mod, vars = "H_occ[1]")
H_coh = combine.mcmc(my_mod, vars = "H_occ[2]")

# Predict new data
phi_coh_pred <- matrix(NA, nrow = length(sd.coh_s), ncol = length(H_coh))
for (i in seq_along(sd.coh_s)) {
  phi_coh_pred[i, ] <- 1/(1 + exp(-1*(H_int + H_coh*sd.coh_s[i])))
}

newdat.sd.phi.coh <- data.frame(phi = apply(phi_coh_pred, 1, mean), 
                            psi1_lci_95 = apply(phi_coh_pred, 1, function(x) quantile(x, 0.025)),
                            psi1_uci_95 = apply(phi_coh_pred, 1, function(x) quantile(x, 0.975)),
                            psi1_lci_80 = apply(phi_coh_pred, 1, function(x) quantile(x, 0.1)),
                            psi1_uci_80 = apply(phi_coh_pred, 1, function(x) quantile(x, 0.9)),
                            sd.coh = sd.coh)

# Plot predictions
(sd.phi.coh.plot <- ggplot(newdat.sd.phi.coh, aes(x = sd.coh, y = phi)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Cohesion of tree cover (%)", 
         y = bquote(Probability ~of ~persistence ~(phi))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_phi_coh_", date, ".png"),
       sd.phi.coh.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)


##### Average winter minimum temperature ----
# Generate new data
sd.wt <- seq(min(easo_df_thres_unsc$avg.winter.tmin, na.rm = TRUE), max(easo_df_thres_unsc$avg.winter.tmin, na.rm = TRUE), length.out = 50)
sd.wt_mean <- mean(easo_df_thres_unsc$avg.winter.tmin, na.rm = TRUE); sd.wt_sd <- sd(easo_df_thres_unsc$avg.winter.tmin, na.rm = TRUE)
sd.wt_s <- (sd.wt - sd.wt_mean) / sd.wt_sd

# Separate out focal mcmc chains
H_int = combine.mcmc(my_mod, vars = "H_occ[1]")
H_wt = combine.mcmc(my_mod, vars = "H_occ[4]")

# Predict new data
phi_wt_pred <- matrix(NA, nrow = length(sd.wt_s), ncol = length(H_wt))
for (i in seq_along(sd.wt_s)) {
  phi_wt_pred[i, ] <- 1/(1 + exp(-1*(H_int + H_wt*sd.wt_s[i])))
}

newdat.sd.wt <- data.frame(phi = apply(phi_wt_pred, 1, mean), 
                           psi1_lci_95 = apply(phi_wt_pred, 1, function(x) quantile(x, 0.025)),
                           psi1_uci_95 = apply(phi_wt_pred, 1, function(x) quantile(x, 0.975)),
                           psi1_lci_80 = apply(phi_wt_pred, 1, function(x) quantile(x, 0.1)),
                           psi1_uci_80 = apply(phi_wt_pred, 1, function(x) quantile(x, 0.9)),
                           sd.wt = sd.wt)

# Plot predictions
(sd.phi.wt.plot <- ggplot(newdat.sd.wt, aes(x = sd.wt, y = phi)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Average winter minimum temperature (°C)", 
         y = bquote(Probability ~of ~persistence ~(phi))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_phi_wt_", date, ".png"),
       sd.phi.wt.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)

##### Average winter minimum temperature - lag 1 year ----
# Generate new data
sd.wtl <- seq(min(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE), max(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE), length.out = 50)
sd.wtl_mean <- mean(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE); sd.wtl_sd <- sd(easo_df_thres_unsc$lag.avg.winter.tmin, na.rm = TRUE)
sd.wtl_s <- (sd.wtl - sd.wtl_mean) / sd.wtl_sd

# Separate out focal mcmc chains
H_int = combine.mcmc(my_mod, vars = "H_occ[1]")
H_wtl = combine.mcmc(my_mod, vars = "H_occ[5]")

# Predict new data
phi_wtl_pred <- matrix(NA, nrow = length(sd.wtl_s), ncol = length(H_wtl))
for (i in seq_along(sd.wtl_s)) {
  phi_wtl_pred[i, ] <- 1/(1 + exp(-1*(H_int + H_wtl*sd.wtl_s[i])))
}

newdat.sd.wtl <- data.frame(phi = apply(phi_wtl_pred, 1, mean), 
                           psi1_lci_95 = apply(phi_wtl_pred, 1, function(x) quantile(x, 0.025)),
                           psi1_uci_95 = apply(phi_wtl_pred, 1, function(x) quantile(x, 0.975)),
                           psi1_lci_80 = apply(phi_wtl_pred, 1, function(x) quantile(x, 0.1)),
                           psi1_uci_80 = apply(phi_wtl_pred, 1, function(x) quantile(x, 0.9)),
                           sd.wtl = sd.wtl)

# Plot predictions
(sd.phi.wtl.plot <- ggplot(newdat.sd.wtl, aes(x = sd.wtl, y = phi)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Average winter minimum temperature (°C) - lag 1 year", 
         y = bquote(Probability ~of ~persistence ~(phi))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_phi_wtl_", date, ".png"),
       sd.phi.wtl.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)


##### Average breeding season temperature ----
# Generate new data
sd.bt <- seq(min(easo_df_thres_unsc$avg.breeding.tmin, na.rm = TRUE), max(easo_df_thres_unsc$avg.breeding.tmin, na.rm = TRUE), length.out = 50)
sd.bt_mean <- mean(easo_df_thres_unsc$avg.breeding.tmin, na.rm = TRUE); sd.bt_sd <- sd(easo_df_thres_unsc$avg.breeding.tmin, na.rm = TRUE)
sd.bt_s <- (sd.bt - sd.bt_mean) / sd.bt_sd

# Separate out focal mcmc chains
H_int = combine.mcmc(my_mod, vars = "H_occ[1]")
H_bt = combine.mcmc(my_mod, vars = "H_occ[6]")

# Predict new data
phi_bt_pred <- matrix(NA, nrow = length(sd.bt_s), ncol = length(H_bt))
for (i in seq_along(sd.bt_s)) {
  phi_bt_pred[i, ] <- 1/(1 + exp(-1*(H_int + H_bt*sd.bt_s[i])))
}

newdat.sd.bt <- data.frame(phi = apply(phi_bt_pred, 1, mean), 
                            psi1_lci_95 = apply(phi_bt_pred, 1, function(x) quantile(x, 0.025)),
                            psi1_uci_95 = apply(phi_bt_pred, 1, function(x) quantile(x, 0.975)),
                            psi1_lci_80 = apply(phi_bt_pred, 1, function(x) quantile(x, 0.1)),
                            psi1_uci_80 = apply(phi_bt_pred, 1, function(x) quantile(x, 0.9)),
                            sd.bt = sd.bt)

# Plot predictions
(sd.phi.bt.plot <- ggplot(newdat.sd.bt, aes(x = sd.bt, y = phi)) + 
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = psi1_lci_80, ymax = psi1_uci_80), alpha = 0.6) + 
    geom_ribbon(aes(ymin = psi1_lci_95, ymax = psi1_uci_95), alpha = 0.2) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
    labs(x = "Average breeding season minimum temperature (°C)", 
         y = bquote(Probability ~of ~persistence ~(phi))) +
    theme_classic() + 
    theme(axis.text = element_text(size = 15,colour = "black"), 
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18)))
ggsave(paste0(root, model_structure,"/", subf, "/plots/easo_", subf_ch, "_pred_phi_bt_", date, ".png"),
       sd.phi.bt.plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600)


