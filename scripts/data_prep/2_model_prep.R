# Purpose: Prepare survey and covariate datasets to fit the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
pacman::p_load(
  here,
  tidyverse
)


# Import data ----
surveys_cov_det <- read.csv(here("data", "processed", "surveys_covs_det.csv"), header = TRUE)
cov_tree <- read.csv(here("data", "processed", "covs_tree.csv"), header = TRUE)
cov_climate <- read.csv(here("data", "processed", "covs_climate.csv"), header = TRUE)
cov_lidar <- read.csv(here("data", "processed", "covs_lidar.csv"), header = TRUE)

# Focal covariates ----
cols_id <- c("Transect", "year", "site_id", "year_id") # non-covariate columns

# Covariate names
cov_names_tree <- setdiff(names(cov_tree), cols_id)
cov_names_lidar <- setdiff(names(cov_lidar), cols_id)
cov_names_climate <- setdiff(names(cov_climate), cols_id)
(cov_names <- data.frame(names_cov = c("intercept", cov_names_tree, cov_names_lidar, cov_names_climate)))

cov_names$names_real <- c("Intercept",
                          "Proportion of tree cover",
                          "Cohesion of tree cover (%)",
                          "Standard deviation in vegetation height (m)",
                          "Breeding season minimum temperature (\u00B0C)",
                          "Breeding season minimum temperature lag 1 year (\u00B0C)",
                          "Winter minimum temperature (\u00B0C)",
                          "Winter minimum temperature lag 1 year (\u00B0C)",
                          "Winter precipitation (mm)"
                          )

## Identify focal covariates by process ----
### Specify covariates
cov_psi1_mod <- c("intercept", "tree.prop", "sd.vegheight.m.2013")
cov_gamma_mod <- c("intercept", "tree.coh", "sum.winter.prcp", "avg.winter.tmin", "lag.avg.winter.tmin", "avg.breeding.tmin")
cov_phi_mod <- c("intercept", "tree.coh", "sum.winter.prcp", "avg.winter.tmin", "lag.avg.winter.tmin", "avg.breeding.tmin", "lag.avg.breeding.tmin")

### Create and format covariate name information (for filtering and plotting later)
#### Initial occupancy
cov_psi1_foc <- cov_names |>
  filter(names_cov %in% cov_psi1_mod) |>
  mutate(names_cov = factor(names_cov, levels = cov_psi1_mod)) %>%
  arrange(names_cov)
cov_psi1_param_name <- paste0("P_occ[", seq(1:nrow(cov_psi1_foc)), "]")
cov_psi1_coeffs <- cbind(cov_psi1_foc,
                         data.frame(param = cov_psi1_param_name,
                                    proc = "Initial occupancy"))

#### Colonization
cov_gamma_foc <- cov_names |>
  filter(names_cov %in% cov_gamma_mod) |>
  mutate(names_cov = factor(names_cov, levels = cov_gamma_mod)) %>%
  arrange(names_cov)
cov_gamma_param_name <- paste0("G_occ[", seq(1:nrow(cov_gamma_foc)), "]")
cov_gamma_coeffs <- cbind(cov_gamma_foc,
                         data.frame(param = cov_gamma_param_name,
                                    proc = "Colonization"))

#### Persistence
cov_phi_foc <- cov_names |>
  filter(names_cov %in% cov_phi_mod) |>
  mutate(names_cov = factor(names_cov, levels = cov_phi_mod)) %>%
  arrange(names_cov)
cov_phi_param_name <- paste0("H_occ[", seq(1:nrow(cov_phi_foc)), "]")
cov_phi_coeffs <- cbind(cov_phi_foc,
                          data.frame(param = cov_phi_param_name,
                                     proc = "Persistence"))

## Coefficients ----
coeff_names <- rbind(cov_psi1_coeffs, cov_gamma_coeffs, cov_phi_coeffs)


# Prepare data for model fitting -----
## Scale covariates ----
### Detection ----
cov_names_det_sc <- c("moon.phase", "date.ord")
cov_names_det_nonsc <- setdiff(names(surveys_cov_det), cov_names_det_sc)

surveys_cov_det_sc <- data.frame(surveys_cov_det[cov_names_det_nonsc], lapply(surveys_cov_det[cov_names_det_sc], scale))

### Occupancy processes ----
#### Tree cover ----
cov_names_tree_sc <- c("tree.prop", "tree.coh")
cov_names_tree_nonsc <- setdiff(names(cov_tree), cov_names_tree_sc)
cov_tree_sc <- data.frame(cov_tree[cov_names_tree_nonsc], lapply(cov_tree[cov_names_tree_sc], scale))

#### Vegetation height ----
cov_lidar_sc <- data.frame(cov_lidar[c("Transect", "year")], lapply(cov_lidar["sd.vegheight.m.2013"], scale))

#### Climate ----
cov_names_climate_sc <- c("avg.breeding.tmin", "lag.avg.breeding.tmin", "avg.winter.tmin", "lag.avg.winter.tmin", "sum.winter.prcp")
cov_names_climate_nonsc <- setdiff(names(cov_climate), cov_names_climate_sc)

cov_climate_sc <- data.frame(cov_climate[cov_names_climate_nonsc], lapply(cov_climate[cov_names_climate_sc], scale))


## Format data for the dynamic occupancy model ----
### Detection/visits -----
easo_obs <- surveys_cov_det_sc |>
  arrange(Transect, year) |>
  mutate(intercept = 1,
         site_id = as.numeric(factor(Transect)),
         year_id = as.numeric(factor(year))) |>
  dplyr::select(Transect, site_id, year, year_id, presence, intercept, moon.phase, date.ord)
cov_det_num <- as.numeric(ncol(dplyr::select(easo_obs, -c("Transect", "site_id", "year", "year_id", "presence"))))
X_det <- as.matrix(dplyr::select(easo_obs, -c("Transect", "site_id", "year", "year_id", "presence")))

easo_surv_sites <- easo_obs |>
  arrange(site_id, year_id) |>
  group_by(site_id, year_id) |>
  summarise(site_yr_pres = max(presence))
easo_surv_sites$site_surv_id <- seq(1, nrow(easo_surv_sites), 1)

### Initial occupancy ----
tree_prop_2013 <- cov_tree_sc |>
  filter(year == 2013) |>
  dplyr::select(-tree.coh)

cov_psi1 <- plyr::join_all(list(tree_prop_2013, cov_lidar_sc), by=c("Transect", "year"), type = 'left')
cov_psi1 <- cov_psi1 |>
  arrange(Transect, year) |>
  mutate(intercept = 1) |>
  dplyr::select(Transect, year, intercept, tree.prop, sd.vegheight.m.2013)

cov_psi1_num <- as.numeric(ncol(dplyr::select(cov_psi1, all_of(cov_psi1_mod))))
cov_psi1_names <- names(dplyr::select(cov_psi1, all_of(cov_psi1_mod)))
cov_psi1_list <- vector("list", cov_psi1_num)

position_lidar <- as.numeric(which(cov_psi1_names == "sd.vegheight.m.2013"))

for(i in 1:cov_psi1_num) {
  # Select focal covariate and create site-year matrix
  cov_psi1_foc <- cov_psi1 |>
    arrange(Transect, year) |>
    pivot_wider(id_cols = Transect, names_from = year, values_from = cov_psi1_names[i]) |>
    dplyr::select(-Transect)

  # Assign df to list and name it as the focal covariate
  cov_psi1_list[[i]] <- cov_psi1_foc
  names(cov_psi1_list)[i] <- cov_psi1_names[i]
}

# Create array
X_psi1 = array(unlist(cov_psi1_list), dim = c(as.numeric(nrow(cov_psi1_foc)), as.numeric(ncol(cov_psi1_foc)), cov_psi1_num))


### Colonization ----
cov_gamma <- cov_tree_sc |>
  arrange(Transect, year) |>
  left_join(cov_climate_sc, by = "year") |>
  mutate(intercept = 1) |>
  dplyr::select(Transect, year, all_of(cov_gamma_mod))
cov_gamma_num <- as.numeric(ncol(dplyr::select(cov_gamma, -c("Transect", "year"))))
cov_gamma_names <- names(dplyr::select(cov_gamma, -c("Transect", "year")))

cov_gamma_list <- vector("list", cov_gamma_num)
for(i in 1:cov_gamma_num) {
  # Select focal covariate and create site-year matrix
  cov_gamma_foc <- cov_gamma |>
    arrange(Transect, year) |>
    pivot_wider(id_cols = Transect, names_from = year, values_from = cov_gamma_names[i]) |>
    dplyr::select(-Transect)

  # Assign df to list and name it as the focal covariate
  cov_gamma_list[[i]] <- cov_gamma_foc
  names(cov_gamma_list)[i] <- cov_gamma_names[i]
}

# Create array
X_gam <- array(unlist(cov_gamma_list), dim = c(as.numeric(nrow(cov_gamma_foc)), as.numeric(ncol(cov_gamma_foc)), cov_gamma_num))


### Persistence ----
cov_phi <- cov_tree_sc |>
  arrange(Transect, year) |>
  left_join(cov_climate_sc, by = "year") |>
  mutate(intercept = 1) |>
  dplyr::select(Transect, year, all_of(cov_phi_mod))
cov_phi_num <- as.numeric(ncol(dplyr::select(cov_phi, -c("Transect", "year"))))
cov_phi_names <- names(dplyr::select(cov_phi, -c("Transect", "year")))

cov_phi_list <- vector("list", cov_phi_num)
for(i in 1:cov_phi_num) {
  # Select focal covariate and create site-year matrix
  cov_phi_foc = cov_phi |>
    arrange(Transect, year) |>
    pivot_wider(id_cols = Transect, names_from = year, values_from = cov_phi_names[i]) |>
    dplyr::select(-Transect)

  # Assign df to list and name it as the focal covariate
  cov_phi_list[[i]] <- cov_phi_foc
  names(cov_phi_list)[i] <- cov_phi_names[i]
}

# Create array
X_phi <- array(unlist(cov_phi_list), dim = c(as.numeric(nrow(cov_phi_foc)), as.numeric(ncol(cov_phi_foc)), cov_phi_num))


# Combine objects for model fitting and save ----
save(easo_obs, cov_det_num, cov_psi1_num, cov_gamma_num, cov_phi_num, X_det, X_psi1, X_phi, X_gam, position_lidar, coeff_names, file = here("data", "processed", "model_data.RData"))

# end script
