# Purpose: Create figure 2 from the dynamic occupancy model analysis for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

# Required packages ----
# install.packages("pacman")
pacman::p_load(
  tidyverse,
  runjags
)

# Import data and model ouput ----
load("data/processed/model_data.RData")
mod_out <- readRDS("output/mod_out_coeffs.rds")


# Calculate posterior means and credible intervals ----
## Extract posterior distributions for coefficient parameters
mod_param_ests = dplyr::select(as.data.frame(as.matrix(mod_out$mcmc)), matches("^G_|^P_|^H_"))

## Calculate means and credible intervals (80% and 95% credible intervals)
mod_mean <- apply(mod_param_ests, 2, mean)
mod_lci_80 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.1)))
mod_uci_80 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.9)))
mod_lci_95 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.025)))
mod_uci_95 <- apply(mod_param_ests, 2, function(x) quantile(x, probs = c(0.975)))

## Use replacement vector to change column names to covariate names
coeff_replacement_vector <- coeff_names |>
  mutate(name_long = paste0(proc, "-", names_real)) |>
  select(param, name_long) |>
  deframe()

mod_param_ests_renamed <- mod_param_ests |>
  rename_with(~ coeff_replacement_vector)


# Format posterior information for plotting ----
mod_cat_dat_df <- data.frame(coeff.names = colnames(mod_param_ests_renamed),
                          mn = mod_mean, 
                          lci.80 = mod_lci_80,
                          uci.80 = mod_uci_80,
                          lci.95 = mod_lci_95,
                          uci.95 = mod_uci_95)
mod_cat_dat <- mod_cat_dat_df |>
  rownames_to_column() |>
  rename(model_param = rowname) |>
  separate(coeff.names, into = c("proc_name", "coeff"), sep = "-") |>
  mutate(proc_name = fct_relevel(proc_name, c("Persistence", "Colonization", "Initial occupancy"))) |>
  group_by(proc_name) |>
  mutate(coeff = factor(coeff, levels = c("Intercept", setdiff(unique(coeff), "Intercept"))),
         coeff = fct_rev(coeff))


# Create caterpillar plot ----
## Facet labels
cat_plot_labels <- as.vector(c("a) Persistence", "b) Colonization", "c) Initial occupancy"))
names(cat_plot_labels) <- c("Persistence", "Colonization", "Initial occupancy")

## Plot
(plot_cat_coeff <- ggplot(mod_cat_dat, aes(x = mn, y = coeff)) + 
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

## Save plot
ggsave(filename = "output/plots/fig2_coeff_param_ests.png",
       plot = plot_cat_coeff,
       width = 23, 
       height = 6.5, 
       units = "in",
       dpi = 600)

# end script