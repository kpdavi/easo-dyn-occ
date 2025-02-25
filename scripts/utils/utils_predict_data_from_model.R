# Purpose: Function to create predictions from the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# This function is called within the function for plotting predictions 
# Author: Kristin P. Davis

predict.cov.data = function(mod_output, cov_df, cov_col, proc_char, cov_position, cov_name_real, proc_name_real, pred_df_name, pred_df_name_sum) {
  
  # Generate and scale a vector of covariate values
  cov_vec <- seq(min(cov_df[[cov_col]], na.rm = TRUE), max(cov_df[[cov_col]], na.rm = TRUE), length.out = 50)
  cov_mean <- mean(cov_df[[cov_col]], na.rm = TRUE) 
  cov_sd <- sd(cov_df[[cov_col]], na.rm = TRUE)
  cov_sc <- (cov_vec - cov_mean) / cov_sd

  # Separate out focal mcmc chains
  proc_int <- combine.mcmc(mod_output, vars = paste0(proc_char, "[1]"))
  proc_cov <- combine.mcmc(mod_output, vars = paste0(proc_char, "[", as.character(cov_position), "]"))

  # Predict new data
  pred_proc_cov <- matrix(NA, nrow = length(cov_sc), ncol = length(proc_cov))
  for (i in seq_along(cov_sc)) {
    pred_proc_cov[i, ] <- 1/(1 + exp(-1*(proc_int + proc_cov * cov_sc[i])))
  }

  # Calculate posterior mean and credible intervals
  pred_proc_cov_sum <- data.frame(
    y = apply(pred_proc_cov, 1, mean), 
    lci_95 = apply(pred_proc_cov, 1, function(x) quantile(x, 0.025)),
    uci_95 = apply(pred_proc_cov, 1, function(x) quantile(x, 0.975)),
    lci_80 = apply(pred_proc_cov, 1, function(x) quantile(x, 0.1)),
    uci_80 = apply(pred_proc_cov, 1, function(x) quantile(x, 0.9)),
    cov.val = cov_vec,
    cov = cov_name_real,
    proc = proc_name_real)
  
  # Assign the new data to a dynamic name in the global environment
  assign(pred_df_name, pred_proc_cov, envir = .GlobalEnv)
  assign(pred_df_name_sum, pred_proc_cov_sum, envir = .GlobalEnv)
  
  }

# end script