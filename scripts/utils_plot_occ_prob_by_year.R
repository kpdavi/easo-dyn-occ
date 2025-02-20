# Purpose: Function for creating ggplots of posterior means, 80%, and 95% credible intervals of model-estimated annual probabilities of occupancy and detecion processes averaged across sites and/or surveys from the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

plot.prob.annual.func <- function(df, param.foc, title.c, label.y.ax, year.range){
  
  # Filter to focal occupancy process / subprocess
  process_df <- df |>
    filter(param == param.foc)
  
  # Create line and ribbon plot
  process_plot <- ggplot(process_df, aes(year, mean)) +
    geom_line(linewidth = 1) + 
    geom_ribbon(aes(ymin = lci.95, ymax = uci.95), alpha = 0.2) +
    geom_ribbon(aes(ymin = lci.80, ymax = uci.80), alpha = 0.6) +
    labs(title = title.c,
         x = "Year",
         y = label.y.ax) +
    scale_x_continuous(breaks = year.range) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 22))
  
  return(process_plot)
  
}

# end script