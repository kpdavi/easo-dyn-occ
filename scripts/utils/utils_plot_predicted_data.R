# Purpose: Function to plot predictions from the dynamic occupancy model for Eastern Screech-Owl in Fort Collins, CO
# Author: Kristin P. Davis

plot.cov.preds <- function(pred_cov_df, label_x, process_text, process_greek_symbol, add_subscript = FALSE) {
  
  plot_pred_cov <- ggplot(pred_cov_df, aes(x = cov.val, y = y)) + 
     geom_line(linewidth = 1) + 
     geom_ribbon(aes(ymin = lci_80, ymax = uci_80), alpha = 0.6) + 
     geom_ribbon(aes(ymin = lci_95, ymax = uci_95), alpha = 0.2) + 
     scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .25)) +
     labs(
       x = label_x, 
       y = if (add_subscript) {
         bquote(.(process_text) ~ "(" * .(as.name(process_greek_symbol))[1] * ")")  # Apply subscript if TRUE
       } else {
         bquote(.(process_text) ~ "(" * .(as.name(process_greek_symbol)) * ")")  # No subscript
       }
     ) +
     theme_classic() + 
     theme(axis.text = element_text(size = 15,colour = "black"), 
           axis.title = element_text(size = 18),
           plot.title = element_text(size = 18))

  return(plot_pred_cov)
  
  }

# end script