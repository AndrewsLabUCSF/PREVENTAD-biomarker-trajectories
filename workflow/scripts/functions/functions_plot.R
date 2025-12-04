
# DIAGNOSTIC PLOTS --------------------------------------------------------

create_diagnostic_plots <- function(model_obj, title, outcome_name) {
  
  # Get residuals
  resids <- residuals(model_obj)
  fitted <- fitted(model_obj)
  
  # 1. Histogram of residuals
  p1 <- ggplot(data.frame(residuals = resids), aes(x = residuals)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Histogram of Residuals"),
         subtitle = title,
         x = "Residuals", y = "Count") +
    theme_minimal()
  
  # 2. Q-Q plot
  p2 <- ggplot(data.frame(residuals = resids), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red", linewidth = 1) +
    labs(title = "Q-Q Plot",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # 3. Residuals vs Fitted
  p3 <- ggplot(data.frame(fitted = fitted, residuals = resids), 
               aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(se = FALSE, color = "blue") +
    labs(title = "Residuals vs Fitted",
         x = "Fitted Values", y = "Residuals") +
    theme_minimal()
  
  # Combine all 4 plots
  combined <- (p1 | p2 | p3) +
    plot_annotation(
      title = paste(outcome_name, "-", title),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  return(combined)
}

