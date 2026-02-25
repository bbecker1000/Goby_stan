# --- 1. Get the Mean of Posterior Draws for Each Parameter ---
mean_params <- colMeans(relevant_params)
mean_params_list <- as.list(mean_params)

# --- 2. Predict Density for Each Raw Data Point (CORRECTED) ---
predicted_densities <- numeric(nrow(df.data))

for (i in 1:nrow(df.data)) {
  current_row_data <- as.list(df.data[i, ])
  
  # This returns the linear predictor for log(COUNT)
  mu_goby_pred_log_count <- calculate_mu_goby_deterministic(
    draw_row = mean_params_list,
    input_data = current_row_data,
    initial_intervention_var = NULL,
    initial_intervention_value = NULL
  )
  
  # --- THE FIX ---
  # Convert log(count) to log(density) by subtracting the log(Area) offset
  mu_goby_pred_log_density <- mu_goby_pred_log_count - current_row_data$Area
  
  # Exponentiate the log(density) to get the predicted density
  predicted_densities[i] <- exp(mu_goby_pred_log_density)
}

df.data$predicted_density <- predicted_densities

# --- 3. Create the Corrected Diagnostic Plot ---
diagnostic_plot <- ggplot(df.data, aes(x = Goby / exp(Area), y = predicted_density)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "CORRECTED Posterior Predictive Check: Observed vs. Predicted Density",
    subtitle = "Points should now cluster around the red 1:1 line.",
    x = "Observed Goby Density (Goby / Area)",
    y = "Predicted Goby Density (from model)"
  ) +
  theme_minimal(base_size = 14) +
  coord_equal(xlim = range(c(df.data$Goby / exp(df.data$Area), df.data$predicted_density), na.rm=TRUE),
              ylim = range(c(df.data$Goby / exp(df.data$Area), df.data$predicted_density), na.rm=TRUE))

print(diagnostic_plot)
