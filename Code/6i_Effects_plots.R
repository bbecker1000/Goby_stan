#effects plots

df.data <- dat

#hist(df.data$Goby/exp(df.data$Area))

# Add this line before your ggplot code
df.data$Goby_Density <- df.data$Goby / exp(df.data$Area)


## and now the effects plots on goby density: 

##updated effects--------

# --- 0. Setup ---
library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(forcats)
library(stringr)
library(cowplot)

# --- 1. Extract Posterior Draws ---
# This assumes your 'mod.SEM.real' and 'df.data' are loaded correctly
draws <- mod.SEM.real$draws(format = "df")
relevant_params <- draws

# --- 2. Deterministic Linear Predictor Function (CORRECTED) ---
calculate_mu_goby_deterministic <- function(draw_row, input_data, initial_intervention_var = NULL, initial_intervention_value = NULL) {
  
  # Initialize current values
  current_values <- as.list(input_data)
  
  # Helper to update squares
  update_squares <- function(vals) {
    if(!is.null(vals$BreachDays)) vals$BreachDays_2 <- vals$BreachDays^2
    if(!is.null(vals$Temp))       vals$Temp_2       <- vals$Temp^2
    if(!is.null(vals$SAV))        vals$SAV_2        <- vals$SAV^2
    if(!is.null(vals$Year))       vals$Year_2       <- vals$Year^2
    return(vals)
  }
  
  # 1. APPLY INITIAL INTERVENTION
  if (!is.null(initial_intervention_var)) {
    current_values[[initial_intervention_var]] <- initial_intervention_value
    # IMMEDIATE FIX: Update the square of the intervened variable immediately
    current_values <- update_squares(current_values)
  }
  
  get_val <- function(var_name) {
    return(current_values[[var_name]])
  }
  
  is_intervened <- function(var) !is.null(initial_intervention_var) && initial_intervention_var == var
  
  # 2. PROPAGATE THROUGH THE DAG
  # IMPORTANT: Every time we calculate a new downstream variable, we MUST update its square
  
  # --- BreachDays ---
  if (!is_intervened("BreachDays")) {
    Breach_nu <- draw_row[["a_BreachDays"]] + 
      draw_row[["beta_Rain"]] * get_val("Rain") + 
      draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$BreachDays <- Breach_nu
    current_values$BreachDays_2 <- Breach_nu^2 # <--- FIX ADDED
  }
  
  # --- Temp (Depends on BreachDays) ---
  if (!is_intervened("Temp")) {
    Temp_nu <- draw_row[["a_Temp"]] + 
      draw_row[["beta_BreachDays_vec[2]"]] * get_val("BreachDays") + 
      draw_row[["beta_Wind_vec[2]"]] * get_val("Wind") + 
      draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$Temp <- Temp_nu
    current_values$Temp_2 <- Temp_nu^2 # <--- FIX ADDED
  }
  
  # --- DO (Depends on Temp) ---
  if (!is_intervened("DO")) {
    DO_nu <- draw_row[["a_DO"]] + 
      draw_row[["beta_Temp_vec[2]"]] * get_val("Temp") + 
      draw_row[["beta_Wind_vec[1]"]] * get_val("Wind") + 
      draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$DO <- DO_nu
    # DO has no square term in your model, so no update needed
  }
  
  # --- SAV (Depends on DO, Temp) ---
  if (!is_intervened("SAV")) {
    SAV_nu <- draw_row[["a_SAV"]] + 
      draw_row[["beta_DO_vec[4]"]] * get_val("DO") + 
      draw_row[["beta_Temp_vec[3]"]] * get_val("Temp") + 
      draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$SAV <- SAV_nu
    current_values$SAV_2 <- SAV_nu^2 # <--- FIX ADDED
  }
  
  # --- SC & SB (Binary) ---
  # Note: For strict Bayesian accuracy, one should integrate over the probability (0/1). 
  # However, using the probability (plogis) as an expectation is a common approximation for effects plots.
  if (!is_intervened("SC_count")) {
    SC_mu_logit <- draw_row[["a_SC"]] + 
      draw_row[["beta_Substrate_vec[2]"]] * get_val("Substrate") + 
      draw_row[["beta_DO_vec[3]"]] * get_val("DO") + 
      draw_row[["beta_SAV_vec[3]"]] * get_val("SAV") + 
      draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$SC_count <- plogis(SC_mu_logit)
  }
  
  if (!is_intervened("SB_count")) {
    SB_mu_logit <- draw_row[["a_SB"]] + 
      draw_row[["beta_DO_vec[2]"]] * get_val("DO") + 
      draw_row[["beta_SAV_vec[2]"]] * get_val("SAV") + 
      draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
    current_values$SB_count <- plogis(SB_mu_logit)
  }
  
  # 3. CALCULATE FINAL GOBY EXPECTATION
  # Now get_val("Temp_2") will return the CORRECT squared value calculated above
  mu_goby_val <- draw_row[[paste0("a_Goby[", current_values$Zone, "]")]] +
    draw_row[["beta_SAV_vec[1]"]] * get_val("SAV") +
    draw_row[["beta_DO_vec[1]"]] * get_val("DO") +
    draw_row[["beta_BreachDays_vec[1]"]] * get_val("BreachDays") +
    draw_row[["beta_Temp_vec[1]"]] * get_val("Temp") +
    draw_row[["beta_SC_count"]] * get_val("SC_count") +
    draw_row[["beta_SB_count"]] * get_val("SB_count") +
    draw_row[["beta_Substrate_vec[1]"]] * get_val("Substrate") +
    draw_row[["beta_Micro"]] * get_val("Micro") +
    draw_row[["beta_Year"]] * get_val("Year") +
    draw_row[["beta_Year_2"]] * get_val("Year_2") +       # This was previously using the static median
    draw_row[["beta_Temp_2"]] * get_val("Temp_2") +       # This was previously using the static median
    draw_row[["beta_BreachDays_2"]] * get_val("BreachDays_2") + # This was previously using the static median
    draw_row[["beta_Goby_lag"]] * get_val("Goby_lag") +
    draw_row[["beta_SAV_2"]] * get_val("SAV_2") +         # Added SAV_2
    current_values$Area 
  
  return(mu_goby_val)
}



# # --- 2. Deterministic Linear Predictor Function (Unchanged) ---
# calculate_mu_goby_deterministic <- function(draw_row, input_data, initial_intervention_var = NULL, initial_intervention_value = NULL) {
#   # This function correctly calculates log(COUNT) and remains unchanged.
#   current_values <- as.list(input_data)
#   if (!is.null(initial_intervention_var)) {
#     current_values[[initial_intervention_var]] <- initial_intervention_value
#   }
#   
#   get_val <- function(var_name) {
#     return(current_values[[var_name]])
#   }
#   
#   is_intervened <- function(var) !is.null(initial_intervention_var) && initial_intervention_var == var
#   
#   if (!is_intervened("BreachDays")) {
#     Breach_nu <- draw_row[["a_BreachDays"]] + draw_row[["beta_Rain"]] * get_val("Rain") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$BreachDays <- Breach_nu
#   }
#   if (!is_intervened("Temp")) {
#     Temp_nu <- draw_row[["a_Temp"]] + draw_row[["beta_BreachDays_vec[2]"]] * get_val("BreachDays") + draw_row[["beta_Wind_vec[2]"]] * get_val("Wind") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$Temp <- Temp_nu
#   }
#   if (!is_intervened("DO")) {
#     DO_nu <- draw_row[["a_DO"]] + draw_row[["beta_Temp_vec[2]"]] * get_val("Temp") + draw_row[["beta_Wind_vec[1]"]] * get_val("Wind") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$DO <- DO_nu
#   }
#   if (!is_intervened("SAV")) {
#     SAV_nu <- draw_row[["a_SAV"]] + draw_row[["beta_DO_vec[4]"]] * get_val("DO") + draw_row[["beta_Temp_vec[3]"]] * get_val("Temp") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SAV <- SAV_nu
#   }
#   if (!is_intervened("SC_count")) {
#     SC_mu_logit <- draw_row[["a_SC"]] + draw_row[["beta_Substrate_vec[2]"]] * get_val("Substrate") + draw_row[["beta_DO_vec[3]"]] * get_val("DO") + draw_row[["beta_SAV_vec[3]"]] * get_val("SAV") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SC_count <- plogis(SC_mu_logit)
#   }
#   if (!is_intervened("SB_count")) {
#     SB_mu_logit <- draw_row[["a_SB"]] + draw_row[["beta_DO_vec[2]"]] * get_val("DO") + draw_row[["beta_SAV_vec[2]"]] * get_val("SAV") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SB_count <- plogis(SB_mu_logit)
#   }
#   
#   mu_goby_val <- draw_row[[paste0("a_Goby[", current_values$Zone, "]")]] +
#     draw_row[["beta_SAV_vec[1]"]] * get_val("SAV") +
#     draw_row[["beta_DO_vec[1]"]] * get_val("DO") +
#     draw_row[["beta_BreachDays_vec[1]"]] * get_val("BreachDays") +
#     draw_row[["beta_Temp_vec[1]"]] * get_val("Temp") +
#     draw_row[["beta_SC_count"]] * get_val("SC_count") +
#     draw_row[["beta_SB_count"]] * get_val("SB_count") +
#     draw_row[["beta_Substrate_vec[1]"]] * get_val("Substrate") +
#     draw_row[["beta_Micro"]] * get_val("Micro") +
#     draw_row[["beta_Year"]] * get_val("Year") +
#     draw_row[["beta_Year_2"]] * get_val("Year_2") +
#     draw_row[["beta_Temp_2"]] * get_val("Temp_2") +
#     draw_row[["beta_BreachDays_2"]] * get_val("BreachDays_2") +
#     draw_row[["beta_Goby_lag"]] * get_val("Goby_lag") +
#     current_values$Area 
#   
#   return(mu_goby_val)
# }

# --- 3. Marginal Effects Calculation Function (CORRECTED) ---
calculate_marginal_effect_goby_by_zone <- function(
    vary_var,                   
    vary_vals,                  
    draws_df,                   
    data_df,                    
    target_zone,                
    n_draws = 40             
) {
  
  message(paste("Calculating marginal effect for:", vary_var, "in Zone", target_zone))
  
  baseline_data_row <- data_df %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~median(.x, na.rm = TRUE))) %>% 
    dplyr::mutate(Zone = as.integer(target_zone)) %>%
    as.data.frame()
  
  total_draws_available <- nrow(draws_df)
  if (n_draws > total_draws_available) {
    n_draws <- total_draws_available
  }
  sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
  selected_draws <- draws_df[sample_indices, ]
  
  results_df <- tibble(
    var_value = rep(vary_vals, each = n_draws),
    draw_idx = rep(1:n_draws, times = length(vary_vals)),
    mu_goby_pred_log_density = numeric(length(vary_vals) * n_draws) # Storing log(DENSITY) now
  )
  
  result_counter <- 1
  for (val in vary_vals) {
    for (d in 1:n_draws) {
      draw_row <- selected_draws[d, ]
      intervened_input_data <- baseline_data_row
      
      if (vary_var == "Year") {
        intervened_input_data[["Year"]] <- val
        if ("Year_2" %in% names(intervened_input_data)) intervened_input_data[["Year_2"]] <- val^2
      } else if (vary_var == "BreachDays") {
        intervened_input_data[["BreachDays"]] <- val
        if ("BreachDays_2" %in% names(intervened_input_data)) intervened_input_data[["BreachDays_2"]] <- val^2
      } 
      
      mu_goby_pred_log_count <- calculate_mu_goby_deterministic(
        draw_row = draw_row,
        input_data = intervened_input_data,
        initial_intervention_var = vary_var, 
        initial_intervention_value = val      
      )
      
      # --- THE FIX ---
      # Convert log(count) to log(density) by subtracting the log(Area) offset
      log_area_offset <- baseline_data_row$Area # This is held at median for predictions
      mu_goby_pred_log_density <- mu_goby_pred_log_count - log_area_offset
      
      results_df$mu_goby_pred_log_density[result_counter] <- mu_goby_pred_log_density
      result_counter <- result_counter + 1
    }
  }
  
  summary_df <- results_df %>%
    dplyr::group_by(var_value) %>%
    dplyr::summarise(
      mean_log_density = mean(mu_goby_pred_log_density, na.rm = TRUE),
      lower_log_density = quantile(mu_goby_pred_log_density, 0.055, na.rm = TRUE),
      upper_log_density = quantile(mu_goby_pred_log_density, 0.945, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      mean_goby_density = exp(mean_log_density), 
      lower_goby_density = exp(lower_log_density), 
      upper_goby_density = exp(upper_log_density),
      Zone = target_zone,
      predictor = vary_var
    )
  
  return(summary_df)
}

# --- 4. Define Variables and Their Ranges for Plots ---
sb_vals <- c(0, 1)
sc_vals <- c(0, 1)
breachdays_vals <- seq(min(df.data$BreachDays, na.rm = TRUE), max(df.data$BreachDays, na.rm = TRUE), length.out = 50)
rain_vals <- seq(min(df.data$Rain, na.rm = TRUE), max(df.data$Rain, na.rm = TRUE), length.out = 50)
do_vals <- seq(min(df.data$DO, na.rm = TRUE), max(df.data$DO, na.rm = TRUE), length.out = 50)
sav_vals <- seq(min(df.data$SAV, na.rm = TRUE), max(df.data$SAV, na.rm = TRUE), length.out = 50)
year_vals <- seq(min(df.data$Year, na.rm = TRUE), max(df.data$Year, na.rm = TRUE), length.out = 50)

# --- 5. Calculate Marginal Effects for Each Variable AND Each Zone ---
all_zones <- sort(unique(df.data$Zone))
all_marginal_effects <- list()

predictors_to_plot <- list(
  "SB_count" = sb_vals,
  "BreachDays" = breachdays_vals,
  "Rain" = rain_vals,
  "DO" = do_vals, 
  "SAV" = sav_vals, 
  "SC_count" = sc_vals,
  "Year" = year_vals
)

for (current_zone in all_zones) {
  for (predictor_name in names(predictors_to_plot)) {
    vals <- predictors_to_plot[[predictor_name]]
    
    me_result <- calculate_marginal_effect_goby_by_zone(
      vary_var = predictor_name,
      vary_vals = vals,
      draws_df = relevant_params,
      data_df = df.data,
      target_zone = current_zone
    )
    me_result$plot_id <- predictor_name
    all_marginal_effects[[length(all_marginal_effects) + 1]] <- me_result
  }
}

full_me_df <- do.call(rbind, all_marginal_effects)

# --- 6. Data Processing and Plot Creation ---
do_data_sb <- full_me_df %>% filter(plot_id == "DO") %>% mutate(plot_id = "DO_SB_Goby")
do_data_sc <- full_me_df %>% filter(plot_id == "DO") %>% mutate(plot_id = "DO_SC_Goby")
sav_data_sc <- full_me_df %>% filter(plot_id == "SAV") %>% mutate(plot_id = "SAV_SC_Goby")

full_me_df_processed <- full_me_df %>%
  filter(!plot_id %in% c("DO", "SAV")) %>% 
  mutate(
    plot_id = case_when(
      plot_id == "BreachDays" ~ "BreachDays_Goby",
      plot_id == "Rain" ~ "Rain_Goby",
      plot_id == "Year" ~ "Year2_Goby",
      TRUE ~ plot_id 
    )
  ) %>%
  bind_rows(do_data_sb, do_data_sc, sav_data_sc) 

zone_labels <- c("1" = "East", "2" = "Northwest", "3" = "West")

plot_titles <- c(
  "SB_count" = "SB \u2192 Goby",
  "BreachDays_Goby" = "BreachDays (Total Effect) \u2192 Goby",
  "Rain_Goby" = "Rain (Total Effect) \u2192 Goby",
  "DO_SB_Goby" = "DO \u2192 SB \u2192 Goby",
  "DO_SC_Goby" = "DO \u2192 SC \u2192 Goby",
  "SAV_SC_Goby" = "SAV \u2192 SC \u2192 Goby",
  "SC_count" = "SC \u2192 Goby",
  "Year2_Goby" = "Year & Year\u00B2 \u2192 Goby"
)

x_axis_labels <- c(
  "SB_count" = "SB (0 or 1)",
  "BreachDays" = "BreachDays",
  "Rain" = "Rain",
  "DO" = "Dissolved Oxygen (DO)",
  "SAV" = "Submerged Aquatic Vegetation (SAV)",
  "SC_count" = "SC (0 or 1)",
  "Year" = "Year"
)

raw_data_alpha <- 0.3
raw_data_color <- "grey30"
raw_data_size <- 0.8
jitter_width <- 0.2

make_marginal_plot <- function(plot_identifier, full_df, raw_df) {
  
  plot_data <- full_df %>% dplyr::filter(plot_id == plot_identifier)
  original_predictor_name <- as.character(plot_data$predictor[1])
  is_binary <- original_predictor_name %in% c("SB_count", "SC_count")
  
  p <- ggplot(plot_data, aes(x = var_value, y = mean_goby_density)) +
    geom_point(data = raw_df, aes_string(x = original_predictor_name, y = "Goby / exp(Area)"),
               alpha = raw_data_alpha, color = raw_data_color, size = raw_data_size,
               position = if (is_binary) position_jitter(width = jitter_width, height = 0) else "identity") +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_ribbon(aes(ymin = lower_goby_density, ymax = upper_goby_density), fill = "#0072B2", alpha = 0.2) +
    labs(
      title = plot_titles[plot_identifier],
      x = x_axis_labels[original_predictor_name],
      y = NULL # --- MODIFIED: Removed y-axis label from individual plots
    ) +
    theme_minimal(base_size = 16) +
    facet_wrap(~ Zone, 
               ncol = 3, labeller = labeller(Zone = zone_labels)) +
    coord_cartesian(ylim = c(0, 20)) # Sets the y-axis limits from 0 to 75
  
  if (is_binary) {
    p <- p + scale_x_continuous(breaks = c(0, 1))
  }
  
  return(p)
}

plot_sb_goby <- make_marginal_plot("SB_count", full_me_df_processed, df.data)
plot_breachdays_goby <- make_marginal_plot("BreachDays_Goby", full_me_df_processed, df.data)
plot_rain_goby <- make_marginal_plot("Rain_Goby", full_me_df_processed, df.data)
plot_do_sb_goby <- make_marginal_plot("DO_SB_Goby", full_me_df_processed, df.data)
plot_do_sc_goby <- make_marginal_plot("DO_SC_Goby", full_me_df_processed, df.data)
plot_sav_sc_goby <- make_marginal_plot("SAV_SC_Goby", full_me_df_processed, df.data)
plot_sc_goby <- make_marginal_plot("SC_count", full_me_df_processed, df.data)
plot_year2_goby <- make_marginal_plot("Year2_Goby", full_me_df_processed, df.data)

# --- 7. Arrange the 8 plots into a single figure and add common y-axis label ---

# First, create the grid of plots WITHOUT the common label

final_grid <- plot_grid(
  plot_rain_goby,
  plot_breachdays_goby,
  plot_sb_goby,
  plot_sc_goby,
  plot_do_sb_goby,
  plot_do_sc_goby,
  plot_sav_sc_goby,
  plot_year2_goby,
  ncol = 2,
  scale = 0.9,
  labels = "AUTO"
)

# Second, create the y-axis label as a separate graphical object
y_axis_label <- ggdraw() +
  draw_label(
    "Expected Goby Density (m\u00B2)", # Using unicode for m-squared
    fontface = 'bold',
    size = 16,
    angle = 90,
    x = 0.5, # Center the label within its own narrow column
    y = 0.5
  )

# Third, combine the label and the plot grid side-by-side
# Use rel_widths to give the label a small amount of space
final_plot_with_ylabel <- plot_grid(
  y_axis_label,          # The label object on the left
  final_grid,            # The 4x2 grid of plots on the right
  ncol = 2,
  scale = 0.9,
  rel_widths = c(0.04, 1) # Give the label 4% of the width, and the plots the rest (100%)
)

# Print the final combined plot
print(final_plot_with_ylabel)

# Optional: Save the final combined figure
ggsave(
  filename = "Output/Plots/MarginalEffects_Grid.png",
  plot = final_plot_with_ylabel, # Make sure to save the final object
  width = 16,
  height = 18,
  dpi = 300,
  bg = "white"
)
