# CORRECTED REAL DATA ANALYSIS SCRIPT
# Updated to match simulation code structure with Year fix and consistent total effects calculation
# Updated 2026-01-11

library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(forcats)
library(stringr)
library(cowplot)
library(writexl)

# Fix namespace conflicts
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# ============================================================================
# HELPER FUNCTIONS FOR CAUSAL INFERENCE (UPDATED FROM SIMULATION CODE)
# ============================================================================

# Function to calculate mu_goby deterministically - FULL DAG propagation
calculate_mu_goby_full_dag <- function(draw_row, input_data, 
                                       intervention_var = NULL, 
                                       intervention_value = NULL) {
  
  current_values <- as.list(input_data)
  
  # Helper to get parameter value - handles both dataframe row and named vector
  get_param <- function(param_name) {
    if (is.data.frame(draw_row)) {
      val <- draw_row[[param_name]]
      if (length(val) == 0) {
        stop(paste("Parameter not found:", param_name))
      }
      return(val[1])
    } else {
      return(draw_row[[param_name]])
    }
  }
  
  # Get zone-specific latent variables
  zone_idx <- input_data$Zone
  U_bio_col <- paste0("U_bio[", zone_idx, "]")
  U_phys_col <- paste0("U_phys[", zone_idx, "]")
  
  # If there's an intervention, set it
  if (!is.null(intervention_var)) {
    current_values[[intervention_var]] <- intervention_value
    
    # Update quadratic terms if intervening on variables with quadratics
    if (intervention_var == "Year") {
      current_values$Year_2 <- intervention_value^2
    }
  }
  
  # CONDITIONAL RECALCULATION based on intervention
  # Only recalculate variables DOWNSTREAM of the intervention
  
  # Determine what to recalculate based on intervention
  # IMPORTANT: For baseline (no intervention), don't recalculate ANYTHING - use actual data!
  no_intervention <- is.null(intervention_var)
  
  # Only recalculate if there IS an intervention AND the variable is downstream
  recalc_breach <- !no_intervention && (intervention_var %in% c("Rain"))
  recalc_temp <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind"))
  recalc_do <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp"))
  recalc_sav <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO"))
  recalc_sc <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO", "SAV", "Substrate"))
  recalc_sb <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO", "SAV"))
  
  # BreachDays (affected by Rain)
  if (!is.na(recalc_breach) && recalc_breach && !identical(intervention_var, "BreachDays")) {
    Breach_nu <- get_param("a_BreachDays") + 
      get_param("beta_Rain") * current_values$Rain + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$BreachDays <- Breach_nu
    current_values$BreachDays_2 <- Breach_nu^2
  } else if (identical(intervention_var, "BreachDays")) {
    # Intervention on BreachDays - update quadratic
    current_values$BreachDays_2 <- current_values$BreachDays^2
  }
  
  # Temp (affected by BreachDays, Wind)
  if (!is.na(recalc_temp) && recalc_temp && !identical(intervention_var, "Temp")) {
    Temp_nu <- get_param("a_Temp") + 
      get_param("beta_BreachDays_vec[2]") * current_values$BreachDays + 
      get_param("beta_Wind_vec[2]") * current_values$Wind + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$Temp <- Temp_nu
    current_values$Temp_2 <- Temp_nu^2
  } else if (identical(intervention_var, "Temp")) {
    # Intervention on Temp - update quadratic
    current_values$Temp_2 <- current_values$Temp^2
  }
  
  # DO (affected by Temp, Wind)
  if (!is.na(recalc_do) && recalc_do && !identical(intervention_var, "DO")) {
    DO_nu <- get_param("a_DO") + 
      get_param("beta_Temp_vec[2]") * current_values$Temp + 
      get_param("beta_Wind_vec[1]") * current_values$Wind + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$DO <- DO_nu
  }
  
  # SAV (affected by DO, Temp)
  if (!is.na(recalc_sav) && recalc_sav && !identical(intervention_var, "SAV")) {
    SAV_nu <- get_param("a_SAV") + 
      get_param("beta_DO_vec[4]") * current_values$DO + 
      get_param("beta_Temp_vec[3]") * current_values$Temp + 
      get_param("k_U_bio") * get_param(U_bio_col)
    current_values$SAV <- SAV_nu
    current_values$SAV_2 <- SAV_nu^2
  } else if (identical(intervention_var, "SAV")) {
    # Intervention on SAV - update quadratic
    current_values$SAV_2 <- current_values$SAV^2
  }
  
  # SC_count (affected by Substrate, DO, SAV)
  if (!is.na(recalc_sc) && recalc_sc && !identical(intervention_var, "SC_count")) {
    SC_mu_logit <- get_param("a_SC") + 
      get_param("beta_Substrate_vec[2]") * current_values$Substrate + 
      get_param("beta_DO_vec[3]") * current_values$DO + 
      get_param("beta_SAV_vec[3]") * current_values$SAV + 
      get_param("k_U_bio") * get_param(U_bio_col)
    current_values$SC_count <- plogis(SC_mu_logit)
  }
  
  # SB_count (affected by DO, SAV)
  if (!is.na(recalc_sb) && recalc_sb && !identical(intervention_var, "SB_count")) {
    SB_mu_logit <- get_param("a_SB") + 
      get_param("beta_DO_vec[2]") * current_values$DO + 
      get_param("beta_SAV_vec[2]") * current_values$SAV + 
      get_param("k_U_bio") * get_param(U_bio_col)
    current_values$SB_count <- plogis(SB_mu_logit)
  }
  
  # Final Goby calculation (same whether baseline or intervention)
  a_Goby_col <- paste0("a_Goby[", current_values$Zone, "]")
  mu_goby_val <- get_param(a_Goby_col) +
    get_param("beta_SAV_vec[1]") * current_values$SAV +
    get_param("beta_SAV_2") * current_values$SAV_2 + 
    get_param("beta_DO_vec[1]") * current_values$DO +
    get_param("beta_BreachDays_vec[1]") * current_values$BreachDays +
    get_param("beta_Temp_vec[1]") * current_values$Temp +
    get_param("beta_SC_count") * current_values$SC_count +
    get_param("beta_SB_count") * current_values$SB_count +
    get_param("beta_Substrate_vec[1]") * current_values$Substrate +
    get_param("beta_Micro") * current_values$Micro +
    get_param("beta_Year") * current_values$Year +
    get_param("beta_Year_2") * current_values$Year_2 +
    get_param("beta_Temp_2") * current_values$Temp_2 +
    get_param("beta_BreachDays_2") * current_values$BreachDays_2 +
    get_param("beta_Goby_lag") * current_values$Goby_lag +
    current_values$Area
  
  return(mu_goby_val)
}

# Function to calculate total effects - UPDATED VERSION WITH YEAR FIX
calculate_total_effect <- function(var_name, draws_df, data_df, sds_df, n_draws = 100, n_obs = 100) {
  
  message(paste("  Calculating total effect for:", var_name))
  
  # SPECIAL HANDLING FOR YEAR
  # For Year, we want the marginal effect at mean Year, not the effect of +1 SD intervention
  # This is because Year has a quadratic term, so the effect depends on the baseline value
  if (var_name == "Year") {
    message("    Special handling: Year has quadratic term, reporting marginal effect at mean")
    
    # Select parameter draws
    total_draws_available <- nrow(draws_df)
    if (n_draws > total_draws_available) n_draws <- total_draws_available
    
    if (total_draws_available == 1) {
      selected_draws <- draws_df
      n_draws <- 1
    } else {
      sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
      selected_draws <- draws_df[sample_indices, , drop = FALSE]
    }
    
    # For centered Year (mean = 0), the marginal effect is:
    # d(beta_Year * Year + beta_Year_2 * Year^2) / dYear |_{Year=0}
    # = beta_Year + 2 * beta_Year_2 * 0
    # = beta_Year
    
    effect_dist <- numeric(n_draws)
    for (d in 1:n_draws) {
      draw_row <- selected_draws[d, , drop = FALSE]
      
      # Get parameter value
      get_param <- function(param_name) {
        if (is.data.frame(draw_row)) {
          val <- draw_row[[param_name]]
          if (length(val) == 0) {
            stop(paste("Parameter not found:", param_name))
          }
          return(val[1])
        } else {
          return(draw_row[[param_name]])
        }
      }
      
      beta_year <- get_param("beta_Year")
      
      # Marginal effect at mean Year = beta_Year
      # (Since Year is centered, Year = 0 at the mean)
      effect_dist[d] <- beta_year
    }
    
    return(effect_dist)
  }
  
  # ORIGINAL CODE FOR ALL OTHER VARIABLES
  # Sample observations from the actual data
  n_available <- nrow(data_df)
  if (n_obs > n_available) n_obs <- n_available
  
  obs_indices <- sample(n_available, n_obs, replace = FALSE)
  sampled_obs <- data_df[obs_indices, ]
  
  # Determine intervention magnitude (1 SD)
  if (var_name %in% names(sds_df)) {
    intervention_delta <- sds_df[[1, var_name]]
  } else {
    intervention_delta <- 1  # Default for derived variables
  }
  
  # Select parameter draws
  total_draws_available <- nrow(draws_df)
  if (n_draws > total_draws_available) n_draws <- total_draws_available
  
  if (total_draws_available == 1) {
    selected_draws <- draws_df
    n_draws <- 1
  } else {
    sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
    selected_draws <- draws_df[sample_indices, , drop = FALSE]
  }
  
  # Calculate effects for each draw and each observation
  effects_matrix <- matrix(NA, nrow = n_draws, ncol = n_obs)
  
  for (d in 1:n_draws) {
    draw_row <- selected_draws[d, , drop = FALSE]
    
    for (i in 1:n_obs) {
      obs_row <- sampled_obs[i, ]
      
      # Baseline prediction (using actual data values)
      mu_goby_baseline <- calculate_mu_goby_full_dag(
        draw_row, 
        obs_row, 
        intervention_var = NULL
      )
      
      # Intervened prediction (increase variable by 1 SD)
      mu_goby_intervene <- calculate_mu_goby_full_dag(
        draw_row, 
        obs_row, 
        intervention_var = var_name,
        intervention_value = obs_row[[var_name]] + intervention_delta
      )
      
      effects_matrix[d, i] <- mu_goby_intervene - mu_goby_baseline
    }
  }
  
  # Average across observations for each draw
  effect_dist <- rowMeans(effects_matrix, na.rm = TRUE)
  
  return(effect_dist)
}

# ============================================================================
# PART 1: FIT THE SCM MODEL ON REAL DATA
# ============================================================================

# Make data into a list
list_df <- lapply(as.list(dat), function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(as.vector(x))
  } else {
    return(x)
  }
})

list_df$N <- length(list_df$Goby)
list_df$J <- length(unique(list_df$Zone))

# Run the model
mod.SEM.real <- mod.SEM$sample(
  data = list_df,
  seed = 123,
  chains = 3,
  iter_warmup = 6000,
  iter_sampling = 1000,
  parallel_chains = 3,
  refresh = 100,
  adapt_delta = 0.97,
  max_treedepth = 12,  # Added based on diagnostics
  output_dir = "stan_output",
  save_warmup = FALSE
)

saveRDS(mod.SEM.real, file = "Output/Models/mod.SEM.real.rds")

# Diagnostics
mod.SEM.real$diagnostic_summary()

# Check convergence
summary_output <- mod.SEM.real$summary()

problematic <- summary_output %>% 
  filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400)

cat("\n=== CONVERGENCE CHECK ===\n")
if (nrow(problematic) == 0) {
  cat("✓ All parameters converged successfully (R-hat < 1.01, ESS > 400)\n")
} else {
  cat("⚠ Some parameters have convergence issues:\n")
  print(problematic)
}
cat("========================\n\n")

# Print summary
print(mod.SEM.real$summary(), n = 130)

# ============================================================================
# PART 2: EXTRACT DIRECT EFFECT COEFFICIENTS
# ============================================================================

coefficient_summary.real <- summarise_draws(
  mod.SEM.real,
  "mean",
  ~quantile(.x, probs = c(0.045, 0.955))
)

beta.joint.real <- coefficient_summary.real %>%
  filter(str_starts(variable, "beta_")) %>%
  mutate(variable = str_remove(variable, "^.*?_")) %>%
  # Filter out Rain and Wind - these are path coefficients, not direct effects on Goby
  filter(!variable %in% c("Rain", "Wind"))
beta.joint.real$model <- "SCM"
beta.joint.real$effect_type <- "Direct"

print(beta.joint.real, n = 50)

# Export direct effects
formatted_summary <- mod.SEM.real$summary() %>%
  slice_head(n = 131) %>%
  dplyr::select(-median, -sd, -mad)

write_xlsx(formatted_summary, path = "Output/Tables/stan_model_summary.xlsx")

# ============================================================================
# PART 3: CALCULATE TOTAL EFFECTS FOR REAL DATA
# ============================================================================

cat("\n=== Calculating Total Effects for Real Data ===\n")

# Prepare data
df.data <- dat
df.data$Goby_Density <- df.data$Goby / exp(df.data$Area)

# Calculate standard deviations
numeric_cols <- names(df.data)[sapply(df.data, is.numeric)]
sds <- data.frame(
  purrr::map_dfc(df.data[numeric_cols], ~ sd(.x, na.rm = TRUE))
)
colnames(sds) <- numeric_cols

# Define variables to test
total_effect_vars <- c("Rain", "Wind", "Substrate", "Micro", "Year",
                       "BreachDays", "Temp", "DO", "SAV")

# Extract posterior draws
draws_scm <- mod.SEM.real$draws(format = "df")

# Calculate total effects
scm_total_effects_list <- purrr::map(
  total_effect_vars,
  function(var_name) {
    tryCatch({
      effect_dist <- calculate_total_effect(
        var_name, 
        draws_scm, 
        df.data, 
        sds, 
        n_draws = 100,
        n_obs = 100
      )
      
      tibble(
        variable = paste0(var_name, " (Total)"),
        mean = mean(effect_dist, na.rm = TRUE),
        `5.5%` = quantile(effect_dist, 0.055, na.rm = TRUE),
        `94.5%` = quantile(effect_dist, 0.945, na.rm = TRUE),
        model = "SCM",
        effect_type = "Total"
      )
    }, error = function(e) {
      message(paste("  ERROR:", var_name, "-", e$message))
      return(NULL)
    })
  }
)

scm_total_effects <- bind_rows(scm_total_effects_list[!sapply(scm_total_effects_list, is.null)])

cat("\n=== TOTAL EFFECTS SUMMARY ===\n")
print(scm_total_effects, n = 20)

# ============================================================================
# PART 4: CREATE TOTAL EFFECTS PLOT
# ============================================================================

# Prepare data for plotting
df_total_plot <- scm_total_effects %>% 
  dplyr::rename(lower_5.5 = `5.5%`, upper_94.5 = `94.5%`) %>%
  mutate(variable = factor(variable))

axis_labels_total <- levels(fct_rev(df_total_plot$variable))
label_colors_total <- ifelse(str_detect(axis_labels_total, "Wind|Rain"), "blue", "black")

p.total.real <- df_total_plot %>%
  ggplot(aes(x = fct_rev(variable), y = mean)) +
  geom_pointrange(aes(ymin = lower_5.5, ymax = upper_94.5), 
                  size = 0.7, color = "red3") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_text(colour = label_colors_total),
    plot.title = element_text(face = "bold", size = 20)
  ) +
  labs(title = "Total Causal Effects on Goby (Real Data)", 
       subtitle = "Direct + all mediated pathways (Year at mean)", 
       y = "Total effect (log scale)", x = "")

print(p.total.real)

# Save plot
ggsave("Output/Plots/real_total_effects.png", 
       plot = p.total.real, 
       width = 20, height = 15, units = "cm", dpi = 300)

# ============================================================================
# PART 5: MARGINAL EFFECTS PLOTS (EXISTING CODE - UPDATED WITH NEW FUNCTION)
# ============================================================================

# Marginal Effects Calculation Function using updated calculate_mu_goby_full_dag
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
  
  # Pre-allocate results
  results_df <- tibble(
    var_value = rep(vary_vals, each = n_draws),
    draw_idx = rep(1:n_draws, times = length(vary_vals)),
    mu_goby_pred_log_density = numeric(length(vary_vals) * n_draws)
  )
  
  result_counter <- 1
  
  for (val in vary_vals) {
    for (d in 1:n_draws) {
      draw_row <- selected_draws[d, ]
      intervened_input_data <- baseline_data_row
      
      # Update the variable being varied and its quadratic (if applicable)
      if (vary_var == "Year") {
        intervened_input_data[["Year"]] <- val
        if ("Year_2" %in% names(intervened_input_data)) {
          intervened_input_data[["Year_2"]] <- val^2
        }
      } else if (vary_var == "BreachDays") {
        intervened_input_data[["BreachDays"]] <- val
        if ("BreachDays_2" %in% names(intervened_input_data)) {
          intervened_input_data[["BreachDays_2"]] <- val^2
        }
      } else if (vary_var == "Temp") {
        intervened_input_data[["Temp"]] <- val
        if ("Temp_2" %in% names(intervened_input_data)) {
          intervened_input_data[["Temp_2"]] <- val^2
        }
      } else if (vary_var == "SAV") {
        intervened_input_data[["SAV"]] <- val
        if ("SAV_2" %in% names(intervened_input_data)) {
          intervened_input_data[["SAV_2"]] <- val^2
        }
      } else {
        # For other variables (DO, SC_count, SB_count, Rain), just set the value
        intervened_input_data[[vary_var]] <- val
      }
      
      # Calculate predicted log(count) using the updated function
      mu_goby_pred_log_count <- calculate_mu_goby_full_dag(
        draw_row = draw_row,
        input_data = intervened_input_data,
        intervention_var = vary_var, 
        intervention_value = val      
      )
      
      # Convert log(count) to log(density) by subtracting log(Area)
      log_area_offset <- baseline_data_row$Area
      mu_goby_pred_log_density <- mu_goby_pred_log_count - log_area_offset
      
      # Store result
      results_df$mu_goby_pred_log_density[result_counter] <- mu_goby_pred_log_density
      result_counter <- result_counter + 1
    }
  }
  
  # Summarize across draws
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

# Define Variables and Their Ranges for Plots
sb_vals <- c(0, 1)
sc_vals <- c(0, 1)
breachdays_vals <- seq(min(df.data$BreachDays, na.rm = TRUE), 
                       max(df.data$BreachDays, na.rm = TRUE), 
                       length.out = 50)
rain_vals <- seq(min(df.data$Rain, na.rm = TRUE), 
                 max(df.data$Rain, na.rm = TRUE), 
                 length.out = 50)
do_vals <- seq(min(df.data$DO, na.rm = TRUE), 
               max(df.data$DO, na.rm = TRUE), 
               length.out = 50)
sav_vals <- seq(min(df.data$SAV, na.rm = TRUE), 
                max(df.data$SAV, na.rm = TRUE), 
                length.out = 50)
year_vals <- seq(min(df.data$Year, na.rm = TRUE), 
                 max(df.data$Year, na.rm = TRUE), 
                 length.out = 50)

# Calculate Marginal Effects for Each Variable AND Each Zone
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

relevant_params <- draws_scm

# Calculate with error handling
for (current_zone in all_zones) {
  for (predictor_name in names(predictors_to_plot)) {
    
    tryCatch({
      vals <- predictors_to_plot[[predictor_name]]
      
      me_result <- calculate_marginal_effect_goby_by_zone(
        vary_var = predictor_name,
        vary_vals = vals,
        draws_df = relevant_params,
        data_df = df.data,
        target_zone = current_zone,
        n_draws = 40
      )
      
      me_result$plot_id <- predictor_name
      all_marginal_effects[[length(all_marginal_effects) + 1]] <- me_result
      
    }, error = function(e) {
      message(paste("Error for", predictor_name, "in zone", current_zone, ":", e$message))
    })
  }
}

full_me_df <- do.call(rbind, all_marginal_effects)

# Data Processing for Plots
do_data_sb <- full_me_df %>% 
  filter(plot_id == "DO") %>% 
  mutate(plot_id = "DO_SB_Goby")

do_data_sc <- full_me_df %>% 
  filter(plot_id == "DO") %>% 
  mutate(plot_id = "DO_SC_Goby")

sav_data_sc <- full_me_df %>% 
  filter(plot_id == "SAV") %>% 
  mutate(plot_id = "SAV_SC_Goby")

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
  "SB_count" = "SB → Goby",
  "BreachDays_Goby" = "BreachDays (Total Effect) → Goby",
  "Rain_Goby" = "Rain (Total Effect) → Goby",
  "DO_SB_Goby" = "DO → SB → Goby",
  "DO_SC_Goby" = "DO → SC → Goby",
  "SAV_SC_Goby" = "SAV → SC → Goby",
  "SC_count" = "SC → Goby",
  "Year2_Goby" = "Year & Year² → Goby"
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
  
  if (nrow(plot_data) == 0) {
    message(paste("Warning: No data for plot", plot_identifier))
    return(NULL)
  }
  
  original_predictor_name <- as.character(plot_data$predictor[1])
  is_binary <- original_predictor_name %in% c("SB_count", "SC_count")
  
  p <- ggplot(plot_data, aes(x = var_value, y = mean_goby_density)) +
    geom_point(
      data = raw_df, 
      aes_string(x = original_predictor_name, y = "Goby / exp(Area)"),
      alpha = raw_data_alpha, 
      color = raw_data_color, 
      size = raw_data_size,
      position = if (is_binary) position_jitter(width = jitter_width, height = 0) else "identity"
    ) +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_ribbon(
      aes(ymin = lower_goby_density, ymax = upper_goby_density), 
      fill = "#0072B2", 
      alpha = 0.2
    ) +
    labs(
      title = plot_titles[plot_identifier],
      x = x_axis_labels[original_predictor_name],
      y = NULL
    ) +
    theme_minimal(base_size = 16) +
    facet_wrap(~ Zone, ncol = 3, labeller = labeller(Zone = zone_labels)) +
    coord_cartesian(ylim = c(0, 25))
  
  if (is_binary) {
    p <- p + scale_x_continuous(breaks = c(0, 1))
  }
  
  return(p)
}

# Create plots
plot_sb_goby <- make_marginal_plot("SB_count", full_me_df_processed, df.data)
plot_breachdays_goby <- make_marginal_plot("BreachDays_Goby", full_me_df_processed, df.data)
plot_rain_goby <- make_marginal_plot("Rain_Goby", full_me_df_processed, df.data)
plot_do_sb_goby <- make_marginal_plot("DO_SB_Goby", full_me_df_processed, df.data)
plot_do_sc_goby <- make_marginal_plot("DO_SC_Goby", full_me_df_processed, df.data)
plot_sav_sc_goby <- make_marginal_plot("SAV_SC_Goby", full_me_df_processed, df.data)
plot_sc_goby <- make_marginal_plot("SC_count", full_me_df_processed, df.data)
plot_year2_goby <- make_marginal_plot("Year2_Goby", full_me_df_processed, df.data)

# Arrange plots
final_grid <- plot_grid(
  plot_sb_goby,
  plot_breachdays_goby,
  plot_rain_goby,
  plot_do_sb_goby,
  plot_do_sc_goby,
  plot_sav_sc_goby,
  plot_sc_goby,
  plot_year2_goby,
  ncol = 2,
  labels = "AUTO"
)

# Create y-axis label
y_axis_label <- ggdraw() +
  draw_label(
    "Expected Goby Density (m²)",
    fontface = 'bold',
    size = 16,
    angle = 90,
    x = 0.5,
    y = 0.5
  )

# Combine
final_plot_with_ylabel <- plot_grid(
  y_axis_label,
  final_grid,
  ncol = 2,
  scale = 0.9,
  rel_widths = c(0.04, 1)
)

print(final_plot_with_ylabel)

# Save
ggsave(
  filename = "Output/Plots/MarginalEffects_Grid_Corrected.png",
  plot = final_plot_with_ylabel,
  width = 16,
  height = 18,
  dpi = 300
)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Direct effects: See beta.joint.real\n")
cat("Total effects: See scm_total_effects\n")
cat("Plots saved to Output/Plots/\n")
cat("Tables saved to Output/Tables/\n\n")



## plot attempt #2
# IMPROVED TOTAL EFFECTS VISUALIZATION
# Addresses the misleading coefficient plot issue
# Updated 2026-01-12

library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(forcats)
library(stringr)

# ============================================================================
# SOLUTION 1: Calculate ACTUAL total effects accounting for quadratics
# ============================================================================

calculate_realistic_total_effect <- function(var_name, draws_df, data_df, sds_df, n_draws = 100, n_obs = 100) {
  
  message(paste("  Calculating realistic total effect for:", var_name))
  
  # For variables with quadratic terms, we need to show the effect
  # over a realistic range, not just at one point
  
  if (var_name == "Year") {
    # For Year: Calculate effect at mean +/- 1 SD to capture the nonlinearity
    
    total_draws_available <- nrow(draws_df)
    if (n_draws > total_draws_available) n_draws <- total_draws_available
    
    if (total_draws_available == 1) {
      selected_draws <- draws_df
      n_draws <- 1
    } else {
      sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
      selected_draws <- draws_df[sample_indices, , drop = FALSE]
    }
    
    # Calculate the average marginal effect across the observed range
    # This captures the nonlinearity better than just the effect at mean
    year_range <- seq(-1, 1, length.out = 10)  # -1SD to +1SD
    
    effect_dist <- numeric(n_draws)
    for (d in 1:n_draws) {
      draw_row <- selected_draws[d, , drop = FALSE]
      
      get_param <- function(param_name) {
        if (is.data.frame(draw_row)) {
          val <- draw_row[[param_name]]
          if (length(val) == 0) stop(paste("Parameter not found:", param_name))
          return(val[1])
        } else {
          return(draw_row[[param_name]])
        }
      }
      
      beta_year <- get_param("beta_Year")
      beta_year2 <- get_param("beta_Year_2")
      
      # Calculate marginal effect at multiple points and average
      marginal_effects <- beta_year + 2 * beta_year2 * year_range
      effect_dist[d] <- mean(marginal_effects)
    }
    
    return(effect_dist)
  }
  
  # For other variables, use the standard calculation
  # (Your existing calculate_total_effect function for non-Year variables)
  n_available <- nrow(data_df)
  if (n_obs > n_available) n_obs <- n_available
  
  obs_indices <- sample(n_available, n_obs, replace = FALSE)
  sampled_obs <- data_df[obs_indices, ]
  
  if (var_name %in% names(sds_df)) {
    intervention_delta <- sds_df[[1, var_name]]
  } else {
    intervention_delta <- 1
  }
  
  total_draws_available <- nrow(draws_df)
  if (n_draws > total_draws_available) n_draws <- total_draws_available
  
  if (total_draws_available == 1) {
    selected_draws <- draws_df
    n_draws <- 1
  } else {
    sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
    selected_draws <- draws_df[sample_indices, , drop = FALSE]
  }
  
  effects_matrix <- matrix(NA, nrow = n_draws, ncol = n_obs)
  
  for (d in 1:n_draws) {
    draw_row <- selected_draws[d, , drop = FALSE]
    
    for (i in 1:n_obs) {
      obs_row <- sampled_obs[i, ]
      
      mu_goby_baseline <- calculate_mu_goby_full_dag(
        draw_row, obs_row, intervention_var = NULL
      )
      
      mu_goby_intervene <- calculate_mu_goby_full_dag(
        draw_row, obs_row, 
        intervention_var = var_name,
        intervention_value = obs_row[[var_name]] + intervention_delta
      )
      
      effects_matrix[d, i] <- mu_goby_intervene - mu_goby_baseline
    }
  }
  
  effect_dist <- rowMeans(effects_matrix, na.rm = TRUE)
  return(effect_dist)
}

# ============================================================================
# SOLUTION 2: Decompose effects into Direct + Indirect components
# ============================================================================

decompose_total_effect <- function(var_name, draws_df, data_df, sds_df, n_draws = 100, n_obs = 100) {
  
  message(paste("  Decomposing effect for:", var_name))
  
  # Get total effect
  total_effect <- calculate_total_effect(var_name, draws_df, data_df, sds_df, n_draws, n_obs)
  
  # Get direct effect (just the beta coefficient if it's in the Goby equation)
  direct_effect <- numeric(length(total_effect))
  
  total_draws_available <- nrow(draws_df)
  actual_n_draws <- length(total_effect)  # Use the actual number from total_effect
  
  if (total_draws_available == 1) {
    selected_draws <- draws_df
  } else {
    # Use the SAME draws that were used for total_effect calculation
    sample_indices <- sample(total_draws_available, actual_n_draws, replace = FALSE)
    selected_draws <- draws_df[sample_indices, , drop = FALSE]
  }
  
  # CRITICAL: Define which variables are ACTUALLY in the Goby equation
  # Based on your Stan model, only these variables directly affect Goby:
  vars_in_goby_equation <- c(
    "SAV", "SAV_2", "DO", "BreachDays", "BreachDays_2", "Temp", "Temp_2",
    "SC_count", "SB_count", "Substrate", "Micro", "Year", "Year_2", "Goby_lag"
  )
  
  # Check if this variable is in the Goby equation
  if (var_name %in% vars_in_goby_equation) {
    
    # Try different possible parameter names
    possible_names <- c(
      paste0("beta_", var_name, "_vec[1]"),  # For variables like SAV, DO, etc.
      paste0("beta_", var_name),              # For variables like Micro, Year
      paste0("beta_", var_name, "_vec.1.")    # Alternative format
    )
    
    beta_name <- NULL
    for (name in possible_names) {
      if (name %in% names(selected_draws)) {
        beta_name <- name
        break
      }
    }
    
    if (!is.null(beta_name)) {
      for (d in 1:actual_n_draws) {
        direct_effect[d] <- selected_draws[d, beta_name][[1]]
      }
      message(paste("    Found direct effect:", beta_name))
    } else {
      # Couldn't find the parameter - this is a problem
      warning(paste("Could not find direct effect parameter for", var_name))
      direct_effect <- rep(0, actual_n_draws)
    }
    
  } else {
    # Variable NOT in Goby equation (Rain, Wind)
    message(paste("    No direct effect (not in Goby equation)"))
    direct_effect <- rep(0, actual_n_draws)
  }
  
  # Indirect effect = Total - Direct
  indirect_effect <- total_effect - direct_effect
  
  return(list(
    total = total_effect,
    direct = direct_effect,
    indirect = indirect_effect
  ))
}

# ============================================================================
# SOLUTION 3: Create an enhanced forest plot with decomposition
# ============================================================================

create_enhanced_effects_plot <- function(draws_scm, df.data, sds, total_effect_vars) {
  
  # Calculate decomposed effects for each variable
  all_effects <- list()
  
  for (var in total_effect_vars) {
    
    tryCatch({
      # Calculate decomposed effects
      effects <- decompose_total_effect(var, draws_scm, df.data, sds, n_draws = 100, n_obs = 100)
      
      # Create dataframe with all components
      df_var <- tibble(
        variable = var,
        effect_type = c("Direct", "Indirect", "Total"),
        mean = c(mean(effects$direct), mean(effects$indirect), mean(effects$total)),
        lower = c(quantile(effects$direct, 0.055), 
                  quantile(effects$indirect, 0.055), 
                  quantile(effects$total, 0.055)),
        upper = c(quantile(effects$direct, 0.945), 
                  quantile(effects$indirect, 0.945), 
                  quantile(effects$total, 0.945))
      )
      
      all_effects[[var]] <- df_var
      
    }, error = function(e) {
      message(paste("Error for", var, ":", e$message))
    })
  }
  
  # Combine all effects
  plot_data <- bind_rows(all_effects)
  
  # Create the plot
  p <- plot_data %>%
    mutate(
      variable = factor(variable, levels = rev(total_effect_vars)),
      effect_type = factor(effect_type, levels = c("Direct", "Indirect", "Total"))
    ) %>%
    ggplot(aes(x = variable, y = mean, color = effect_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(
      aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 0.5),
      size = 0.6
    ) +
    coord_flip() +
    scale_color_manual(
      values = c("Direct" = "#E69F00", "Indirect" = "#56B4E9", "Total" = "#D55E00"),
      name = "Effect Type"
    ) +
    labs(
      title = "Decomposed Causal Effects on Goby",
      subtitle = "Direct effects + Indirect (mediated) effects = Total effects",
      x = "",
      y = "Effect on log(Goby density)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16)
    )
  
  return(list(plot = p, data = plot_data))
}

# ============================================================================
# SOLUTION 4: Alternative - Show effect at different points for Year
# ============================================================================

create_year_effect_breakdown <- function(draws_scm) {
  
  message("Creating Year effect breakdown...")
  
  # Extract Year coefficients
  beta_year_draws <- draws_scm$beta_Year
  beta_year2_draws <- draws_scm$beta_Year_2
  
  # Calculate marginal effects at different points
  year_points <- c(-1, -0.5, 0, 0.5, 1)  # -1SD, -0.5SD, mean, +0.5SD, +1SD
  
  year_effects <- tibble(
    year_value = rep(year_points, each = length(beta_year_draws)),
    draw = rep(1:length(beta_year_draws), times = length(year_points)),
    marginal_effect = rep(beta_year_draws, times = length(year_points)) + 
      2 * rep(beta_year2_draws, times = length(year_points)) * 
      rep(year_points, each = length(beta_year_draws))
  )
  
  year_summary <- year_effects %>%
    group_by(year_value) %>%
    summarise(
      mean = mean(marginal_effect),
      lower = quantile(marginal_effect, 0.055),
      upper = quantile(marginal_effect, 0.945),
      .groups = "drop"
    )
  
  # Create plot
  p_year <- ggplot(year_summary, aes(x = year_value, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#D55E00") +
    geom_line(color = "#D55E00", size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    labs(
      title = "Year Effect on Goby (Marginal Effect)",
      subtitle = "Shows how Year effect varies across the observed range",
      x = "Year (standardized: 0 = mean)",
      y = "Marginal effect on log(Goby density)"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  return(list(plot = p_year, data = year_summary))
}

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# Using your real data analysis:

# 1. Create decomposed effects plot
cat("\n=== Creating Enhanced Effects Plot ===\n")
enhanced_results <- create_enhanced_effects_plot(
  draws_scm = draws_scm,
  df.data = df.data,
  sds = sds,
  total_effect_vars = c("Rain", "Wind", "Substrate", "Micro", "Year",
                        "BreachDays", "Temp", "DO", "SAV")
)

print(enhanced_results$plot)

ggsave("Output/Plots/decomposed_effects.png", 
       plot = enhanced_results$plot,
       width = 10, height = 8, dpi = 300)

# 2. Create Year-specific breakdown
cat("\n=== Creating Year Effect Breakdown ===\n")
year_results <- create_year_effect_breakdown(draws_scm)

print(year_results$plot)

ggsave("Output/Plots/year_effect_breakdown.png",
       plot = year_results$plot,
       width = 8, height = 6, dpi = 300)

# 3. Print summary table
cat("\n=== Decomposed Effects Summary ===\n")
enhanced_results$data %>%
  arrange(variable, effect_type) %>%
  print(n = 50)

# 4. Identify which variables have strong indirect effects
cat("\n=== Variables with Substantial Indirect Effects ===\n")
enhanced_results$data %>%
  filter(effect_type == "Indirect") %>%
  mutate(abs_indirect = abs(mean)) %>%
  arrange(desc(abs_indirect)) %>%
  select(variable, mean, lower, upper) %>%
  print()

# ============================================================================
# INTERPRETATION GUIDE
# ============================================================================

cat("\n=== HOW TO INTERPRET THE DECOMPOSED PLOT ===\n")
cat("
1. DIRECT effects (orange): 
   - The coefficient in the Goby equation
   - What you see in standard regression

2. INDIRECT effects (blue):
   - Effects through mediated pathways
   - Rain → BreachDays → Temp → ... → Goby
   - Often larger than direct effects!

3. TOTAL effects (red):
   - What actually matters for prediction/intervention
   - Direct + Indirect = Total

4. For variables like Rain and Wind:
   - Direct = 0 (not in Goby equation)
   - All effect is Indirect (through pathways)
   - This is why they look 'small' in coefficient plot but show patterns in marginal plots

5. For Year:
   - The effect shown is averaged across the observed range
   - See separate Year breakdown plot for full nonlinearity
   - The quadratic term creates the inverted-U shape you see in marginal plots
")


##---

#zone coefficient diagnostic plot

# ZONE-SPECIFIC TOTAL EFFECTS DIAGNOSTIC
# Investigates why Rain/Wind show patterns in marginal plots but small total effects
# Updated 2026-01-12

library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)
library(forcats)

# ============================================================================
# FUNCTION: Calculate total effects separately for each zone
# ============================================================================

calculate_total_effect_by_zone <- function(var_name, draws_df, data_df, sds_df, 
                                           target_zone, n_draws = 100, n_obs = 50) {
  
  message(paste("  Calculating total effect for:", var_name, "in Zone", target_zone))
  
  # SPECIAL HANDLING FOR YEAR (same as before)
  if (var_name == "Year") {
    total_draws_available <- nrow(draws_df)
    if (n_draws > total_draws_available) n_draws <- total_draws_available
    
    if (total_draws_available == 1) {
      selected_draws <- draws_df
      n_draws <- 1
    } else {
      sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
      selected_draws <- draws_df[sample_indices, , drop = FALSE]
    }
    
    effect_dist <- numeric(n_draws)
    for (d in 1:n_draws) {
      draw_row <- selected_draws[d, , drop = FALSE]
      
      get_param <- function(param_name) {
        if (is.data.frame(draw_row)) {
          val <- draw_row[[param_name]]
          if (length(val) == 0) stop(paste("Parameter not found:", param_name))
          return(val[1])
        } else {
          return(draw_row[[param_name]])
        }
      }
      
      beta_year <- get_param("beta_Year")
      effect_dist[d] <- beta_year
    }
    
    return(effect_dist)
  }
  
  # FILTER DATA TO ONLY THIS ZONE
  zone_data <- data_df %>% filter(Zone == target_zone)
  
  n_available <- nrow(zone_data)
  if (n_available == 0) {
    warning(paste("No data for zone", target_zone))
    return(rep(NA, n_draws))
  }
  
  if (n_obs > n_available) n_obs <- n_available
  
  obs_indices <- sample(n_available, n_obs, replace = FALSE)
  sampled_obs <- zone_data[obs_indices, ]
  
  # Determine intervention magnitude (1 SD from FULL dataset, not just this zone)
  if (var_name %in% names(sds_df)) {
    intervention_delta <- sds_df[[1, var_name]]
  } else {
    intervention_delta <- 1
  }
  
  # Select parameter draws
  total_draws_available <- nrow(draws_df)
  if (n_draws > total_draws_available) n_draws <- total_draws_available
  
  if (total_draws_available == 1) {
    selected_draws <- draws_df
    n_draws <- 1
  } else {
    sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
    selected_draws <- draws_df[sample_indices, , drop = FALSE]
  }
  
  # Calculate effects for each draw and each observation
  effects_matrix <- matrix(NA, nrow = n_draws, ncol = n_obs)
  
  for (d in 1:n_draws) {
    draw_row <- selected_draws[d, , drop = FALSE]
    
    for (i in 1:n_obs) {
      obs_row <- sampled_obs[i, ]
      
      # Baseline prediction
      mu_goby_baseline <- calculate_mu_goby_full_dag(
        draw_row, obs_row, intervention_var = NULL
      )
      
      # Intervened prediction
      mu_goby_intervene <- calculate_mu_goby_full_dag(
        draw_row, obs_row, 
        intervention_var = var_name,
        intervention_value = obs_row[[var_name]] + intervention_delta
      )
      
      effects_matrix[d, i] <- mu_goby_intervene - mu_goby_baseline
    }
  }
  
  # Average across observations for each draw
  effect_dist <- rowMeans(effects_matrix, na.rm = TRUE)
  return(effect_dist)
}

# ============================================================================
# FUNCTION: Calculate zone-specific effects for all variables
# ============================================================================

calculate_all_zone_specific_effects <- function(draws_scm, df.data, sds, 
                                                total_effect_vars,
                                                n_draws = 100, n_obs = 50) {
  
  all_zones <- sort(unique(df.data$Zone))
  all_results <- list()
  
  for (zone in all_zones) {
    cat("\n=== Processing Zone", zone, "===\n")
    
    for (var in total_effect_vars) {
      tryCatch({
        effect_dist <- calculate_total_effect_by_zone(
          var_name = var,
          draws_df = draws_scm,
          data_df = df.data,
          sds_df = sds,
          target_zone = zone,
          n_draws = n_draws,
          n_obs = n_obs
        )
        
        result_df <- tibble(
          variable = var,
          zone = zone,
          mean = mean(effect_dist, na.rm = TRUE),
          median = median(effect_dist, na.rm = TRUE),
          lower = quantile(effect_dist, 0.055, na.rm = TRUE),
          upper = quantile(effect_dist, 0.945, na.rm = TRUE),
          sd = sd(effect_dist, na.rm = TRUE)
        )
        
        all_results[[paste(var, zone, sep = "_")]] <- result_df
        
      }, error = function(e) {
        message(paste("  ERROR:", var, "Zone", zone, "-", e$message))
      })
    }
  }
  
  return(bind_rows(all_results))
}

# ============================================================================
# FUNCTION: Create zone-specific effects plot
# ============================================================================

create_zone_specific_plot <- function(zone_effects_df, highlight_vars = c("Rain", "Wind")) {
  
  zone_labels <- c("1" = "East", "2" = "Northwest", "3" = "West")
  
  # Filter to highlight variables
  plot_data <- zone_effects_df %>%
    filter(variable %in% highlight_vars) %>%
    mutate(
      zone_label = factor(zone, levels = c(1, 2, 3), labels = zone_labels),
      variable = factor(variable, levels = rev(highlight_vars))
    )
  
  p <- ggplot(plot_data, aes(x = variable, y = mean, color = zone_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(
      aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 0.5),
      size = 0.8
    ) +
    coord_flip() +
    scale_color_manual(
      values = c("East" = "#E69F00", "Northwest" = "#56B4E9", "West" = "#009E73"),
      name = "Zone"
    ) +
    labs(
      title = "Zone-Specific Total Effects on Goby",
      subtitle = "Rain and Wind effects by zone",
      x = "",
      y = "Total effect on log(Goby density)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16)
    )
  
  return(p)
}

# ============================================================================
# FUNCTION: Create heatmap of all zone-specific effects
# ============================================================================

create_zone_effects_heatmap <- function(zone_effects_df) {
  
  zone_labels <- c("1" = "East", "2" = "Northwest", "3" = "West")
  
  plot_data <- zone_effects_df %>%
    mutate(
      zone_label = factor(zone, levels = c(1, 2, 3), labels = zone_labels),
      significant = sign(lower) == sign(upper)  # CI doesn't cross zero
    )
  
  p <- ggplot(plot_data, aes(x = zone_label, y = variable, fill = mean)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 3) +
    geom_point(
      data = plot_data %>% filter(!significant),
      aes(x = zone_label, y = variable),
      shape = 4, size = 3, color = "white", stroke = 1.5
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      name = "Effect",
      limits = c(-max(abs(plot_data$mean)), max(abs(plot_data$mean)))
    ) +
    labs(
      title = "Zone-Specific Total Effects Heatmap",
      subtitle = "X marks indicate effect not significantly different from zero",
      x = "Zone",
      y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  return(p)
}

# ============================================================================
# FUNCTION: Analyze why effects are small (decompose by baseline values)
# ============================================================================

analyze_baseline_values_by_zone <- function(df.data, vars_of_interest = c("Rain", "Wind", "BreachDays", "Temp", "DO")) {
  
  zone_labels <- c("1" = "East", "2" = "Northwest", "3" = "West")
  
  baseline_summary <- df.data %>%
    group_by(Zone) %>%
    summarise(
      across(all_of(vars_of_interest), 
             list(mean = ~mean(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE),
                  min = ~min(.x, na.rm = TRUE),
                  max = ~max(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(zone_label = zone_labels[as.character(Zone)])
  
  cat("\n=== BASELINE VALUES BY ZONE ===\n")
  print(baseline_summary, width = Inf)
  
  # Check if Rain varies much by zone
  rain_by_zone <- df.data %>%
    group_by(Zone) %>%
    summarise(
      Rain_mean = mean(Rain, na.rm = TRUE),
      Rain_sd = sd(Rain, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n=== RAIN DISTRIBUTION BY ZONE ===\n")
  print(rain_by_zone)
  
  # Check correlations by zone
  cat("\n=== CORRELATIONS BY ZONE ===\n")
  for (z in sort(unique(df.data$Zone))) {
    cat(paste("\nZone", z, "-", zone_labels[as.character(z)], ":\n"))
    zone_data <- df.data %>% filter(Zone == z)
    
    cor_rain_breach <- cor(zone_data$Rain, zone_data$BreachDays, use = "complete.obs")
    cor_rain_goby <- cor(zone_data$Rain, zone_data$Goby_Density, use = "complete.obs")
    
    cat(sprintf("  Rain ~ BreachDays: %.3f\n", cor_rain_breach))
    cat(sprintf("  Rain ~ Goby_Density: %.3f\n", cor_rain_goby))
  }
  
  return(baseline_summary)
}

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

cat("\n=== ZONE-SPECIFIC EFFECTS DIAGNOSTIC ===\n")
cat("This will take a few minutes...\n\n")

# 1. Calculate zone-specific effects for all variables
zone_effects <- calculate_all_zone_specific_effects(
  draws_scm = draws_scm,
  df.data = df.data,
  sds = sds,
  total_effect_vars = c("Rain", "Wind", "Substrate", "Micro", "Year",
                        "BreachDays", "Temp", "DO", "SAV"),
  n_draws = 100,
  n_obs = 50
)

# 2. Print summary
cat("\n=== ZONE-SPECIFIC EFFECTS SUMMARY ===\n")
print(zone_effects, n = 50)

# 3. Focus on Rain and Wind
cat("\n=== RAIN AND WIND BY ZONE ===\n")
rain_wind_by_zone <- zone_effects %>%
  filter(variable %in% c("Rain", "Wind")) %>%
  arrange(variable, zone)

print(rain_wind_by_zone)

# 4. Create zone-specific plot for Rain and Wind
p_zone_rain_wind <- create_zone_specific_plot(zone_effects, highlight_vars = c("Rain", "Wind"))
print(p_zone_rain_wind)

ggsave("Output/Plots/zone_specific_rain_wind.png",
       plot = p_zone_rain_wind,
       width = 8, height = 5, dpi = 300)

# 5. Create heatmap of all effects
p_heatmap <- create_zone_effects_heatmap(zone_effects)
print(p_heatmap)

ggsave("Output/Plots/zone_effects_heatmap.png",
       plot = p_heatmap,
       width = 10, height = 7, dpi = 300)

# 6. Analyze baseline values by zone
baseline_analysis <- analyze_baseline_values_by_zone(
  df.data = df.data,
  vars_of_interest = c("Rain", "Wind", "BreachDays", "Temp", "DO", "SAV")
)

# 7. Compare zone-specific vs. overall effects
cat("\n=== COMPARISON: ZONE-SPECIFIC VS OVERALL ===\n")

overall_effects <- scm_total_effects %>%
  mutate(variable = str_remove(variable, " \\(Total\\)")) %>%
  select(variable, overall_mean = mean, overall_lower = `5.5%`, overall_upper = `94.5%`)

comparison <- zone_effects %>%
  left_join(overall_effects, by = "variable") %>%
  mutate(
    diff_from_overall = mean - overall_mean,
    zone_stronger = abs(mean) > abs(overall_mean)
  )

cat("\nVariables with stronger effects in at least one zone:\n")
comparison %>%
  filter(zone_stronger) %>%
  arrange(variable, desc(abs(mean))) %>%
  select(variable, zone, mean, overall_mean, diff_from_overall) %>%
  print(n = 50)

# 8. Statistical test: Do zones differ significantly?
cat("\n=== ZONE DIFFERENCES TEST ===\n")

for (var in c("Rain", "Wind", "BreachDays", "Temp", "DO", "SAV")) {
  var_data <- zone_effects %>% filter(variable == var)
  
  # Check if confidence intervals overlap
  overlap_12 <- max(var_data$lower[var_data$zone == 1]) >= min(var_data$upper[var_data$zone == 2]) &&
    min(var_data$upper[var_data$zone == 1]) >= max(var_data$lower[var_data$zone == 2])
  overlap_13 <- max(var_data$lower[var_data$zone == 1]) >= min(var_data$upper[var_data$zone == 3]) &&
    min(var_data$upper[var_data$zone == 1]) >= max(var_data$lower[var_data$zone == 3])
  overlap_23 <- max(var_data$lower[var_data$zone == 2]) >= min(var_data$upper[var_data$zone == 3]) &&
    min(var_data$upper[var_data$zone == 2]) >= max(var_data$lower[var_data$zone == 3])
  
  all_overlap <- overlap_12 & overlap_13 & overlap_23
  
  cat(sprintf("\n%s: Zones %s\n", var, 
              if(all_overlap) "do NOT differ significantly" else "MAY differ"))
  
  if (!all_overlap) {
    cat("  Zone means:\n")
    for (z in 1:3) {
      zone_mean <- var_data$mean[var_data$zone == z]
      zone_ci <- paste0("[", 
                        sprintf("%.3f", var_data$lower[var_data$zone == z]), ", ",
                        sprintf("%.3f", var_data$upper[var_data$zone == z]), "]")
      cat(sprintf("    Zone %d: %.3f %s\n", z, zone_mean, zone_ci))
    }
  }
}

# 9. Create summary table for export
summary_table <- zone_effects %>%
  select(variable, zone, mean, lower, upper) %>%
  pivot_wider(
    names_from = zone,
    values_from = c(mean, lower, upper),
    names_sep = "_Zone"
  )

write.csv(summary_table, "Output/Tables/zone_specific_effects.csv", row.names = FALSE)

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
cat("Key findings:\n")
cat("1. Check if Rain/Wind effects differ by zone (plot and heatmap)\n")
cat("2. Check if baseline values differ by zone (printed above)\n")
cat("3. Check if averaging across zones is masking zone-specific effects\n")
cat("4. Tables saved to Output/Tables/zone_specific_effects.csv\n\n")




## old below

# CORRECTED REAL DATA ANALYSIS SCRIPT
# Fixed to use U_bio and U_phys correctly

#run the real data

# #make data into a list
# list_df <- as.list(dat)
# 
# # Need to Convert dataframe to list and flatten any 1-column matrices
# list_df <- lapply(as.list(dat), function(x) {
#   if (is.matrix(x) && ncol(x) == 1) {
#     return(as.vector(x))
#   } else {
#     return(x)
#   }
# })
# 
# #need an N and J for mod.SEM modelindexing
# list_df$N <- length(list_df$Goby)
# 
# # Calculate J from the number of unique zones
# list_df$J <- length(unique(list_df$Zone))
# 
# #run the model
# library(cmdstanr)
# 
# mod.SEM.real <- mod.SEM$sample(
#   data = list_df,
#   seed = 123,
#   chains = 3,
#   iter_warmup = 6000,
#   iter_sampling = 1000,
#   parallel_chains = 3,
#   refresh = 100,
#   adapt_delta = 0.97,
#   output_dir = "stan_output",
#   save_warmup = FALSE
# )
# 
# saveRDS(mod.SEM.real, file = "Output/Models/mod.SEM.real.rds")
# 
# mod.SEM.real$diagnostic_summary()
# 
# #more diagnostics: 
# # Check which parameters are problematic
# mod.SEM.real$cmdstan_diagnose()
# 
# # Look at pairs plot for correlations
# library(bayesplot)
# mcmc_pairs(mod.SEM.real$draws(), pars = c("beta_Year", "beta_Year_2", "k_U_bio", "k_U_phys"))
# 
# 
# # Even with chain 3's warnings, check convergence
# summary_output <- mod.SEM.real$summary()
# 
# # Look for R-hat problems
# problematic <- summary_output %>% 
#   filter(rhat > 1.01 | ess_bulk < 400 | ess_tail < 400)
# 
# print(problematic)
# 
# # If this is empty or minimal, you're actually fine!
# 
# 
# # for saving and seeing
# print(mod.SEM.real$summary(), n = 1000)
# 
# ## get summary
# coefficient_summary.real <- summarise_draws(
#   mod.SEM.real,
#   "mean",
#   ~quantile(.x, probs = c(0.045, 0.955))
# )
# 
# beta.joint.real <- coefficient_summary.real %>%
#   filter(str_starts(variable, "beta_")) %>%
#   mutate(variable = str_remove(variable, "^.*?_"))
# beta.joint.real$model <- "Causal SEM"
# print(beta.joint.real, n = 50) 
# 
# # Export results
# library(dplyr)
# library(writexl)
# 
# SEM.real.coefficients <- mod.SEM.real$summary()
# 
# formatted_summary <- SEM.real.coefficients %>%
#   slice_head(n = 131) %>%
#   dplyr::select(-median, -sd, -mad)
# 
# output_filename <- "Output/Tables/stan_model_summary.xlsx"
# write_xlsx(formatted_summary, path = output_filename)
# 
# ### Plots ----------
# 
# df.data <- dat
# df.data$Goby_Density <- df.data$Goby / exp(df.data$Area)
# 
# ## and now the effects plots on goby density: 
# 
# # --- 0. Setup ---
# library(cmdstanr)
# library(posterior)
# library(dplyr)
# library(ggplot2)
# library(purrr)
# library(tibble)
# library(forcats)
# library(stringr)
# library(cowplot)
# 
# # --- 1. Extract Posterior Draws ---
# draws <- mod.SEM.real$draws(format = "df")
# relevant_params <- draws
# 
# # --- 2. CORRECTED Deterministic Linear Predictor Function ---
# calculate_mu_goby_deterministic <- function(draw_row, input_data, 
#                                             initial_intervention_var = NULL, 
#                                             initial_intervention_value = NULL) {
#   
#   current_values <- as.list(input_data)
#   if (!is.null(initial_intervention_var)) {
#     current_values[[initial_intervention_var]] <- initial_intervention_value
#   }
#   
#   get_val <- function(var_name) {
#     return(current_values[[var_name]])
#   }
#   
#   is_intervened <- function(var) {
#     !is.null(initial_intervention_var) && initial_intervention_var == var
#   }
#   
#   # CRITICAL: Get zone-specific latent variables
#   zone_idx <- input_data$Zone
#   U_bio_col <- paste0("U_bio[", zone_idx, "]")
#   U_phys_col <- paste0("U_phys[", zone_idx, "]")
#   
#   # BreachDays (PHYSICAL variable - use U_phys)
#   if (!is_intervened("BreachDays")) {
#     Breach_nu <- draw_row[["a_BreachDays"]] + 
#       draw_row[["beta_Rain"]] * get_val("Rain") + 
#       draw_row[["k_U_phys"]] * draw_row[[U_phys_col]]  # FIXED
#     current_values$BreachDays <- Breach_nu
#     current_values$BreachDays_2 <- Breach_nu^2  # Update quadratic
#   }
#   
#   # Temp (PHYSICAL variable - use U_phys)
#   if (!is_intervened("Temp")) {
#     Temp_nu <- draw_row[["a_Temp"]] + 
#       draw_row[["beta_BreachDays_vec[2]"]] * get_val("BreachDays") + 
#       draw_row[["beta_Wind_vec[2]"]] * get_val("Wind") + 
#       draw_row[["k_U_phys"]] * draw_row[[U_phys_col]]  # FIXED
#     current_values$Temp <- Temp_nu
#     current_values$Temp_2 <- Temp_nu^2  # Update quadratic
#   }
#   
#   # DO (PHYSICAL variable - use U_phys)
#   if (!is_intervened("DO")) {
#     DO_nu <- draw_row[["a_DO"]] + 
#       draw_row[["beta_Temp_vec[2]"]] * get_val("Temp") + 
#       draw_row[["beta_Wind_vec[1]"]] * get_val("Wind") + 
#       draw_row[["k_U_phys"]] * draw_row[[U_phys_col]]  # FIXED
#     current_values$DO <- DO_nu
#   }
#   
#   # SAV (BIOLOGICAL variable - use U_bio)
#   if (!is_intervened("SAV")) {
#     SAV_nu <- draw_row[["a_SAV"]] + 
#       draw_row[["beta_DO_vec[4]"]] * get_val("DO") + 
#       draw_row[["beta_Temp_vec[3]"]] * get_val("Temp") + 
#       draw_row[["k_U_bio"]] * draw_row[[U_bio_col]]  # CORRECT
#     current_values$SAV <- SAV_nu
#     current_values$SAV_2 <- SAV_nu^2  # Update quadratic
#   }
#   
#   # SC_count (BIOLOGICAL variable - use U_bio)
#   if (!is_intervened("SC_count")) {
#     SC_mu_logit <- draw_row[["a_SC"]] + 
#       draw_row[["beta_Substrate_vec[2]"]] * get_val("Substrate") + 
#       draw_row[["beta_DO_vec[3]"]] * get_val("DO") + 
#       draw_row[["beta_SAV_vec[3]"]] * get_val("SAV") + 
#       draw_row[["k_U_bio"]] * draw_row[[U_bio_col]]  # CORRECT
#     current_values$SC_count <- plogis(SC_mu_logit)
#   }
#   
#   # SB_count (BIOLOGICAL variable - use U_bio)
#   if (!is_intervened("SB_count")) {
#     SB_mu_logit <- draw_row[["a_SB"]] + 
#       draw_row[["beta_DO_vec[2]"]] * get_val("DO") + 
#       draw_row[["beta_SAV_vec[2]"]] * get_val("SAV") + 
#       draw_row[["k_U_bio"]] * draw_row[[U_bio_col]]  # CORRECT
#     current_values$SB_count <- plogis(SB_mu_logit)
#   }
#   
#   # Final Goby model (does NOT include U directly - correct)
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
# 
# # --- 3. Marginal Effects Calculation Function (CORRECTED) ---
# calculate_marginal_effect_goby_by_zone <- function(
#     vary_var,                   
#     vary_vals,                  
#     draws_df,                   
#     data_df,                    
#     target_zone,                
#     n_draws = 40             
# ) {
#   
#   message(paste("Calculating marginal effect for:", vary_var, "in Zone", target_zone))
#   
#   baseline_data_row <- data_df %>%
#     dplyr::summarise(dplyr::across(dplyr::everything(), ~median(.x, na.rm = TRUE))) %>% 
#     dplyr::mutate(Zone = as.integer(target_zone)) %>%
#     as.data.frame()
#   
#   total_draws_available <- nrow(draws_df)
#   if (n_draws > total_draws_available) {
#     n_draws <- total_draws_available
#   }
#   sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
#   selected_draws <- draws_df[sample_indices, ]
#   
#   # Pre-allocate results
#   results_df <- tibble(
#     var_value = rep(vary_vals, each = n_draws),
#     draw_idx = rep(1:n_draws, times = length(vary_vals)),
#     mu_goby_pred_log_density = numeric(length(vary_vals) * n_draws)
#   )
#   
#   result_counter <- 1
#   
#   for (val in vary_vals) {
#     for (d in 1:n_draws) {
#       draw_row <- selected_draws[d, ]
#       intervened_input_data <- baseline_data_row
#       
#       # Update the variable being varied and its quadratic (if applicable)
#       if (vary_var == "Year") {
#         intervened_input_data[["Year"]] <- val
#         if ("Year_2" %in% names(intervened_input_data)) {
#           intervened_input_data[["Year_2"]] <- val^2
#         }
#       } else if (vary_var == "BreachDays") {
#         intervened_input_data[["BreachDays"]] <- val
#         if ("BreachDays_2" %in% names(intervened_input_data)) {
#           intervened_input_data[["BreachDays_2"]] <- val^2
#         }
#       } else if (vary_var == "Temp") {
#         intervened_input_data[["Temp"]] <- val
#         if ("Temp_2" %in% names(intervened_input_data)) {
#           intervened_input_data[["Temp_2"]] <- val^2
#         }
#       } else if (vary_var == "SAV") {
#         intervened_input_data[["SAV"]] <- val
#         if ("SAV_2" %in% names(intervened_input_data)) {
#           intervened_input_data[["SAV_2"]] <- val^2
#         }
#       } else {
#         # For other variables (DO, SC_count, SB_count, Rain), just set the value
#         intervened_input_data[[vary_var]] <- val
#       }
#       
#       # Calculate predicted log(count)
#       mu_goby_pred_log_count <- calculate_mu_goby_deterministic(
#         draw_row = draw_row,
#         input_data = intervened_input_data,
#         initial_intervention_var = vary_var, 
#         initial_intervention_value = val      
#       )
#       
#       # Convert log(count) to log(density) by subtracting log(Area)
#       log_area_offset <- baseline_data_row$Area
#       mu_goby_pred_log_density <- mu_goby_pred_log_count - log_area_offset
#       
#       # Store result
#       results_df$mu_goby_pred_log_density[result_counter] <- mu_goby_pred_log_density
#       result_counter <- result_counter + 1
#     }
#   }
#   
#   # Summarize across draws
#   summary_df <- results_df %>%
#     dplyr::group_by(var_value) %>%
#     dplyr::summarise(
#       mean_log_density = mean(mu_goby_pred_log_density, na.rm = TRUE),
#       lower_log_density = quantile(mu_goby_pred_log_density, 0.055, na.rm = TRUE),
#       upper_log_density = quantile(mu_goby_pred_log_density, 0.945, na.rm = TRUE),
#       .groups = 'drop'
#     ) %>%
#     dplyr::mutate(
#       mean_goby_density = exp(mean_log_density), 
#       lower_goby_density = exp(lower_log_density), 
#       upper_goby_density = exp(upper_log_density),
#       Zone = target_zone,
#       predictor = vary_var
#     )
#   
#   return(summary_df)
# }
# 
# # --- 4. Define Variables and Their Ranges for Plots ---
# sb_vals <- c(0, 1)
# sc_vals <- c(0, 1)
# breachdays_vals <- seq(min(df.data$BreachDays, na.rm = TRUE), 
#                        max(df.data$BreachDays, na.rm = TRUE), 
#                        length.out = 50)
# rain_vals <- seq(min(df.data$Rain, na.rm = TRUE), 
#                  max(df.data$Rain, na.rm = TRUE), 
#                  length.out = 50)
# do_vals <- seq(min(df.data$DO, na.rm = TRUE), 
#                max(df.data$DO, na.rm = TRUE), 
#                length.out = 50)
# sav_vals <- seq(min(df.data$SAV, na.rm = TRUE), 
#                 max(df.data$SAV, na.rm = TRUE), 
#                 length.out = 50)
# year_vals <- seq(min(df.data$Year, na.rm = TRUE), 
#                  max(df.data$Year, na.rm = TRUE), 
#                  length.out = 50)
# 
# # --- 5. Calculate Marginal Effects for Each Variable AND Each Zone ---
# all_zones <- sort(unique(df.data$Zone))
# all_marginal_effects <- list()
# 
# predictors_to_plot <- list(
#   "SB_count" = sb_vals,
#   "BreachDays" = breachdays_vals,
#   "Rain" = rain_vals,
#   "DO" = do_vals, 
#   "SAV" = sav_vals, 
#   "SC_count" = sc_vals,
#   "Year" = year_vals
# )
# 
# # Calculate with error handling
# for (current_zone in all_zones) {
#   for (predictor_name in names(predictors_to_plot)) {
#     
#     tryCatch({
#       vals <- predictors_to_plot[[predictor_name]]
#       
#       me_result <- calculate_marginal_effect_goby_by_zone(
#         vary_var = predictor_name,
#         vary_vals = vals,
#         draws_df = relevant_params,
#         data_df = df.data,
#         target_zone = current_zone,
#         n_draws = 40
#       )
#       
#       me_result$plot_id <- predictor_name
#       all_marginal_effects[[length(all_marginal_effects) + 1]] <- me_result
#       
#     }, error = function(e) {
#       message(paste("Error for", predictor_name, "in zone", current_zone, ":", e$message))
#     })
#   }
# }
# 
# full_me_df <- do.call(rbind, all_marginal_effects)
# 
# # --- 6. Data Processing and Plot Creation ---
# do_data_sb <- full_me_df %>% 
#   filter(plot_id == "DO") %>% 
#   mutate(plot_id = "DO_SB_Goby")
# 
# do_data_sc <- full_me_df %>% 
#   filter(plot_id == "DO") %>% 
#   mutate(plot_id = "DO_SC_Goby")
# 
# sav_data_sc <- full_me_df %>% 
#   filter(plot_id == "SAV") %>% 
#   mutate(plot_id = "SAV_SC_Goby")
# 
# full_me_df_processed <- full_me_df %>%
#   filter(!plot_id %in% c("DO", "SAV")) %>% 
#   mutate(
#     plot_id = case_when(
#       plot_id == "BreachDays" ~ "BreachDays_Goby",
#       plot_id == "Rain" ~ "Rain_Goby",
#       plot_id == "Year" ~ "Year2_Goby",
#       TRUE ~ plot_id 
#     )
#   ) %>%
#   bind_rows(do_data_sb, do_data_sc, sav_data_sc) 
# 
# zone_labels <- c("1" = "East", "2" = "Northwest", "3" = "West")
# 
# plot_titles <- c(
#   "SB_count" = "SB → Goby",
#   "BreachDays_Goby" = "BreachDays (Total Effect) → Goby",
#   "Rain_Goby" = "Rain (Total Effect) → Goby",
#   "DO_SB_Goby" = "DO → SB → Goby",
#   "DO_SC_Goby" = "DO → SC → Goby",
#   "SAV_SC_Goby" = "SAV → SC → Goby",
#   "SC_count" = "SC → Goby",
#   "Year2_Goby" = "Year & Year² → Goby"
# )
# 
# x_axis_labels <- c(
#   "SB_count" = "SB (0 or 1)",
#   "BreachDays" = "BreachDays",
#   "Rain" = "Rain",
#   "DO" = "Dissolved Oxygen (DO)",
#   "SAV" = "Submerged Aquatic Vegetation (SAV)",
#   "SC_count" = "SC (0 or 1)",
#   "Year" = "Year"
# )
# 
# raw_data_alpha <- 0.3
# raw_data_color <- "grey30"
# raw_data_size <- 0.8
# jitter_width <- 0.2
# 
# make_marginal_plot <- function(plot_identifier, full_df, raw_df) {
#   
#   plot_data <- full_df %>% dplyr::filter(plot_id == plot_identifier)
#   
#   if (nrow(plot_data) == 0) {
#     message(paste("Warning: No data for plot", plot_identifier))
#     return(NULL)
#   }
#   
#   original_predictor_name <- as.character(plot_data$predictor[1])
#   is_binary <- original_predictor_name %in% c("SB_count", "SC_count")
#   
#   p <- ggplot(plot_data, aes(x = var_value, y = mean_goby_density)) +
#     geom_point(
#       data = raw_df, 
#       aes_string(x = original_predictor_name, y = "Goby / exp(Area)"),
#       alpha = raw_data_alpha, 
#       color = raw_data_color, 
#       size = raw_data_size,
#       position = if (is_binary) position_jitter(width = jitter_width, height = 0) else "identity"
#     ) +
#     geom_line(color = "#0072B2", linewidth = 1) +
#     geom_ribbon(
#       aes(ymin = lower_goby_density, ymax = upper_goby_density), 
#       fill = "#0072B2", 
#       alpha = 0.2
#     ) +
#     labs(
#       title = plot_titles[plot_identifier],
#       x = x_axis_labels[original_predictor_name],
#       y = NULL
#     ) +
#     theme_minimal(base_size = 16) +
#     facet_wrap(~ Zone, ncol = 3, labeller = labeller(Zone = zone_labels)) +
#     coord_cartesian(ylim = c(0, 25))
#   
#   if (is_binary) {
#     p <- p + scale_x_continuous(breaks = c(0, 1))
#   }
#   
#   return(p)
# }
# 
# # Create plots
# plot_sb_goby <- make_marginal_plot("SB_count", full_me_df_processed, df.data)
# plot_breachdays_goby <- make_marginal_plot("BreachDays_Goby", full_me_df_processed, df.data)
# plot_rain_goby <- make_marginal_plot("Rain_Goby", full_me_df_processed, df.data)
# plot_do_sb_goby <- make_marginal_plot("DO_SB_Goby", full_me_df_processed, df.data)
# plot_do_sc_goby <- make_marginal_plot("DO_SC_Goby", full_me_df_processed, df.data)
# plot_sav_sc_goby <- make_marginal_plot("SAV_SC_Goby", full_me_df_processed, df.data)
# plot_sc_goby <- make_marginal_plot("SC_count", full_me_df_processed, df.data)
# plot_year2_goby <- make_marginal_plot("Year2_Goby", full_me_df_processed, df.data)
# 
# # --- 7. Arrange plots ---
# final_grid <- plot_grid(
#   plot_sb_goby,
#   plot_breachdays_goby,
#   plot_rain_goby,
#   plot_do_sb_goby,
#   plot_do_sc_goby,
#   plot_sav_sc_goby,
#   plot_sc_goby,
#   plot_year2_goby,
#   ncol = 2,
#   labels = "AUTO"
# )
# 
# # Create y-axis label
# y_axis_label <- ggdraw() +
#   draw_label(
#     "Expected Goby Density (m²)",
#     fontface = 'bold',
#     size = 16,
#     angle = 90,
#     x = 0.5,
#     y = 0.5
#   )
# 
# # Combine
# final_plot_with_ylabel <- plot_grid(
#   y_axis_label,
#   final_grid,
#   ncol = 2,
#   scale = 0.9,
#   rel_widths = c(0.04, 1)
# )
# 
# print(final_plot_with_ylabel)
# 
# # Save
# ggsave(
#   filename = "Output/Plots/MarginalEffects_Grid_Corrected_Labels.png",
#   plot = final_plot_with_ylabel,
#   width = 16,
#   height = 18,
#   dpi = 300
# )
# 
# 
