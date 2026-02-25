#updated 2025-12-28

# OPTIMIZED CAUSAL INFERENCE CODE - PARALLELIZED VERSION
# Much faster using furrr (parallel purrr) and optimized data structures

# --- 0. Setup ---
library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(forcats)
library(stringr)
library(MASS)
library(furrr)  # Parallel version of purrr
library(future) # For parallel processing

# Set up parallel processing
# Use all available cores minus 1 (to keep system responsive)
plan(multisession, workers = parallel::detectCores() - 1)

df.data$Goby_Density <- df.data$Goby / exp(df.data$Area)

# --- 1. Extract Posterior Draws ---
draws <- mod.SEM.real$draws(format = "df")
relevant_params <- draws

# Calculate standard deviations
numeric_cols <- names(df.data)[sapply(df.data, is.numeric)]
sds <- data.frame(
  purrr::map_dfc(df.data[numeric_cols], ~ sd(.x, na.rm = TRUE))
)
colnames(sds) <- numeric_cols

# --- 2. Define Causal Paths ---
causal_paths <- list(
  # Direct Effects
  "Micro -> Goby" = c("Micro", "Goby"),
  "Substrate -> Goby" = c("Substrate", "Goby"),
  "Year -> Goby" = c("Year", "Goby"),
  "Year_2 -> Goby" = c("Year_2", "Goby"),
  "Temp_2 -> Goby" = c("Temp_2", "Goby"),
  "BreachDays_2 -> Goby" = c("BreachDays_2", "Goby"),
  "Goby_lag -> Goby" = c("Goby_lag", "Goby"),
  "SAV -> Goby" = c("SAV", "Goby"),
  "SAV_2 -> Goby" = c("SAV_2", "Goby"),
  "DO -> Goby" = c("DO", "Goby"),
  "BreachDays -> Goby" = c("BreachDays", "Goby"),
  "Temp -> Goby" = c("Temp", "Goby"),
  "SC_count -> Goby" = c("SC_count", "Goby"),
  "SB_count -> Goby" = c("SB_count", "Goby"),
  
  # Mediated Effects
  "Rain -> BreachDays -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Goby"),
  "Rain -> BreachDays -> Temp -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2", "Goby"),
  "Rain -> BreachDays -> Temp -> DO -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2", "DO", "Goby"),
  "Rain -> BreachDays -> Temp -> DO -> SAV -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2", "DO", "SAV", "Goby"),
  "Wind -> Temp -> Goby" = c("Wind", "Temp", "Temp_2", "Goby"),
  "Wind -> DO -> Goby" = c("Wind", "DO", "Goby"),
  "Wind -> Temp -> DO -> Goby" = c("Wind", "Temp", "Temp_2", "DO", "Goby"),
  "Wind -> DO -> SAV -> Goby" = c("Wind", "DO", "SAV", "Goby"),
  "Wind -> Temp -> SAV -> Goby" = c("Wind", "Temp", "Temp_2", "SAV", "Goby"),
  "BreachDays -> Temp -> Goby" = c("BreachDays", "BreachDays_2", "Temp", "Temp_2", "Goby"),
  "BreachDays -> Temp -> DO -> Goby" = c("BreachDays", "BreachDays_2", "Temp", "Temp_2", "DO", "Goby"),
  "Temp -> DO -> Goby" = c("Temp", "DO", "Goby"),
  "Temp -> SAV -> Goby" = c("Temp", "SAV", "Goby"),
  "Substrate -> SC_count -> Goby" = c("Substrate", "SC_count", "Goby"),
  "DO -> SAV -> Goby" = c("DO", "SAV", "Goby"),
  "DO -> SC_count -> Goby" = c("DO", "SC_count", "Goby"),
  "DO -> SB_count -> Goby" = c("DO", "SB_count", "Goby"),
  "SAV -> SC_count -> Goby" = c("SAV", "SC_count", "Goby"),
  "SAV -> SB_count -> Goby" = c("SAV", "SB_count", "Goby")
)

# --- 3. OPTIMIZED: Vectorized calculation function ---
# Key optimization: Process all draws at once using matrix operations where possible
calculate_mu_goby_vectorized <- function(draws_matrix, input_data, path, 
                                         initial_intervention_var = NULL, 
                                         initial_intervention_value = NULL) {
  
  n_draws <- nrow(draws_matrix)
  
  # Pre-allocate result vector
  mu_goby_results <- numeric(n_draws)
  
  # Convert input_data to named list once
  current_values <- as.list(input_data)
  if (!is.null(initial_intervention_var)) {
    current_values[[initial_intervention_var]] <- initial_intervention_value
  }
  
  # Helper function
  get_val <- function(var_name) {
    if (!is.null(initial_intervention_var) && var_name == initial_intervention_var) {
      return(current_values[[var_name]])
    } else if (var_name %in% path) {
      return(current_values[[var_name]])
    } else {
      return(input_data[[var_name]])
    }
  }
  
  is_intervened <- function(var) !is.null(initial_intervention_var) && initial_intervention_var == var
  
  # Pre-extract zone index
  zone_idx <- input_data$Zone
  U_bio_col <- paste0("U_bio[", zone_idx, "]")
  U_phys_col <- paste0("U_phys[", zone_idx, "]")
  a_Goby_col <- paste0("a_Goby[", zone_idx, "]")
  
  # VECTORIZED: Extract all parameters at once
  k_U_phys_vec <- draws_matrix[, "k_U_phys"]
  k_U_bio_vec <- draws_matrix[, "k_U_bio"]
  U_bio_vec <- draws_matrix[, U_bio_col]
  U_phys_vec <- draws_matrix[, U_phys_col]
  
  # Calculate mediators (still need some conditional logic, but vectorized within each condition)
  
  # BreachDays (physical)
  if ("BreachDays" %in% path || is_intervened("BreachDays") || is_intervened("Rain")) {
    Breach_nu <- draws_matrix[, "a_BreachDays"] + 
      draws_matrix[, "beta_Rain"] * get_val("Rain") + 
      k_U_phys_vec * U_phys_vec
    current_values$BreachDays <- mean(Breach_nu)  # Store scalar for get_val
    current_values$BreachDays_2 <- mean(Breach_nu^2)
    BreachDays_vec <- Breach_nu  # Store vector for final calculation
    BreachDays_2_vec <- Breach_nu^2
  } else {
    BreachDays_vec <- rep(get_val("BreachDays"), n_draws)
    BreachDays_2_vec <- rep(get_val("BreachDays_2"), n_draws)
  }
  
  # Temp (physical)
  if ("Temp" %in% path || is_intervened("Temp") || is_intervened("BreachDays") || is_intervened("Wind")) {
    Temp_nu <- draws_matrix[, "a_Temp"] + 
      draws_matrix[, "beta_BreachDays_vec[2]"] * BreachDays_vec + 
      draws_matrix[, "beta_Wind_vec[2]"] * get_val("Wind") + 
      k_U_phys_vec * U_phys_vec
    current_values$Temp <- mean(Temp_nu)
    current_values$Temp_2 <- mean(Temp_nu^2)
    Temp_vec <- Temp_nu
    Temp_2_vec <- Temp_nu^2
  } else {
    Temp_vec <- rep(get_val("Temp"), n_draws)
    Temp_2_vec <- rep(get_val("Temp_2"), n_draws)
  }
  
  # DO (physical)
  if ("DO" %in% path || is_intervened("DO") || is_intervened("Temp") || is_intervened("Wind")) {
    DO_nu <- draws_matrix[, "a_DO"] + 
      draws_matrix[, "beta_Temp_vec[2]"] * Temp_vec + 
      draws_matrix[, "beta_Wind_vec[1]"] * get_val("Wind") + 
      k_U_phys_vec * U_phys_vec
    current_values$DO <- mean(DO_nu)
    DO_vec <- DO_nu
  } else {
    DO_vec <- rep(get_val("DO"), n_draws)
  }
  
  # SAV (biological)
  if ("SAV" %in% path || is_intervened("SAV") || is_intervened("DO") || is_intervened("Temp")) {
    SAV_nu <- draws_matrix[, "a_SAV"] + 
      draws_matrix[, "beta_DO_vec[4]"] * DO_vec + 
      draws_matrix[, "beta_Temp_vec[3]"] * Temp_vec + 
      k_U_bio_vec * U_bio_vec
    current_values$SAV <- mean(SAV_nu)
    current_values$SAV_2 <- mean(SAV_nu^2)
    SAV_vec <- SAV_nu
    SAV_2_vec <- SAV_nu^2
  } else {
    SAV_vec <- rep(get_val("SAV"), n_draws)
    SAV_2_vec <- rep(get_val("SAV_2"), n_draws)
  }
  
  # SC_count (biological)
  if ("SC_count" %in% path || is_intervened("SC_count") || is_intervened("Substrate") || is_intervened("DO") || is_intervened("SAV")) {
    SC_mu_logit <- draws_matrix[, "a_SC"] + 
      draws_matrix[, "beta_Substrate_vec[2]"] * get_val("Substrate") + 
      draws_matrix[, "beta_DO_vec[3]"] * DO_vec + 
      draws_matrix[, "beta_SAV_vec[3]"] * SAV_vec + 
      k_U_bio_vec * U_bio_vec
    SC_vec <- plogis(SC_mu_logit)
    current_values$SC_count <- mean(SC_vec)
  } else {
    SC_vec <- rep(get_val("SC_count"), n_draws)
  }
  
  # SB_count (biological)
  if ("SB_count" %in% path || is_intervened("SB_count") || is_intervened("DO") || is_intervened("SAV")) {
    SB_mu_logit <- draws_matrix[, "a_SB"] + 
      draws_matrix[, "beta_DO_vec[2]"] * DO_vec + 
      draws_matrix[, "beta_SAV_vec[2]"] * SAV_vec + 
      k_U_bio_vec * U_bio_vec
    SB_vec <- plogis(SB_mu_logit)
    current_values$SB_count <- mean(SB_vec)
  } else {
    SB_vec <- rep(get_val("SB_count"), n_draws)
  }
  
  # VECTORIZED: Final Goby calculation for all draws at once
  mu_goby_results <- draws_matrix[, a_Goby_col] +
    draws_matrix[, "beta_SAV_vec[1]"] * SAV_vec +
    draws_matrix[, "beta_SAV_2"] * SAV_2_vec + 
    draws_matrix[, "beta_DO_vec[1]"] * DO_vec +
    draws_matrix[, "beta_BreachDays_vec[1]"] * BreachDays_vec +
    draws_matrix[, "beta_Temp_vec[1]"] * Temp_vec +
    draws_matrix[, "beta_SC_count"] * SC_vec +
    draws_matrix[, "beta_SB_count"] * SB_vec +
    draws_matrix[, "beta_Substrate_vec[1]"] * get_val("Substrate") +
    draws_matrix[, "beta_Micro"] * get_val("Micro") +
    draws_matrix[, "beta_Year"] * get_val("Year") +
    draws_matrix[, "beta_Year_2"] * get_val("Year_2") +
    draws_matrix[, "beta_Temp_2"] * Temp_2_vec +
    draws_matrix[, "beta_BreachDays_2"] * BreachDays_2_vec +
    draws_matrix[, "beta_Goby_lag"] * get_val("Goby_lag") +
    get_val("Area")
  
  return(mu_goby_results)
}

# --- 4. OPTIMIZED: Calculate path effect with vectorization ---
calculate_path_effect_optimized <- function(path_name, path, draws_df, data_df, sds_df, n_draws = 100) {
  message(paste("Calculating effect for path:", path_name))
  start_node <- path[1]
  
  # 1. Prepare Baseline Data
  baseline_data_row <- data_df %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~median(.x, na.rm = TRUE))) %>%
    dplyr::mutate(Zone = as.integer(round(Zone))) %>%
    as.data.frame()
  
  # 2. Prepare Intervened Data
  intervened_data_row <- baseline_data_row
  
  if (start_node %in% names(sds_df)) {
    intervened_data_row[[start_node]] <- baseline_data_row[[start_node]] + sds_df[[1, start_node]]
  } else if (start_node %in% c("SC_count", "SB_count")) {
    intervened_data_row[[start_node]] <- 1
  } else {
    intervened_data_row[[start_node]] <- baseline_data_row[[start_node]] + sds_df[[1, start_node]]
  }
  
  # 3. Select Draws and Convert to Matrix (FASTER than repeated row access)
  total_draws_available <- nrow(draws_df)
  if (n_draws > total_draws_available) n_draws <- total_draws_available
  sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
  selected_draws <- draws_df[sample_indices, ]
  
  # Convert to matrix for faster access
  draws_matrix <- as.matrix(selected_draws)
  
  # 4. VECTORIZED: Calculate both baseline and intervention for all draws at once
  mu_goby_baseline <- calculate_mu_goby_vectorized(
    draws_matrix, 
    baseline_data_row, 
    path = c()
  )
  
  mu_goby_intervene <- calculate_mu_goby_vectorized(
    draws_matrix, 
    baseline_data_row, 
    path = path, 
    initial_intervention_var = start_node, 
    initial_intervention_value = intervened_data_row[[start_node]]
  )
  
  # Return the effect distribution
  return(mu_goby_intervene - mu_goby_baseline)
}

# --- 5. PARALLELIZED: Loop over paths ---
cat("\n=== Computing Posterior Path Effects (Parallelized) ===\n")

# Use future_map instead of regular map for parallel processing
all_path_effects_posterior <- future_map(
  names(causal_paths),
  function(path_name) {
    current_path <- causal_paths[[path_name]]
    start_node_for_sd <- current_path[1]
    
    # Check if we have SD
    if (!(start_node_for_sd %in% c("SB_count", "SC_count")) && 
        !(start_node_for_sd %in% names(sds))) {
      return(NULL)
    }
    
    calculate_path_effect_optimized(
      path_name, 
      current_path, 
      relevant_params, 
      df.data, 
      sds, 
      n_draws = 50  # Increased since it's faster now
    )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE  # Show progress bar
)

# Name the list
names(all_path_effects_posterior) <- names(causal_paths)

# Remove NULL entries
all_path_effects_posterior <- all_path_effects_posterior[!sapply(all_path_effects_posterior, is.null)]

# Summarize
effects_summary <- map_df(all_path_effects_posterior, ~{
  tibble(
    mean = mean(.x, na.rm = TRUE), 
    lower = quantile(.x, 0.055, na.rm = TRUE), 
    upper = quantile(.x, 0.945, na.rm = TRUE)
  )
}, .id = "path") %>%
  mutate(type = "Posterior")

# ==============================================================================
# 6. GENERATE PRIORS (Can also be parallelized)
# ==============================================================================

cat("\n=== Generating Prior Predictive Distributions ===\n")

n_prior_draws <- 1000

# Pre-allocate tibble with all draws at once (FASTER)
prior_draws_df <- tibble(.chain = 1, .iteration = 1:n_prior_draws, .draw = 1:n_prior_draws)

# VECTORIZED: Generate all priors at once
prior_draws_df$beta_Rain <- rnorm(n_prior_draws, 0.25, 0.25)
prior_draws_df$beta_Micro <- rnorm(n_prior_draws, 0.25, 0.25)
prior_draws_df$beta_Year <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$beta_Year_2 <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$beta_Temp_2 <- rnorm(n_prior_draws, -0.1, 0.25)
prior_draws_df$beta_BreachDays_2 <- rnorm(n_prior_draws, -0.10, 0.25)
prior_draws_df$beta_SAV_2 <- rnorm(n_prior_draws, -0.1, 0.25)
prior_draws_df$beta_Goby_lag <- rnorm(n_prior_draws, 0.25, 0.25)
prior_draws_df$beta_SC_count <- rnorm(n_prior_draws, -0.1, 0.25)
prior_draws_df$beta_SB_count <- rnorm(n_prior_draws, 0, 0.5)

# Hierarchical coefficients
prior_draws_df$`beta_DO_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_DO_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_DO_vec[3]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_DO_vec[4]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$`beta_SAV_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_SAV_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_SAV_vec[3]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$`beta_Temp_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_Temp_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_Temp_vec[3]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$`beta_BreachDays_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_BreachDays_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$`beta_Wind_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_Wind_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$`beta_Substrate_vec[1]` <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$`beta_Substrate_vec[2]` <- rnorm(n_prior_draws, 0, 0.5)

prior_draws_df$a_BreachDays <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$a_Temp <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$a_DO <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$a_SAV <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$a_SC <- rnorm(n_prior_draws, 0, 0.5)
prior_draws_df$a_SB <- rnorm(n_prior_draws, 0, 0.5)

# Hierarchical Goby intercepts
a_bar_Goby_prior <- rnorm(n_prior_draws, 0, 0.5)
sigma_Goby_prior <- rexp(n_prior_draws, 1)

for(z in 1:3) {
  prior_draws_df[[paste0("a_Goby[", z, "]")]] <- rnorm(n_prior_draws, a_bar_Goby_prior, sigma_Goby_prior)
  prior_draws_df[[paste0("U_bio[", z, "]")]] <- rnorm(n_prior_draws, 0, 1)
  prior_draws_df[[paste0("U_phys[", z, "]")]] <- rnorm(n_prior_draws, 0, 1)
}

prior_draws_df$k_U_phys <- rexp(n_prior_draws, rate = 1)
prior_draws_df$k_U_bio <- rexp(n_prior_draws, rate = 1)

# PARALLELIZED: Prior path effects
all_path_effects_prior <- future_map(
  names(causal_paths),
  function(path_name) {
    current_path <- causal_paths[[path_name]]
    start_node_for_sd <- current_path[1]
    
    if (!(start_node_for_sd %in% c("SB_count", "SC_count")) && 
        !(start_node_for_sd %in% names(sds))) {
      return(NULL)
    }
    
    calculate_path_effect_optimized(
      path_name, 
      current_path, 
      draws_df = prior_draws_df,
      data_df = df.data, 
      sds_df = sds, 
      n_draws = n_prior_draws
    )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

names(all_path_effects_prior) <- names(causal_paths)
all_path_effects_prior <- all_path_effects_prior[!sapply(all_path_effects_prior, is.null)]

priors_summary <- map_df(all_path_effects_prior, ~{
  tibble(
    mean = mean(.x, na.rm = TRUE), 
    lower = quantile(.x, 0.055, na.rm = TRUE), 
    upper = quantile(.x, 0.945, na.rm = TRUE)
  )
}, .id = "path") %>%
  mutate(type = "Prior")

# --- 7. Combine and Plot (Same as before) ---
combined_summary <- bind_rows(effects_summary, priors_summary) %>%
  mutate(
    path_label = str_replace_all(path, " -> ", " → "),
    path_label = str_replace_all(path_label, "SB_count", "SB"),
    path_label = str_replace_all(path_label, "SC_count", "SC"),
    path_label = str_replace_all(path_label, "_2", "²")
  ) %>%
  mutate(
    significant_color = case_when(
      type == "Posterior" & lower > 0 & upper > 0 ~ "Positive",
      type == "Posterior" & lower < 0 & upper < 0 ~ "Negative",
      type == "Posterior" ~ "Overlap Zero",
      TRUE ~ NA_character_
    )
  )

path_order <- combined_summary %>%
  filter(type == "Posterior") %>%
  arrange(mean) %>%
  pull(path_label) %>%
  unique()

combined_summary$path_label <- factor(combined_summary$path_label, levels = path_order)
combined_summary_for_plot <- combined_summary %>% filter(!is.na(path_label))

# Create plot
causal_path_forest_plot <- ggplot(
  combined_summary_for_plot,
  aes(x = path_label, y = mean, ymin = lower, ymax = upper)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  geom_pointrange(
    aes(shape = type, color = significant_color, size = type, linewidth = type),
    position = position_nudge(x = ifelse(combined_summary_for_plot$type == "Posterior", 0.1, -0.1))
  ) +
  coord_flip() +
  scale_shape_manual(name = NULL, values = c("Prior" = 1, "Posterior" = 16)) +
  scale_color_manual(
    name = NULL,
    values = c("Positive" = "#0072B2", "Negative" = "#D55E00", "Overlap Zero" = "black"),
    na.value = "grey70"
  ) +
  scale_size_manual(name = NULL, values = c("Prior" = 0.6, "Posterior" = 1.0), guide = "none") +
  scale_linewidth_manual(name = NULL, values = c("Prior" = 0.6, "Posterior" = 1.0), guide = "none") +
  labs(
    x = "Causal Pathway",
    y = "log Δ in Goby Density per SD (β × ΔX)",
    title = "Causal Path Effects on Goby Density",
    subtitle = "Posterior estimates with 89% credible intervals (parallelized computation)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.22, 0.45),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.5)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(color = c("grey70", "black")))
  )

print(causal_path_forest_plot)

ggsave("Output/Plots/causal_path_forest_plot_optimized.png", 
       plot = causal_path_forest_plot, 
       width = 22, height = 30, units = "cm", dpi = 300, bg = "white")

# Clean up parallel workers
plan(sequential)

cat("\n=== Optimized Analysis Complete ===\n")
cat("Results saved to: Output/Plots/causal_path_forest_plot_optimized.png\n")


# #Old below
# 
# # Add this line before your ggplot code
# df.data$Goby_Density <- df.data$Goby / exp(df.data$Area)
# 
# # --- 0. Setup ---
# library(cmdstanr)
# library(posterior)
# library(dplyr)
# library(ggplot2)
# library(purrr)
# library(tibble)
# library(forcats)
# library(stringr) # Added for string manipulation
# 
# # Assuming 'mod.SEM.real' is your cmdstanr fitted object
# # and 'df.data' is the original dataframe you used to fit the model
# 
# # --- 1. Extract Posterior Draws ---
# draws <- mod.SEM.real$draws(format = "df")
# relevant_params <- draws
# 
# # Calculate standard deviations for all continuous variables
# numeric_cols <- names(df.data)[sapply(df.data, is.numeric)]
# sds <- data.frame(
#   purrr::map_dfc(df.data[numeric_cols], ~ sd(.x, na.rm = TRUE))
# )
# colnames(sds) <- numeric_cols
# 
# # --- 2. Define Causal Paths ---
# # !!! IMPORTANT: Assuming SAV_2 is now a direct predictor in your model !!!
# # You might need to add SAV_2 calculation to your data simulation/prep if not already done.
# # Assuming df.data contains SAV_2 = SAV^2
# causal_paths <- list(
#   # Direct Effects
#   "Micro -> Goby" = c("Micro", "Goby"),
#   "Substrate -> Goby" = c("Substrate", "Goby"),
#   "Year -> Goby" = c("Year", "Goby"),
#   "Year_2 -> Goby" = c("Year_2", "Goby"),
#   "Temp_2 -> Goby" = c("Temp_2", "Goby"),
#   "BreachDays_2 -> Goby" = c("BreachDays_2", "Goby"),
#   "Goby_lag -> Goby" = c("Goby_lag", "Goby"),
#   "SAV -> Goby" = c("SAV", "Goby"),
#   "SAV_2 -> Goby" = c("SAV_2", "Goby"), # <-- ADDED SAV_2 path
#   "DO -> Goby" = c("DO", "Goby"),
#   "BreachDays -> Goby" = c("BreachDays", "Goby"),
#   "Temp -> Goby" = c("Temp", "Goby"),
#   "SC_count -> Goby" = c("SC_count", "Goby"),
#   "SB_count -> Goby" = c("SB_count", "Goby"),
#   # Mediated Effects
#   "Rain -> BreachDays -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Goby"),
#   "Rain -> BreachDays -> Temp -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2",  "Goby"),
#   "Rain -> BreachDays -> Temp -> DO -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2",  "DO", "Goby"),
#   "Rain -> BreachDays -> Temp -> DO -> SAV -> Goby" = c("Rain", "BreachDays", "BreachDays_2", "Temp", "Temp_2", "DO", "SAV", "Goby"),
#   "Wind -> Temp -> Goby" = c("Wind", "Temp", "Temp_2", "Goby"),
#   "Wind -> DO -> Goby" = c("Wind", "DO", "Goby"),
#   "Wind -> Temp -> DO -> Goby" = c("Wind", "Temp", "Temp_2", "DO", "Goby"),
#   "Wind -> DO -> SAV -> Goby" = c("Wind", "DO", "SAV", "Goby"),
#   "Wind -> Temp -> SAV -> Goby" = c("Wind", "Temp", "Year_2", "SAV", "Goby"),
#   "BreachDays -> Temp -> Goby" = c("BreachDays", "BreachDays_2", "Temp", "Temp_2",  "Goby"),
#   "BreachDays -> Temp -> DO -> Goby" = c("BreachDays", "BreachDays_2", "Temp", "Temp_2",  "DO", "Goby"),
#   "Temp -> DO -> Goby" = c("Temp", "DO", "Goby"),
#   "Temp -> SAV -> Goby" = c("Temp", "SAV", "Goby"),
#   "Substrate -> SC_count -> Goby" = c("Substrate", "SC_count", "Goby"),
#   "DO -> SAV -> Goby" = c("DO", "SAV", "Goby"),
#   "DO -> SC_count -> Goby" = c("DO", "SC_count", "Goby"),
#   "DO -> SB_count -> Goby" = c("DO", "SB_count", "Goby"),
#   "SAV -> SC_count -> Goby" = c("SAV", "SC_count", "Goby"),
#   "SAV -> SB_count -> Goby" = c("SAV", "SB_count", "Goby")
# )
# 
# # --- 3. Deterministic Path-Specific Linear Predictor Function ---
# calculate_mu_goby_deterministic <- function(draw_row, input_data, path, initial_intervention_var = NULL, initial_intervention_value = NULL) {
#   current_values <- as.list(input_data)
#   if (!is.null(initial_intervention_var)) {
#     current_values[[initial_intervention_var]] <- initial_intervention_value
#   }
#   
#   get_val <- function(var_name) {
#     if (!is.null(initial_intervention_var) && var_name == initial_intervention_var) return(current_values[[var_name]])
#     else if (var_name %in% path) return(current_values[[var_name]])
#     # For baseline, pull from input_data if not intervened/on path
#     else return(input_data[[var_name]])
#   }
#   
#   is_intervened <- function(var) !is.null(initial_intervention_var) && initial_intervention_var == var
#   
#   # --- Deterministically Calculate Mediators in DAG order ---
#   if ("BreachDays" %in% path || is_intervened("BreachDays") || is_intervened("Rain")) { # Recalc if Rain changed
#     Breach_nu <- draw_row[["a_BreachDays"]] + draw_row[["beta_Rain"]] * get_val("Rain") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$BreachDays <- Breach_nu
#   }
#   if ("Temp" %in% path || is_intervened("Temp") || is_intervened("BreachDays") || is_intervened("Wind")) { # Recalc if parents changed
#     Temp_nu <- draw_row[["a_Temp"]] + draw_row[["beta_BreachDays_vec[2]"]] * get_val("BreachDays") + draw_row[["beta_Wind_vec[2]"]] * get_val("Wind") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$Temp <- Temp_nu
#   }
#   if ("DO" %in% path || is_intervened("DO") || is_intervened("Temp") || is_intervened("Wind")) { # Recalc if parents changed
#     DO_nu <- draw_row[["a_DO"]] + draw_row[["beta_Temp_vec[2]"]] * get_val("Temp") + draw_row[["beta_Wind_vec[1]"]] * get_val("Wind") + draw_row[["k_U"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$DO <- DO_nu
#   }
#   if ("SAV" %in% path || is_intervened("SAV") || is_intervened("DO") || is_intervened("Temp")) { # Recalc if parents changed
#     SAV_nu <- draw_row[["a_SAV"]] + draw_row[["beta_DO_vec[4]"]] * get_val("DO") + draw_row[["beta_Temp_vec[3]"]] * get_val("Temp") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SAV <- SAV_nu
#     # Update SAV_2 if SAV changes
#     if ("SAV_2" %in% names(current_values)) {
#       current_values$SAV_2 <- current_values$SAV^2
#     }
#   }
#   if ("SC_count" %in% path || is_intervened("SC_count") || is_intervened("Substrate") || is_intervened("DO") || is_intervened("SAV")) { # Recalc if parents changed
#     SC_mu_logit <- draw_row[["a_SC"]] + draw_row[["beta_Substrate_vec[2]"]] * get_val("Substrate") + draw_row[["beta_DO_vec[3]"]] * get_val("DO") + draw_row[["beta_SAV_vec[3]"]] * get_val("SAV") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SC_count <- plogis(SC_mu_logit)
#   }
#   if ("SB_count" %in% path || is_intervened("SB_count") || is_intervened("DO") || is_intervened("SAV")) { # Recalc if parents changed
#     SB_mu_logit <- draw_row[["a_SB"]] + draw_row[["beta_DO_vec[2]"]] * get_val("DO") + draw_row[["beta_SAV_vec[2]"]] * get_val("SAV") + draw_row[["k_U_bio"]] * draw_row[[paste0("U[", input_data$Zone, "]")]]
#     current_values$SB_count <- plogis(SB_mu_logit)
#   }
#   
#   mu_goby_val <- draw_row[[paste0("a_Goby[", current_values$Zone, "]")]] +
#     draw_row[["beta_SAV_vec[1]"]] * get_val("SAV") +
#     draw_row[["beta_SAV_2"]] * get_val("SAV_2") + 
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
#   return(mu_goby_val)
# }
# 
# # --- 4. Main Function to Calculate Path Effects (Final Robust Version) ---
# calculate_path_effect <- function(path_name, path, draws_df, data_df, sds_df, n_draws = 100) {
#   message(paste("Calculating effect for path:", path_name))
#   start_node <- path[1]
#   
#   # 1. Prepare Baseline Data (Median)
#   baseline_data_row <- data_df %>%
#     dplyr::summarise(dplyr::across(dplyr::everything(), ~median(.x, na.rm = TRUE))) %>%
#     dplyr::mutate(Zone = as.integer(round(Zone))) %>%
#     as.data.frame()
#   
#   # 2. Prepare Intervened Data (Median + SD)
#   intervened_data_row <- baseline_data_row
#   
#   # Calculate delta_X based on variable type
#   if (start_node %in% names(sds_df)) {
#     # Continuous: Add 1 SD
#     intervened_data_row[[start_node]] <- baseline_data_row[[start_node]] + sds_df[[1, start_node]]
#     
#     # Note: We do NOT need to manually update squared terms here (e.g. Year_2).
#     # The 'calculate_mu_goby_deterministic' function now handles that update 
#     # automatically when 'initial_intervention_var' is passed.
#     
#   } else if (start_node %in% c("SC_count", "SB_count")) {
#     # Binary: Toggle to 1
#     intervened_data_row[[start_node]] <- 1
#   } else {
#     # Fallback for squared terms if intervened directly (rare)
#     intervened_data_row[[start_node]] <- baseline_data_row[[start_node]] + sds_df[[1, start_node]]
#   }
#   
#   # 3. Select Draws
#   total_draws_available <- nrow(draws_df)
#   if (n_draws > total_draws_available) n_draws <- total_draws_available
#   sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
#   selected_draws <- draws_df[sample_indices, ]
#   
#   effect_dist <- numeric(n_draws)
#   
#   # 4. Simulation Loop (Used for ALL paths now)
#   # We removed the 'if (length(path) == 2)' block because it failed to account 
#   # for quadratic terms in direct effects.
#   
#   for (d in 1:n_draws) {
#     draw_row <- selected_draws[d, ]
#     
#     # Baseline Prediction
#     mu_goby_baseline <- calculate_mu_goby_deterministic(
#       draw_row, 
#       baseline_data_row, 
#       path = c() # Empty path means just predict using baseline values
#     )
#     
#     # Intervened Prediction
#     # This correctly updates linear AND quadratic terms for the start_node
#     mu_goby_intervene <- calculate_mu_goby_deterministic(
#       draw_row, 
#       baseline_data_row, 
#       path = path, 
#       initial_intervention_var = start_node, 
#       initial_intervention_value = intervened_data_row[[start_node]]
#     )
#     
#     effect_dist[d] <- mu_goby_intervene - mu_goby_baseline
#   }
#   
#   return(effect_dist)
# }
# 
# # --- 5. Loop and Summarize POSTERIORS ---
# all_path_effects_posterior <- list()
# 
# for (path_name in names(causal_paths)) {
#   current_path <- causal_paths[[path_name]]
#   # Error handling for missing SDs (e.g., if SAV_2 wasn't in numeric_cols)
#   start_node_for_sd <- current_path[1]
#   if (!(start_node_for_sd %in% c("SB_count", "SC_count")) && !(start_node_for_sd %in% names(sds))) {
#     message(paste("Skipping path:", path_name, "- SD not found for", start_node_for_sd))
#     next # Skip this path
#   }
#   all_path_effects_posterior[[path_name]] <- calculate_path_effect(path_name, current_path, relevant_params, df.data, sds, n_draws = 50) # Reduced n_draws
# }
# 
# effects_summary <- purrr::map_df(all_path_effects_posterior, ~{
#   tibble(mean = mean(.x, na.rm = TRUE), lower = quantile(.x, 0.055, na.rm = TRUE), upper = quantile(.x, 0.945, na.rm = TRUE))
# }, .id = "path") %>%
#   mutate(type = "Posterior")
# 
# 
# # ==============================================================================
# # 6. GENERATE INDUCED PRIORS FOR ALL PATHS (Simulation Method)
# # ==============================================================================
# # Instead of multiplying means, we simulate the full causal chain using 
# # random draws from the Prior Distributions. This correctly captures 
# # how uncertainty propagates through mediation and quadratic terms.
# 
# n_prior_draws <- 1000 # Enough to get a smooth distribution
# 
# # --- 6a. Define the Priors (MUST MATCH YOUR STAN MODEL) ---
# # You need to fill these in with the actual priors you used in your Stan code.
# # I have populated them based on your simulation code snippet.
# 
# prior_draws_df <- tibble(
#   .chain = 1,
#   .iteration = 1:n_prior_draws,
#   .draw = 1:n_prior_draws
# )
# 
# # Helper for normal priors
# rnorm_vec <- function(n, mean, sd) rnorm(n, mean, sd)
# 
# # --- Generate Priors for Main Goby Model ---
# prior_draws_df$beta_Rain <- rnorm_vec(n_prior_draws, 0.5, 0.5)
# prior_draws_df$beta_Micro <- rnorm_vec(n_prior_draws, 0.2, 0.25)
# prior_draws_df$beta_Year <- rnorm_vec(n_prior_draws, 0, 0.5)
# prior_draws_df$beta_Year_2 <- rnorm_vec(n_prior_draws, 0, 0.5)
# prior_draws_df$beta_Temp_2 <- rnorm_vec(n_prior_draws, -0.1, 0.25)
# prior_draws_df$beta_BreachDays_2 <- rnorm_vec(n_prior_draws, -0.1, 0.25)
# prior_draws_df$beta_SAV_2 <- rnorm_vec(n_prior_draws, -0.1, 0.25)
# prior_draws_df$beta_Goby_lag <- rnorm_vec(n_prior_draws, 0.1, 0.25)
# prior_draws_df$beta_SC_count <- rnorm_vec(n_prior_draws, -0.5, 0.5)
# prior_draws_df$beta_SB_count <- rnorm_vec(n_prior_draws, 0, 0.5)
# 
# # --- Generate Priors for Vectorized Parameters (Submodels + Goby) ---
# # NOTE: You must check your Stan code to see which index [1], [2], etc. 
# # corresponds to which variable. I am inferring from your DAG structure.
# 
# # DO Vector: [1]=Goby, [2]=SB, [3]=SC, [4]=SAV
# prior_draws_df$`beta_DO_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5) # DO -> Goby
# prior_draws_df$`beta_DO_vec[2]` <- rnorm_vec(n_prior_draws, 0, 0.5) # DO -> SB
# prior_draws_df$`beta_DO_vec[3]` <- rnorm_vec(n_prior_draws, 0, 0.5) # DO -> SC
# prior_draws_df$`beta_DO_vec[4]` <- rnorm_vec(n_prior_draws, 0.6, 0.5) # DO -> SAV
# 
# # SAV Vector: [1]=Goby, [2]=SB, [3]=SC
# prior_draws_df$`beta_SAV_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5) # SAV -> Goby
# prior_draws_df$`beta_SAV_vec[2]` <- rnorm_vec(n_prior_draws, 0, 0.5) # SAV -> SB
# prior_draws_df$`beta_SAV_vec[3]` <- rnorm_vec(n_prior_draws, 0, 0.5) # SAV -> SC
# 
# # Temp Vector: [1]=Goby, [2]=DO, [3]=SAV
# prior_draws_df$`beta_Temp_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Temp -> Goby
# prior_draws_df$`beta_Temp_vec[2]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Temp -> DO
# prior_draws_df$`beta_Temp_vec[3]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Temp -> SAV
# 
# # BreachDays Vector: [1]=Goby, [2]=Temp
# prior_draws_df$`beta_BreachDays_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Breach -> Goby
# prior_draws_df$`beta_BreachDays_vec[2]` <- rnorm_vec(n_prior_draws, 0.3, 0.5) # Breach -> Temp
# 
# # Wind Vector: [1]=DO, [2]=Temp
# prior_draws_df$`beta_Wind_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Wind -> DO
# prior_draws_df$`beta_Wind_vec[2]` <- rnorm_vec(n_prior_draws, 0, 0.5) # Wind -> Temp
# 
# # Substrate Vector: [1]=Goby, [2]=SC
# prior_draws_df$`beta_Substrate_vec[1]` <- rnorm_vec(n_prior_draws, 0, 0.5)
# prior_draws_df$`beta_Substrate_vec[2]` <- rnorm_vec(n_prior_draws, 0, 0.5)
# 
# # --- Generate Intercepts & Noise ---
# # These are needed for the function to run, even if they cancel out in the 
# # difference calculation. We sample them to be safe.
# for(z in 1:3) { # Assuming 3 zones
#   prior_draws_df[[paste0("a_Goby[", z, "]")]] <- rnorm(n_prior_draws, 2.6, 1)
#   prior_draws_df[[paste0("U[", z, "]")]] <- rnorm(n_prior_draws, 0, 1)
# }
# prior_draws_df$a_BreachDays <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$a_Temp <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$a_DO <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$a_SAV <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$a_SC <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$a_SB <- rnorm(n_prior_draws, 0, 1)
# prior_draws_df$k_U <- rnorm(n_prior_draws, 0.1, 0.1) 
# prior_draws_df$k_U_bio <- rnorm(n_prior_draws, 0.1, 0.1)
# 
# 
# # --- 6b. Run the Path Calculation on the PRIOR Data ---
# all_path_effects_prior <- list()
# 
# for (path_name in names(causal_paths)) {
#   current_path <- causal_paths[[path_name]]
#   
#   # Check if we have SD for the start node
#   start_node_for_sd <- current_path[1]
#   if (!(start_node_for_sd %in% c("SB_count", "SC_count")) && !(start_node_for_sd %in% names(sds))) {
#     next 
#   }
#   
#   # !!! HERE IS THE TRICK !!!
#   # We call the exact same function, but pass 'prior_draws_df' 
#   # instead of 'relevant_params' (the posterior).
#   all_path_effects_prior[[path_name]] <- calculate_path_effect(
#     path_name, 
#     current_path, 
#     draws_df = prior_draws_df,  # <--- PASSING PRIORS HERE
#     data_df = df.data, 
#     sds_df = sds, 
#     n_draws = n_prior_draws
#   )
# }
# 
# # --- 6c. Summarize Priors ---
# priors_summary <- purrr::map_df(all_path_effects_prior, ~{
#   tibble(
#     mean = mean(.x, na.rm = TRUE), 
#     lower = quantile(.x, 0.055, na.rm = TRUE), 
#     upper = quantile(.x, 0.945, na.rm = TRUE)
#   )
# }, .id = "path") %>%
#   mutate(type = "Prior")
# 
# 
# # ... (Continue with your existing plot code) ...
# 
# 
# # # --- 6. Generate and Summarize PRIORS for DIRECT paths (ON BETA * deltaX SCALE) ---
# # n_prior_draws <- 400 # was 4000
# # prior_effects_list <- list()
# # 
# # # Define the direct effect parameters and their specific prior distributions
# # direct_goby_params_priors <- list(
# #   "Micro" = list(param_name = "beta_Micro", prior_mean = 0.25, prior_sd = 0.25),
# #   "Substrate" = list(param_name = "beta_Substrate_vec[1]", prior_mean = 0.25, prior_sd = 0.25), # From b_bar_Substrate
# #   "Year" = list(param_name = "beta_Year", prior_mean = 0, prior_sd = 0.5),
# #   "Year_2" = list(param_name = "beta_Year_2", prior_mean = 0, prior_sd = 0.5),
# #   "Temp_2" = list(param_name = "beta_Temp_2", prior_mean = -0.1, prior_sd = 0.25),
# #   "BreachDays_2" = list(param_name = "beta_BreachDays_2", prior_mean = -0.10, prior_sd = 0.25),
# #   "Goby_lag" = list(param_name = "beta_Goby_lag", prior_mean = 0.25, prior_sd = 0.25),
# #   "SAV" = list(param_name = "beta_SAV_vec[1]", prior_mean = 0.00, prior_sd = 0.25), # From b_bar_SAV
# #   "SAV_2" = list(param_name = "beta_SAV_2", prior_mean = -0.1, prior_sd = 0.25),
# #   "DO" = list(param_name = "beta_DO_vec[1]", prior_mean = 0.25, prior_sd = 0.25), # From b_bar_DO
# #   "BreachDays" = list(param_name = "beta_BreachDays_vec[1]", prior_mean = 0.25, prior_sd = 0.25), # From b_bar_Breach
# #   "Temp" = list(param_name = "beta_Temp_vec[1]", prior_mean = 0, prior_sd = 0.5), # From b_bar_Temp
# #   "SC_count" = list(param_name = "beta_SC_count", prior_mean = -0.10, prior_sd = 0.25),
# #   "SB_count" = list(param_name = "beta_SB_count", prior_mean = 0, prior_sd = 0.5)
# # )
# # 
# # # Iterate through direct paths to generate prior predictive effects
# # for (path_start_node in names(direct_goby_params_priors)) {
# #   path_name <- paste0(path_start_node, " -> Goby")
# #   
# #   if (path_name %in% names(causal_paths)) {
# #     prior_info <- direct_goby_params_priors[[path_start_node]]
# #     
# #     # Generate draws directly from the prior for beta
# #     prior_draws_beta <- rnorm(n_prior_draws, mean = prior_info$prior_mean, sd = prior_info$prior_sd)
# #     
# #     # --- CHANGE: Scale prior draws by delta_X BEFORE summarizing ---
# #     delta_X <- if (path_start_node %in% c("SC_count", "SB_count")) 1 else sds[[1, path_start_node]]
# #     if (is.null(delta_X) || is.na(delta_X)) {
# #       message(paste("Warning: SD for", path_start_node, "not found. Assuming delta_X = 1 for prior scaling."))
# #       delta_X <- 1
# #     }
# #     scaled_prior_dist <- prior_draws_beta * delta_X # Now scaling again
# #     
# #     # Calculate summary stats on the SCALED prior distribution
# #     prior_effects_list[[path_name]] <- tibble(
# #       path = path_name,
# #       mean = mean(scaled_prior_dist),
# #       lower = quantile(scaled_prior_dist, 0.055), # 89% CI for beta * deltaX
# #       upper = quantile(scaled_prior_dist, 0.945), # 89% CI for beta * deltaX
# #       type = "Prior"
# #     )
# #   }
# # }
# # priors_summary <- bind_rows(prior_effects_list)
# 
# # --- 7. Combine Posterior and Prior Summaries and Finalize Labels ---
# # combines scaled posteriors and scaled priors
# combined_summary <- bind_rows(
#   effects_summary,
#   priors_summary
# ) %>%
#   mutate(
#     path_label = stringr::str_replace_all(path, " -> ", " \u2192 "),
#     path_label = stringr::str_replace_all(path_label, "SB_count", "SB"),
#     path_label = stringr::str_replace_all(path_label, "SC_count", "SC"),
#     path_label = stringr::str_replace_all(path_label, "_2", "\u00B2")
#   )
# 
# combined_summary <- combined_summary %>%
#   mutate(
#     significant_color = case_when(
#       type == "Posterior" & lower > 0 & upper > 0 ~ "Positive",
#       type == "Posterior" & lower < 0 & upper < 0 ~ "Negative",
#       type == "Posterior" ~ "Overlap Zero",
#       TRUE ~ NA_character_
#     )
#   )
# 
# path_order <- combined_summary %>%
#   filter(type == "Posterior") %>%
#   arrange(mean) %>%
#   pull(path_label) %>%
#   unique()
# 
# combined_summary$path_label <- factor(combined_summary$path_label, levels = path_order)
# 
# combined_summary_for_plot <- combined_summary %>%
#   filter(!is.na(path_label))
# 
# # --- 8. Create the Combined Forest Plot  ---
# causal_path_forest_plot <- ggplot(
#   combined_summary_for_plot,
#   aes(x = path_label, y = mean, ymin = lower, ymax = upper)
# ) +
#   # Vertical line at zero effect (after coord_flip)
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
#   # Pointranges for Prior and Posterior
#   geom_pointrange(
#     aes(
#       shape = type,
#       color = significant_color,
#       size = type,
#       linewidth = type
#     ),
#     position = position_nudge(x = ifelse(combined_summary_for_plot$type == "Posterior", 0.1, -0.1))
#   ) +
#   coord_flip() +
#   
#   # --- coord_cartesian(ylim = ...) REMOVED for now ---
#   
#   # Scales
#   scale_shape_manual(
#     name = NULL,
#     values = c("Prior" = 1, "Posterior" = 16),
#     labels = c("Prior", "Posterior")
#   ) +
#   scale_color_manual(
#     name = NULL,
#     values = c("Positive" = "#0072B2", "Negative" = "#D55E00", "Overlap Zero" = "black"),
#     breaks = c("Positive", "Negative", "Overlap Zero"),
#     labels = c("Positive Effect", "Negative Effect", "Overlaps Zero"),
#     na.value = "grey70"
#   ) +
#   scale_size_manual(
#     name = NULL,
#     values = c("Prior" = 0.6, "Posterior" = 1.0),
#     guide = "none"
#   ) +
#   scale_linewidth_manual(
#     name = NULL,
#     values = c("Prior" = 0.6, "Posterior" = 1.0),
#     guide = "none"
#   ) +
#   labs(
#     x = "Causal Pathway", # Label for the VERTICAL axis (after flip)
#     y = "log Δ in Goby Density per SD (β × ΔX)" # Label for the HORIZONTAL axis (after flip)
#   ) +
#   theme_minimal(base_size = 15) +
#   theme(
#     plot.title.position = "plot",
#     axis.text.y = element_text(size = 10), # Path labels on vertical axis
#     axis.text.x = element_text(size = 10), # Effect sizes on horizontal axis
#     legend.position = c(0.22, 0.45),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.5),
#     legend.margin = margin(t = 2, r = 5, b = 5, l = 5, unit = "pt"),
#     legend.key.size = unit(0.8, "lines"),
#     legend.box = "vertical"
#   ) +
#   guides(
#     color = guide_legend(order = 1, override.aes = list(shape = 16)),
#     shape = guide_legend(order = 2, override.aes = list(color = c("grey70", "black")))
#   )
# 
# print(causal_path_forest_plot)
# 
# # Save the plot
# ggsave("Output/Plots/causal_path_forest_plot_with_priors.png", 
#        plot = causal_path_forest_plot, 
#        width = 22, height = 30, units = "cm", dpi = 300,
#        bg = "white")
