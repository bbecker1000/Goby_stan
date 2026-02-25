# ZONE-LEVEL LATENT ENVIRONMENT MODEL SIMULATION ANALYSIS
# Updated to match zone_latent_environment_SEM.stan
# Includes measurement error correction
# 2026-01-12

library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(forcats)
library(stringr)
library(tidyr)
library(patchwork)

# Fix namespace conflicts
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# ============================================================================
# LOAD SIMULATED DATA AND TRUE PARAMETERS
# ============================================================================

# Load the simulated data with measurement error
load("simulated_data_zone_latent.RData")  # Loads sim_data, true_params, true_total_effects

cat("=== DATA LOADED ===\n")
cat("Sample size:", nrow(sim_data), "\n")
cat("Zones:", length(unique(sim_data$Zone)), "\n")
cat("Measurement error correlation (Temp):", 
    round(cor(sim_data$Temp_measured, sim_data$Temp_true), 3), "\n")
cat("Measurement error correlation (DO):", 
    round(cor(sim_data$DO_measured, sim_data$DO_true), 3), "\n\n")

# ============================================================================
# HELPER FUNCTIONS FOR CAUSAL INFERENCE - ZONE LATENT MODEL
# ============================================================================

# Function to calculate mu_goby deterministically - ZONE LATENT MODEL
calculate_mu_goby_full_dag <- function(draw_row, input_data, 
                                       intervention_var = NULL, 
                                       intervention_value = NULL) {
  
  current_values <- as.list(input_data)
  
  # Helper to get parameter value
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
  
  # Get zone-specific latent variable and TRUE environment
  zone_idx <- input_data$Zone
  U_phys_col <- paste0("U_phys[", zone_idx, "]")
  Temp_true_col <- paste0("Temp_true_zone[", zone_idx, "]")
  DO_true_col <- paste0("DO_true_zone[", zone_idx, "]")
  
  # If there's an intervention, set it
  if (!is.null(intervention_var)) {
    current_values[[intervention_var]] <- intervention_value
    
    # Update quadratic terms
    if (intervention_var == "Year") {
      current_values$Year_2 <- intervention_value^2
    }
  }
  
  # CONDITIONAL RECALCULATION based on intervention
  no_intervention <- is.null(intervention_var)
  
  # Recalculate downstream variables
  recalc_breach <- !no_intervention && (intervention_var %in% c("Rain"))
  recalc_temp <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind"))
  recalc_do <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp_true"))
  
  # BreachDays (affected by Rain)
  if (!is.na(recalc_breach) && recalc_breach && !identical(intervention_var, "BreachDays")) {
    Breach_nu <- get_param("a_BreachDays") + 
      get_param("beta_Rain_Breach") * current_values$Rain + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$BreachDays <- Breach_nu
    current_values$BreachDays_2 <- Breach_nu^2
  } else if (identical(intervention_var, "BreachDays")) {
    current_values$BreachDays_2 <- current_values$BreachDays^2
  }
  
  # TRUE Temp (zone-level, affected by BreachDays, Wind)
  if (!is.na(recalc_temp) && recalc_temp && !identical(intervention_var, "Temp_true")) {
    Temp_nu <- get_param("a_Temp") + 
      get_param("beta_Breach_Temp") * current_values$BreachDays + 
      get_param("beta_Wind_Temp") * current_values$Wind + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$Temp_true <- Temp_nu
    current_values$Temp_2 <- Temp_nu^2
  } else if (identical(intervention_var, "Temp_true")) {
    current_values$Temp_2 <- current_values$Temp_true^2
  } else {
    # Use zone-level true temp from model
    current_values$Temp_true <- get_param(Temp_true_col)
    current_values$Temp_2 <- current_values$Temp_true^2
  }
  
  # TRUE DO (zone-level, affected by Temp, Wind)
  if (!is.na(recalc_do) && recalc_do && !identical(intervention_var, "DO_true")) {
    DO_nu <- get_param("a_DO") + 
      get_param("beta_Temp_DO") * current_values$Temp_true + 
      get_param("beta_Wind_DO") * current_values$Wind + 
      get_param("k_U_phys") * get_param(U_phys_col)
    current_values$DO_true <- DO_nu
  } else if (!identical(intervention_var, "DO_true")) {
    # Use zone-level true DO from model
    current_values$DO_true <- get_param(DO_true_col)
  }
  
  # SAV: measured (never recalculate)
  if (identical(intervention_var, "SAV")) {
    current_values$SAV_2 <- current_values$SAV^2
  }
  
  # Final Goby calculation (responds to TRUE environment)
  a_Goby_col <- paste0("a_Goby[", current_values$Zone, "]")
  mu_goby_val <- get_param(a_Goby_col) +
    get_param("beta_SAV_Goby") * current_values$SAV +
    get_param("beta_SAV_2") * current_values$SAV_2 + 
    get_param("beta_DO_Goby") * current_values$DO_true +  # TRUE DO
    get_param("beta_BreachDays_Goby") * current_values$BreachDays +
    get_param("beta_BreachDays_2") * current_values$BreachDays_2 +
    get_param("beta_Temp_Goby") * current_values$Temp_true +  # TRUE Temp
    get_param("beta_Temp_2") * current_values$Temp_2 +
    get_param("beta_SC_count") * current_values$SC_count +
    get_param("beta_SB_count") * current_values$SB_count +
    get_param("beta_Substrate_Goby") * current_values$Substrate +  # FIXED
    get_param("beta_Micro") * current_values$Micro +
    get_param("beta_Year") * current_values$Year +
    get_param("beta_Year_2") * current_values$Year_2 +
    get_param("beta_Goby_lag") * current_values$Goby_lag +
    current_values$Area
  
  return(mu_goby_val)
}

# Function to calculate total effects
calculate_total_effect <- function(var_name, draws_df, data_df, sds_df, n_draws = 100, n_obs = 100) {
  
  message(paste("  Calculating total effect for:", var_name))
  
  # SPECIAL HANDLING FOR YEAR
  if (var_name == "Year") {
    message("    Special handling: Year has quadratic term, reporting marginal effect at mean")
    
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
          if (length(val) == 0) {
            stop(paste("Parameter not found:", param_name))
          }
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
  
  # ORIGINAL CODE FOR ALL OTHER VARIABLES
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
# PART 1: FIT THE ZONE LATENT SCM MODEL ON SIMULATED DATA
# ============================================================================

# Prepare data - use MEASURED Temp and DO
sim_data_for_stan <- sim_data
sim_data_for_stan$Temp <- sim_data$Temp_measured  # Model sees measured
sim_data_for_stan$DO <- sim_data$DO_measured      # Model sees measured

list_df_sim <- lapply(as.list(sim_data_for_stan), function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(as.vector(x))
  } else {
    return(x)
  }
})

# Remove the extra columns we don't need for Stan
list_df_sim$Temp_measured <- NULL
list_df_sim$Temp_true <- NULL
list_df_sim$DO_measured <- NULL
list_df_sim$DO_true <- NULL
list_df_sim$Goby_Density <- NULL

list_df_sim$N <- length(list_df_sim$Goby)
list_df_sim$J <- length(unique(list_df_sim$Zone))

# Rename for Stan's expected variable names
list_df_sim$Temp_measured <- list_df_sim$Temp
list_df_sim$DO_measured <- list_df_sim$DO
list_df_sim$Temp <- NULL
list_df_sim$DO <- NULL

cat("\n=== Fitting Zone Latent Environment SCM ===\n")

# Run the zone latent SCM model
mod.SEM.sim <- mod.SEM$sample(
  data = list_df_sim,
  seed = 123,
  chains = 3,
  iter_warmup = 6000,
  iter_sampling = 1000,
  parallel_chains = 3,
  refresh = 100,
  adapt_delta = 0.97,
  max_treedepth = 13,
  output_dir = "stan_output",
  save_warmup = FALSE
)

saveRDS(mod.SEM.sim, file = "Output/Models/mod.SEM.zone_latent.sim.rds")
mod.SEM.sim$diagnostic_summary()
print(mod.SEM.sim$summary(), n = 150)

# Check measurement error estimates
cat("\n=== MEASUREMENT ERROR RECOVERY ===\n")
meas_error_summary <- mod.SEM.sim$summary(variables = c("measurement_error_Temp", 
                                                        "measurement_error_DO",
                                                        "within_zone_var_Temp",
                                                        "within_zone_var_DO"))
print(meas_error_summary)

cat("\nTrue values:\n")
cat("sigma_Temp_measure:", true_params$sigma_Temp_measure, "\n")
cat("sigma_DO_measure:", true_params$sigma_DO_measure, "\n")
cat("sigma_Temp_within:", true_params$sigma_Temp_within, "\n")
cat("sigma_DO_within:", true_params$sigma_DO_within, "\n\n")

# ============================================================================
# PART 2: EXTRACT DIRECT EFFECT COEFFICIENTS
# ============================================================================

coefficient_summary.sim <- summarise_draws(
  mod.SEM.sim,
  "mean",
  ~quantile(.x, probs = c(0.045, 0.955))
)

# SCM betas
beta.joint.sim <- coefficient_summary.sim %>%
  dplyr::filter(str_starts(variable, "beta_")) %>%
  dplyr::mutate(variable = str_remove(variable, "^.*?_"))
beta.joint.sim$model <- "SCM"
beta.joint.sim$effect_type <- "Direct"

# ============================================================================
# PART 3: FIT GLMM FOR COMPARISON (uses MEASURED Temp/DO)
# ============================================================================

cat("\n=== Fitting GLMM for comparison ===\n")

library(brms)

# GLMM uses MEASURED Temp and DO (suffers from attenuation)
model.brms.sim <- brm(
  Goby ~ Rain + Wind + BreachDays + BreachDays_2 + Temp_measured + Temp_2 + 
    DO_measured + SAV + SAV_2 + SC_count + SB_count + Substrate + Micro + 
    Year + Year_2 + Goby_lag + offset(Area) + (1 | Zone),
  data = sim_data,
  family = negbinomial(),
  prior = c(
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sd),
    prior(gamma(0.01, 0.01), class = shape)
  ),
  chains = 3,
  iter = 4000,
  warmup = 3000,
  cores = 3,
  seed = 123,
  control = list(adapt_delta = 0.97)
)

summary(model.brms.sim)
saveRDS(model.brms.sim, file = "Output/Models/model.brms.zone_latent.sim.rds")

# GLMM betas
draws_glmm <- as_draws_df(model.brms.sim)

coefficient_summary.sim.brm <- summarise_draws(
  model.brms.sim,
  "mean",
  ~quantile(.x, probs = c(0.045, 0.955))
)

beta.sim.brm <- coefficient_summary.sim.brm %>%
  dplyr::filter(str_starts(variable, "b_")) %>%
  dplyr::filter(variable != "b_Intercept") %>%
  dplyr::mutate(variable = str_remove(variable, "^.*?_"))
beta.sim.brm$model <- "GLMM"
beta.sim.brm$effect_type <- "Direct"

# ============================================================================
# PART 4: TRUE DIRECT EFFECTS
# ============================================================================

betas_true_goby <- c(
  # Pathway coefficients
  Rain_Breach = true_params$beta_Rain_Breach,
  Breach_Temp = true_params$beta_Breach_Temp,
  Temp_DO = true_params$beta_Temp_DO,
  Wind_Temp = true_params$beta_Wind_Temp,
  Wind_DO = true_params$beta_Wind_DO,
  
  # Direct effects on Goby
  BreachDays_Goby = true_params$beta_BreachDays_Goby,
  BreachDays_2 = true_params$beta_BreachDays_2,
  Temp_Goby = true_params$beta_Temp_Goby,
  Temp_2 = true_params$beta_Temp_2,
  DO_Goby = true_params$beta_DO_Goby,
  SAV_Goby = true_params$beta_SAV_Goby,
  SAV_2 = true_params$beta_SAV_2,
  SC_count = true_params$beta_SC_count,
  SB_count = true_params$beta_SB_count,
  Substrate_Goby = true_params$beta_Substrate,  # Store with correct name
  Micro = true_params$beta_Micro,
  Year = true_params$beta_Year,
  Year_2 = true_params$beta_Year_2,
  Goby_lag = true_params$beta_Goby_lag
)

df.true.direct <- tibble(
  variable = names(betas_true_goby),
  mean = betas_true_goby,
  `4.5%` = NA,
  `95.5%` = NA,
  model = "True",
  effect_type = "Direct"
)

# ============================================================================
# PART 5: TRUE TOTAL EFFECTS
# ============================================================================

cat("\n=== Using pre-calculated TRUE total effects ===\n")
print(true_total_effects)

true_total_effects_df <- true_total_effects %>%
  mutate(
    variable = paste0(Variable, " (Total)"),
    mean = Total_Effect,
    `4.5%` = NA,
    `95.5%` = NA,
    model = "True",
    effect_type = "Total"
  ) %>%
  select(variable, mean, `4.5%`, `95.5%`, model, effect_type)

# ============================================================================
# PART 6: CALCULATE ESTIMATED TOTAL EFFECTS
# ============================================================================

cat("\n=== Calculating SCM Total Effects ===\n")

# Calculate SDs (use measured values for interventions)
sim_data$Goby_Density <- sim_data$Goby / exp(sim_data$Area)
numeric_cols <- names(sim_data)[sapply(sim_data, is.numeric)]
sds <- data.frame(
  purrr::map_dfc(sim_data[numeric_cols], ~ sd(.x, na.rm = TRUE))
)
colnames(sds) <- numeric_cols

# Variables to test
total_effect_vars <- c("Rain", "Wind", "Substrate", "Micro", "Year",
                       "BreachDays", "SAV")

# Extract SCM draws
draws_scm <- mod.SEM.sim$draws(format = "df")

# Calculate SCM total effects
scm_total_effects_list <- purrr::map(
  total_effect_vars,
  function(var_name) {
    cat(sprintf("\n=== Processing %s ===\n", var_name))
    tryCatch({
      effect_dist <- calculate_total_effect(
        var_name, 
        draws_scm, 
        sim_data, 
        sds, 
        n_draws = 100,
        n_obs = 100
      )
      
      cat(sprintf("  SUCCESS: Mean effect = %.4f\n", mean(effect_dist, na.rm = TRUE)))
      
      tibble(
        variable = paste0(var_name, " (Total)"),
        mean = mean(effect_dist, na.rm = TRUE),
        `5.5%` = quantile(effect_dist, 0.055, na.rm = TRUE),
        `94.5%` = quantile(effect_dist, 0.945, na.rm = TRUE),
        model = "SCM",
        effect_type = "Total"
      )
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      cat(sprintf("  Full error: %s\n", toString(e)))
      return(NULL)
    })
  }
)

scm_total_effects <- bind_rows(scm_total_effects_list[!sapply(scm_total_effects_list, is.null)])

# GLMM total effects
glmm_total_effects_list <- purrr::map(
  total_effect_vars,
  function(var_name) {
    coef_name <- paste0("b_", var_name)
    
    if (coef_name %in% names(draws_glmm)) {
      coef_draws <- draws_glmm[[coef_name]]
      
      tibble(
        variable = paste0(var_name, " (Total)"),
        mean = mean(coef_draws, na.rm = TRUE),
        `4.5%` = quantile(coef_draws, 0.045, na.rm = TRUE),
        `95.5%` = quantile(coef_draws, 0.955, na.rm = TRUE),
        model = "GLMM",
        effect_type = "Total"
      )
    } else {
      return(NULL)
    }
  }
)

glmm_total_effects <- bind_rows(glmm_total_effects_list[!sapply(glmm_total_effects_list, is.null)])

# ============================================================================
# PART 7: COMBINE AND PLOT
# ============================================================================

# Combine direct effects
df.direct <- bind_rows(
  df.true.direct %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`),
  beta.joint.sim %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`),
  beta.sim.brm %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)
)

# Combine total effects
true_total_effects_renamed <- true_total_effects_df %>% 
  dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)

scm_total_effects_renamed <- scm_total_effects %>% 
  dplyr::rename(lower_4.5 = `5.5%`, upper_95.5 = `94.5%`)

glmm_total_effects_renamed <- glmm_total_effects %>% 
  dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)

# DEBUG: Check what we have before combining
cat("\n=== DEBUGGING TOTAL EFFECTS DATA ===\n")
cat("True effects:\n")
print(true_total_effects_renamed)
cat("\nSCM effects:\n")
print(scm_total_effects_renamed)
cat("\nGLMM effects:\n")
print(glmm_total_effects_renamed)

df.total <- bind_rows(
  true_total_effects_renamed,
  scm_total_effects_renamed,
  glmm_total_effects_renamed
)

cat("\nCombined total effects:\n")
print(df.total)
cat("\nUnique variables:", unique(df.total$variable), "\n")
cat("Unique models:", unique(df.total$model), "\n")
cat("=================================\n\n")

# Relabel for plotting
df.direct_relabeled <- df.direct %>%
  mutate(
    variable = str_replace(variable, "_2$", "²"),
    variable = str_replace(variable, "^Goby_lag$", "Goby lag"),
    variable = str_replace(variable, "Rain_Breach", "Rain → BreachDays"),
    variable = str_replace(variable, "Breach_Temp", "BreachDays → Temp"),
    variable = str_replace(variable, "Temp_DO", "Temp → DO"),
    variable = str_replace(variable, "Wind_DO", "Wind → DO"),
    variable = str_replace(variable, "Wind_Temp", "Wind → Temp")
  )

# Create plots
df_direct_plot <- df.direct_relabeled %>% mutate(variable = factor(variable))

p.direct <- df_direct_plot %>%
  ggplot(aes(x = fct_rev(variable), y = mean, color = model)) +
  geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), 
                  size = 0.7, position = position_dodge(width = 0.7)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("True" = "black", "SCM" = "red3", "GLMM" = "blue")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.85, 0.96),
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_text(face = "bold", size = 20)
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(title = "A. Direct Effects on Goby (Zone Latent Model)", 
       y = "Beta estimate", x = "")

# Total effects plot
df_total_plot <- df.total %>% mutate(variable = factor(variable))

p.total <- df_total_plot %>%
  ggplot(aes(x = fct_rev(variable), y = mean, color = model)) +
  geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), 
                  size = 0.7, position = position_dodge(width = 0.7)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("True" = "black", "SCM" = "red3", "GLMM" = "blue"),
                     breaks = c("True", "SCM", "GLMM")) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.85, 0.96),
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_text(face = "bold", size = 20)
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(title = "B. Total Causal Effects on Goby (Zone Latent Model)", 
       subtitle = "With measurement error correction (~44 params, 7.1 obs/param; GLMM suffers from attenuation)", 
       y = "Total effect (log scale)", x = "")

p.combined <- p.direct / p.total + plot_layout(heights = c(2, 1))
print(p.combined)

# Save plots
ggsave("Output/Plots/sim.forest.combined.zone_latent.png", 
       plot = p.combined, width = 25, height = 35, units = "cm", dpi = 300)
ggsave("Output/Plots/sim.forest.direct.zone_latent.png", 
       plot = p.direct, width = 20, height = 30, units = "cm")
ggsave("Output/Plots/sim.forest.total.zone_latent.png", 
       plot = p.total, width = 20, height = 15, units = "cm")

# Calculate coverage
df.all <- bind_rows(df.direct_relabeled, df.total)
coverage_summary <- df.all %>%
  dplyr::filter(model != "True") %>%
  left_join(df.all %>% dplyr::filter(model == "True") %>% 
              dplyr::select(variable, effect_type, true_mean = mean),
            by = c("variable", "effect_type")) %>%
  mutate(covered = (true_mean >= lower_4.5) & (true_mean <= upper_95.5),
         abs_error = abs(mean - true_mean),
         rel_error = abs_error / abs(true_mean))

coverage_by_model <- coverage_summary %>%
  group_by(model, effect_type) %>%
  summarise(n_effects = n(), 
            coverage_rate = mean(covered, na.rm = TRUE),
            mean_abs_error = mean(abs_error, na.rm = TRUE),
            mean_rel_error = mean(rel_error, na.rm = TRUE),
            .groups = "drop")

cat("\n=== PARAMETER RECOVERY SUMMARY ===\n")
print(coverage_by_model, n = 20)

cat("\n=== TOTAL EFFECT COMPARISON ===\n")
total_comparison <- df.total %>%
  dplyr::select(variable, model, mean, lower_4.5, upper_95.5) %>%
  arrange(variable, model)
print(total_comparison, n = 30)

cat("\n=== KEY INSIGHT: MEASUREMENT ERROR CORRECTION ===\n")
cat("Zone Latent Model (~44 params, 7.1 obs/param):\n")
cat("  - Models TRUE zone-level Temp and DO\n")
cat("  - Measured Temp/DO are noisy observations\n")
cat("  - Strong informative priors prevent attenuation\n")
cat("  - SC/SB treated as measured (no submodels)\n\n")

cat("Expected benefits over maximal simplification:\n")
cat("  1. Preserves environmental cascade (no omitted variable bias)\n")
cat("  2. Corrects for measurement error\n")
cat("  3. Should recover TRUE effects despite noisy measurements\n")
cat("  4. GLMM should show attenuation (smaller effects due to noise)\n\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Plots saved to Output/Plots/\n\n")

#old below--------


# # SIMULATION COMPARISON: TOTAL EFFECTS ONLY
# # Compares True vs SCM vs GLMM for Direct Effects + Total Causal Effects
# # Much cleaner and more accurate than individual mediated paths!
# 
# 
# 
# #2026-01-11
# # SIMULATION COMPARISON: TOTAL EFFECTS ONLY
# # Compares True vs SCM vs GLMM for Direct Effects + Total Causal Effects
# # Much cleaner and more accurate than individual mediated paths!
# # UPDATED: Fixed Year effect calculation to report marginal effect at mean
# 
# #2026-01-09 (Updated 2026-01-11)
# 
# # SIMULATION COMPARISON: TOTAL EFFECTS ONLY
# # Compares True vs SCM vs GLMM for Direct Effects + Total Causal Effects
# # Much cleaner and more accurate than individual mediated paths!
# # UPDATED: Fixed Year effect calculation to report marginal effect at mean
# 
# #(Updated 2026-01-11)
# 
# library(cmdstanr)
# library(posterior)
# library(dplyr)
# library(ggplot2)
# library(purrr)
# library(tibble)
# library(forcats)
# library(stringr)
# library(tidyr)
# library(patchwork)
# 
# # Fix namespace conflicts
# select <- dplyr::select
# rename <- dplyr::rename
# filter <- dplyr::filter
# 
# # ============================================================================
# # LOAD SIMULATED DATA AND TRUE PARAMETERS
# # ============================================================================
# 
# # Load the simulated data and true parameters from your simulation
# load("simulated_data_with_latents.RData")  # Loads sim_data and true_params
# 
# # ============================================================================
# # HELPER FUNCTIONS FOR CAUSAL INFERENCE
# # ============================================================================
# 
# # Function to calculate mu_goby deterministically - FULL DAG propagation
# calculate_mu_goby_full_dag <- function(draw_row, input_data, 
#                                        intervention_var = NULL, 
#                                        intervention_value = NULL) {
#   
#   current_values <- as.list(input_data)
#   
#   # Helper to get parameter value - handles both dataframe row and named vector
#   get_param <- function(param_name) {
#     if (is.data.frame(draw_row)) {
#       val <- draw_row[[param_name]]
#       if (length(val) == 0) {
#         stop(paste("Parameter not found:", param_name))
#       }
#       return(val[1])
#     } else {
#       return(draw_row[[param_name]])
#     }
#   }
#   
#   # Get zone-specific latent variables
#   zone_idx <- input_data$Zone
#   U_bio_col <- paste0("U_bio[", zone_idx, "]")
#   U_phys_col <- paste0("U_phys[", zone_idx, "]")
#   
#   # If there's an intervention, set it
#   if (!is.null(intervention_var)) {
#     current_values[[intervention_var]] <- intervention_value
#     
#     # Update quadratic terms if intervening on variables with quadratics
#     if (intervention_var == "Year") {
#       current_values$Year_2 <- intervention_value^2
#     }
#     # Note: Micro, Wind, Rain, Substrate don't have quadratics in Goby equation
#   }
#   
#   # CONDITIONAL RECALCULATION based on intervention
#   # Only recalculate variables DOWNSTREAM of the intervention
#   
#   # Determine what to recalculate based on intervention
#   # IMPORTANT: For baseline (no intervention), don't recalculate ANYTHING - use actual data!
#   no_intervention <- is.null(intervention_var)
#   
#   # Only recalculate if there IS an intervention AND the variable is downstream
#   recalc_breach <- !no_intervention && (intervention_var %in% c("Rain"))
#   recalc_temp <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind"))
#   recalc_do <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp"))
#   recalc_sav <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO"))
#   recalc_sc <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO", "SAV", "Substrate"))
#   recalc_sb <- !no_intervention && (intervention_var %in% c("Rain", "BreachDays", "Wind", "Temp", "DO", "SAV"))
#   
#   # BreachDays (affected by Rain)
#   if (!is.na(recalc_breach) && recalc_breach && !identical(intervention_var, "BreachDays")) {
#     Breach_nu <- get_param("a_BreachDays") + 
#       get_param("beta_Rain") * current_values$Rain + 
#       get_param("k_U_phys") * get_param(U_phys_col)
#     current_values$BreachDays <- Breach_nu
#     current_values$BreachDays_2 <- Breach_nu^2
#   } else if (identical(intervention_var, "BreachDays")) {
#     # Intervention on BreachDays - update quadratic
#     current_values$BreachDays_2 <- current_values$BreachDays^2
#   }
#   
#   # Temp (affected by BreachDays, Wind)
#   if (!is.na(recalc_temp) && recalc_temp && !identical(intervention_var, "Temp")) {
#     Temp_nu <- get_param("a_Temp") + 
#       get_param("beta_BreachDays_vec[2]") * current_values$BreachDays + 
#       get_param("beta_Wind_vec[2]") * current_values$Wind + 
#       get_param("k_U_phys") * get_param(U_phys_col)
#     current_values$Temp <- Temp_nu
#     current_values$Temp_2 <- Temp_nu^2
#   } else if (identical(intervention_var, "Temp")) {
#     # Intervention on Temp - update quadratic
#     current_values$Temp_2 <- current_values$Temp^2
#   }
#   
#   # DO (affected by Temp, Wind)
#   if (!is.na(recalc_do) && recalc_do && !identical(intervention_var, "DO")) {
#     DO_nu <- get_param("a_DO") + 
#       get_param("beta_Temp_vec[2]") * current_values$Temp + 
#       get_param("beta_Wind_vec[1]") * current_values$Wind + 
#       get_param("k_U_phys") * get_param(U_phys_col)
#     current_values$DO <- DO_nu
#   }
#   
#   # SAV (affected by DO, Temp)
#   if (!is.na(recalc_sav) && recalc_sav && !identical(intervention_var, "SAV")) {
#     SAV_nu <- get_param("a_SAV") + 
#       get_param("beta_DO_vec[4]") * current_values$DO + 
#       get_param("beta_Temp_vec[3]") * current_values$Temp + 
#       get_param("k_U_bio") * get_param(U_bio_col)
#     current_values$SAV <- SAV_nu
#     current_values$SAV_2 <- SAV_nu^2
#   } else if (identical(intervention_var, "SAV")) {
#     # Intervention on SAV - update quadratic
#     current_values$SAV_2 <- current_values$SAV^2
#   }
#   
#   # SC_count (affected by Substrate, DO, SAV)
#   if (!is.na(recalc_sc) && recalc_sc && !identical(intervention_var, "SC_count")) {
#     SC_mu_logit <- get_param("a_SC") + 
#       get_param("beta_Substrate_vec[2]") * current_values$Substrate + 
#       get_param("beta_DO_vec[3]") * current_values$DO + 
#       get_param("beta_SAV_vec[3]") * current_values$SAV + 
#       get_param("k_U_bio") * get_param(U_bio_col)
#     current_values$SC_count <- plogis(SC_mu_logit)
#   }
#   
#   # SB_count (affected by DO, SAV)
#   if (!is.na(recalc_sb) && recalc_sb && !identical(intervention_var, "SB_count")) {
#     SB_mu_logit <- get_param("a_SB") + 
#       get_param("beta_DO_vec[2]") * current_values$DO + 
#       get_param("beta_SAV_vec[2]") * current_values$SAV + 
#       get_param("k_U_bio") * get_param(U_bio_col)
#     current_values$SB_count <- plogis(SB_mu_logit)
#   }
#   
#   # Final Goby calculation (same whether baseline or intervention)
#   a_Goby_col <- paste0("a_Goby[", current_values$Zone, "]")
#   mu_goby_val <- get_param(a_Goby_col) +
#     get_param("beta_SAV_vec[1]") * current_values$SAV +
#     get_param("beta_SAV_2") * current_values$SAV_2 + 
#     get_param("beta_DO_vec[1]") * current_values$DO +
#     get_param("beta_BreachDays_vec[1]") * current_values$BreachDays +
#     get_param("beta_Temp_vec[1]") * current_values$Temp +
#     get_param("beta_SC_count") * current_values$SC_count +
#     get_param("beta_SB_count") * current_values$SB_count +
#     get_param("beta_Substrate_vec[1]") * current_values$Substrate +
#     get_param("beta_Micro") * current_values$Micro +
#     get_param("beta_Year") * current_values$Year +
#     get_param("beta_Year_2") * current_values$Year_2 +
#     get_param("beta_Temp_2") * current_values$Temp_2 +
#     get_param("beta_BreachDays_2") * current_values$BreachDays_2 +
#     get_param("beta_Goby_lag") * current_values$Goby_lag +
#     current_values$Area
#   
#   return(mu_goby_val)
# }
# 
# # Function to calculate total effects - UPDATED VERSION WITH YEAR FIX
# calculate_total_effect <- function(var_name, draws_df, data_df, sds_df, n_draws = 100, n_obs = 100) {
#   
#   message(paste("  Calculating total effect for:", var_name))
#   
#   # SPECIAL HANDLING FOR YEAR
#   # For Year, we want the marginal effect at mean Year, not the effect of +1 SD intervention
#   # This is because Year has a quadratic term, so the effect depends on the baseline value
#   if (var_name == "Year") {
#     message("    Special handling: Year has quadratic term, reporting marginal effect at mean")
#     
#     # Select parameter draws
#     total_draws_available <- nrow(draws_df)
#     if (n_draws > total_draws_available) n_draws <- total_draws_available
#     
#     if (total_draws_available == 1) {
#       selected_draws <- draws_df
#       n_draws <- 1
#     } else {
#       sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
#       selected_draws <- draws_df[sample_indices, , drop = FALSE]
#     }
#     
#     # For centered Year (mean = 0), the marginal effect is:
#     # d(beta_Year * Year + beta_Year_2 * Year^2) / dYear |_{Year=0}
#     # = beta_Year + 2 * beta_Year_2 * 0
#     # = beta_Year
#     
#     effect_dist <- numeric(n_draws)
#     for (d in 1:n_draws) {
#       draw_row <- selected_draws[d, , drop = FALSE]
#       
#       # Get parameter value
#       get_param <- function(param_name) {
#         if (is.data.frame(draw_row)) {
#           val <- draw_row[[param_name]]
#           if (length(val) == 0) {
#             stop(paste("Parameter not found:", param_name))
#           }
#           return(val[1])
#         } else {
#           return(draw_row[[param_name]])
#         }
#       }
#       
#       beta_year <- get_param("beta_Year")
#       
#       # Marginal effect at mean Year = beta_Year
#       # (Since Year is centered, Year = 0 at the mean)
#       effect_dist[d] <- beta_year
#     }
#     
#     return(effect_dist)
#   }
#   
#   # ORIGINAL CODE FOR ALL OTHER VARIABLES
#   # Sample observations from the actual data
#   n_available <- nrow(data_df)
#   if (n_obs > n_available) n_obs <- n_available
#   
#   obs_indices <- sample(n_available, n_obs, replace = FALSE)
#   sampled_obs <- data_df[obs_indices, ]
#   
#   # Determine intervention magnitude (1 SD)
#   if (var_name %in% names(sds_df)) {
#     intervention_delta <- sds_df[[1, var_name]]
#   } else {
#     intervention_delta <- 1  # Default for derived variables
#   }
#   
#   # Select parameter draws
#   total_draws_available <- nrow(draws_df)
#   if (n_draws > total_draws_available) n_draws <- total_draws_available
#   
#   if (total_draws_available == 1) {
#     selected_draws <- draws_df
#     n_draws <- 1
#   } else {
#     sample_indices <- sample(total_draws_available, n_draws, replace = FALSE)
#     selected_draws <- draws_df[sample_indices, , drop = FALSE]
#   }
#   
#   # Calculate effects for each draw and each observation
#   effects_matrix <- matrix(NA, nrow = n_draws, ncol = n_obs)
#   
#   for (d in 1:n_draws) {
#     draw_row <- selected_draws[d, , drop = FALSE]
#     
#     for (i in 1:n_obs) {
#       obs_row <- sampled_obs[i, ]
#       
#       # Baseline prediction (using actual data values)
#       mu_goby_baseline <- calculate_mu_goby_full_dag(
#         draw_row, 
#         obs_row, 
#         intervention_var = NULL
#       )
#       
#       # Intervened prediction (increase variable by 1 SD)
#       mu_goby_intervene <- calculate_mu_goby_full_dag(
#         draw_row, 
#         obs_row, 
#         intervention_var = var_name,
#         intervention_value = obs_row[[var_name]] + intervention_delta
#       )
#       
#       effects_matrix[d, i] <- mu_goby_intervene - mu_goby_baseline
#     }
#   }
#   
#   # Average across observations for each draw
#   effect_dist <- rowMeans(effects_matrix, na.rm = TRUE)
#   
#   return(effect_dist)
# }
# 
# # ============================================================================
# # PART 1: FIT THE SCM MODEL ON SIMULATED DATA
# # ============================================================================
# 
# # Make data into a list
# list_df_sim <- lapply(as.list(sim_data), function(x) {
#   if (is.matrix(x) && ncol(x) == 1) {
#     return(as.vector(x))
#   } else {
#     return(x)
#   }
# })
# 
# list_df_sim$N <- length(list_df_sim$Goby)
# list_df_sim$J <- length(unique(list_df_sim$Zone))
# 
# # Run the SCM model
# mod.SEM.sim <- mod.SEM$sample(
#   data = list_df_sim,
#   seed = 123,
#   chains = 3,
#   iter_warmup = 6000,
#   iter_sampling = 1000,
#   parallel_chains = 3,
#   refresh = 100,
#   adapt_delta = 0.95,
#   output_dir = "stan_output",
#   save_warmup = FALSE
# )
# 
# saveRDS(mod.SEM.sim, file = "Output/Models/mod.SEM.sim.rds")
# mod.SEM.sim$diagnostic_summary()
# print(mod.SEM.sim$summary(), n = 130)
# 
# # ============================================================================
# # PART 2: EXTRACT DIRECT EFFECT COEFFICIENTS
# # ============================================================================
# 
# coefficient_summary.sim <- summarise_draws(
#   mod.SEM.sim,
#   "mean",
#   ~quantile(.x, probs = c(0.045, 0.955))
# )
# 
# # SCM betas
# beta.joint.sim <- coefficient_summary.sim %>%
#   dplyr::filter(str_starts(variable, "beta_")) %>%
#   dplyr::mutate(variable = str_remove(variable, "^.*?_")) %>%
#   # Filter out Rain and Wind - these are path coefficients, not direct effects on Goby
#   dplyr::filter(!variable %in% c("Rain", "Wind"))
# beta.joint.sim$model <- "SCM"
# beta.joint.sim$effect_type <- "Direct"
# 
# # GLMM betas (assuming you have model.brms.sim fitted)
# coefficient_summary.sim.brm <- summarise_draws(
#   model.brms.sim,
#   "mean",
#   ~quantile(.x, probs = c(0.045, 0.955))
# )
# beta.sim.brm <- coefficient_summary.sim.brm %>%
#   dplyr::filter(str_starts(variable, "b_")) %>%
#   dplyr::filter(variable != "b_Intercept") %>%
#   dplyr::mutate(variable = str_remove(variable, "^.*?_"))
# beta.sim.brm$model <- "GLMM"
# beta.sim.brm$effect_type <- "Direct"
# 
# # True direct effects from simulation - Use ACTUAL values
# betas_true_goby <- c(
#   # Rain and Wind have NO direct effect on Goby (only indirect through pathways)
#   # Show as zero to contrast with GLMM's incorrect estimates
#   Rain = 0,
#   Wind = 0,
#   
#   Substrate = true_params$beta_Substrate_vec[1],
#   Temp = true_params$beta_Temp_vec[1],
#   Temp_2 = true_params$beta_Temp_2,
#   DO = true_params$beta_DO_vec[1],
#   SAV = true_params$beta_SAV_vec[1],
#   SAV_2 = true_params$beta_SAV_2,
#   SC_count = true_params$beta_SC_count,
#   SB_count = true_params$beta_SB_count,
#   BreachDays = true_params$beta_BreachDays_vec[1],
#   BreachDays_2 = true_params$beta_BreachDays_2,
#   Micro = true_params$beta_Micro,
#   Year = true_params$beta_Year,
#   Year_2 = true_params$beta_Year_2,
#   Goby_lag = true_params$beta_Goby_lag
# )
# 
# # Path coefficients from submodels
# betas_true_paths <- c(
#   "DO_vec[2]" = true_params$beta_DO_vec[2],
#   "DO_vec[3]" = true_params$beta_DO_vec[3],
#   "DO_vec[4]" = true_params$beta_DO_vec[4],
#   "SAV_vec[2]" = true_params$beta_SAV_vec[2],
#   "SAV_vec[3]" = true_params$beta_SAV_vec[3],
#   "Temp_vec[2]" = true_params$beta_Temp_vec[2],
#   "Temp_vec[3]" = true_params$beta_Temp_vec[3],
#   "BreachDays_vec[2]" = true_params$beta_BreachDays_vec[2],
#   "Wind_vec[1]" = true_params$beta_Wind_vec[1],
#   "Wind_vec[2]" = true_params$beta_Wind_vec[2],
#   "Substrate_vec[2]" = true_params$beta_Substrate_vec[2]
# )
# 
# # Combine all true direct effects
# betas_true_direct <- c(betas_true_goby, betas_true_paths)
# 
# df.true.direct <- tibble(
#   variable = names(betas_true_direct),
#   mean = betas_true_direct,
#   `4.5%` = NA,
#   `95.5%` = NA,
#   model = "True",
#   effect_type = "Direct"
# )
# 
# # ============================================================================
# # PART 3: CALCULATE TRUE TOTAL EFFECTS VIA SIMULATION
# # ============================================================================
# 
# cat("\n=== Calculating TRUE total effects via simulation (with Year fix) ===\n")
# 
# # Calculate standard deviations
# sim_data$Goby_Density <- sim_data$Goby / exp(sim_data$Area)
# numeric_cols <- names(sim_data)[sapply(sim_data, is.numeric)]
# sds <- data.frame(
#   purrr::map_dfc(sim_data[numeric_cols], ~ sd(.x, na.rm = TRUE))
# )
# colnames(sds) <- numeric_cols
# 
# # Define exogenous/early variables to test PLUS endogenous mediators
# total_effect_vars <- c("Rain", "Wind", "Substrate", "Micro", "Year",
#                        "BreachDays", "Temp", "DO", "SAV")  # Added endogenous variables
# 
# # CREATE DATAFRAME WITH TRUE PARAMETER VALUES FROM SIMULATION
# true_params_df <- data.frame(
#   beta_Rain = true_params$beta_Rain,
#   beta_Wind_vec.1. = true_params$beta_Wind_vec[1],
#   beta_Wind_vec.2. = true_params$beta_Wind_vec[2],
#   beta_Substrate_vec.1. = true_params$beta_Substrate_vec[1],
#   beta_Substrate_vec.2. = true_params$beta_Substrate_vec[2],
#   beta_Temp_vec.1. = true_params$beta_Temp_vec[1],
#   beta_Temp_vec.2. = true_params$beta_Temp_vec[2],
#   beta_Temp_vec.3. = true_params$beta_Temp_vec[3],
#   beta_Temp_2 = true_params$beta_Temp_2,
#   beta_DO_vec.1. = true_params$beta_DO_vec[1],
#   beta_DO_vec.2. = true_params$beta_DO_vec[2],
#   beta_DO_vec.3. = true_params$beta_DO_vec[3],
#   beta_DO_vec.4. = true_params$beta_DO_vec[4],
#   beta_SAV_vec.1. = true_params$beta_SAV_vec[1],
#   beta_SAV_vec.2. = true_params$beta_SAV_vec[2],
#   beta_SAV_vec.3. = true_params$beta_SAV_vec[3],
#   beta_SAV_2 = true_params$beta_SAV_2,
#   beta_SC_count = true_params$beta_SC_count,
#   beta_SB_count = true_params$beta_SB_count,
#   beta_BreachDays_vec.1. = true_params$beta_BreachDays_vec[1],
#   beta_BreachDays_vec.2. = true_params$beta_BreachDays_vec[2],
#   beta_BreachDays_2 = true_params$beta_BreachDays_2,
#   beta_Micro = true_params$beta_Micro,
#   beta_Year = true_params$beta_Year,
#   beta_Year_2 = true_params$beta_Year_2,
#   beta_Goby_lag = true_params$beta_Goby_lag,
#   a_BreachDays = 0,
#   a_Temp = 0,
#   a_DO = 0,
#   a_SAV = 0,
#   a_SC = -0.5,
#   a_SB = -0.2,
#   a_Goby.1. = true_params$a_Goby_zone[1],
#   a_Goby.2. = true_params$a_Goby_zone[2],
#   a_Goby.3. = true_params$a_Goby_zone[3],
#   U_bio.1. = true_params$U_bio[1],
#   U_bio.2. = true_params$U_bio[2],
#   U_bio.3. = true_params$U_bio[3],
#   U_phys.1. = true_params$U_phys[1],
#   U_phys.2. = true_params$U_phys[2],
#   U_phys.3. = true_params$U_phys[3],
#   k_U_bio = true_params$k_U_bio,
#   k_U_phys = true_params$k_U_phys
# )
# 
# # Fix column names
# names(true_params_df) <- gsub("\\.", "[", names(true_params_df))
# names(true_params_df) <- gsub("\\[([0-9]+)\\[$", "[\\1]", names(true_params_df))
# 
# # Verify key parameters
# cat("\n=== VERIFYING TRUE PARAMETERS ===\n")
# cat(sprintf("beta_Micro: %.4f (should be 0.2)\n", true_params_df$beta_Micro[1]))
# cat(sprintf("beta_Year: %.4f (should be 0.2)\n", true_params_df$beta_Year[1]))
# cat(sprintf("beta_Year_2: %.4f (should be -0.2)\n", true_params_df$beta_Year_2[1]))
# cat(sprintf("beta_Rain: %.4f (should be 0.5)\n", true_params_df$beta_Rain[1]))
# cat("===================================\n\n")
# 
# # Calculate TRUE total effects using the updated function
# true_total_effects_list <- purrr::map(
#   total_effect_vars,
#   function(var_name) {
#     message(paste("  Calculating TRUE total effect for:", var_name))
#     
#     tryCatch({
#       effect_value <- calculate_total_effect(
#         var_name, 
#         draws_df = true_params_df,
#         data_df = sim_data, 
#         sds_df = sds, 
#         n_draws = 1,
#         n_obs = 100
#       )
#       
#       tibble(
#         variable = paste0(var_name, " (Total)"),
#         mean = effect_value[1],
#         `4.5%` = NA,
#         `95.5%` = NA,
#         model = "True",
#         effect_type = "Total"
#       )
#     }, error = function(e) {
#       message(paste("  ERROR:", var_name, "-", e$message))
#       return(NULL)
#     })
#   }
# )
# 
# true_total_effects <- bind_rows(true_total_effects_list[!sapply(true_total_effects_list, is.null)])
# 
# # ============================================================================
# # PART 4: CALCULATE ESTIMATED TOTAL EFFECTS FROM SCM
# # ============================================================================
# 
# cat("\n=== Calculating SCM Total Effects (with Year fix) ===\n")
# 
# draws_scm <- mod.SEM.sim$draws(format = "df")
# 
# scm_total_effects_list <- purrr::map(
#   total_effect_vars,
#   function(var_name) {
#     tryCatch({
#       effect_dist <- calculate_total_effect(
#         var_name, 
#         draws_scm, 
#         sim_data, 
#         sds, 
#         n_draws = 100,
#         n_obs = 100
#       )
#       
#       tibble(
#         variable = paste0(var_name, " (Total)"),
#         mean = mean(effect_dist, na.rm = TRUE),
#         `5.5%` = quantile(effect_dist, 0.055, na.rm = TRUE),
#         `94.5%` = quantile(effect_dist, 0.945, na.rm = TRUE),
#         model = "SCM",
#         effect_type = "Total"
#       )
#     }, error = function(e) {
#       message(paste("  ERROR:", var_name, "-", e$message))
#       return(NULL)
#     })
#   }
# )
# 
# scm_total_effects <- bind_rows(scm_total_effects_list[!sapply(scm_total_effects_list, is.null)])
# 
# # ============================================================================
# # PART 4B: CALCULATE GLMM "TOTAL EFFECTS" (ACTUALLY CONFOUNDED)
# # ============================================================================
# 
# cat("\n=== Calculating GLMM 'Total Effects' (Confounded) ===\n")
# 
# # For GLMM, the "total effect" is just the direct coefficient because
# # GLMM doesn't separate direct from mediated effects - it's all confounded!
# # This will show how GLMM fails to properly estimate causal effects
# 
# draws_glmm <- as_draws_df(model.brms.sim)
# 
# # Extract GLMM coefficients for our variables of interest
# glmm_total_effects_list <- purrr::map(
#   total_effect_vars,
#   function(var_name) {
#     # Get the coefficient name in brms format
#     coef_name <- paste0("b_", var_name)
#     
#     # Check if this coefficient exists in the GLMM
#     if (coef_name %in% names(draws_glmm)) {
#       coef_draws <- draws_glmm[[coef_name]]
#       
#       # SPECIAL HANDLING FOR YEAR
#       if (var_name == "Year") {
#         message("  GLMM Year: Reporting beta_Year only (marginal effect at mean)")
#         # For Year, just report beta_Year (not including Year_2)
#         # This is the marginal effect at mean Year, comparable to SCM
#       }
#       
#       tibble(
#         variable = paste0(var_name, " (Total)"),
#         mean = mean(coef_draws, na.rm = TRUE),
#         `4.5%` = quantile(coef_draws, 0.045, na.rm = TRUE),
#         `95.5%` = quantile(coef_draws, 0.955, na.rm = TRUE),
#         model = "GLMM",
#         effect_type = "Total"
#       )
#     } else {
#       # Variable not in GLMM (e.g., endogenous variables)
#       return(NULL)
#     }
#   }
# )
# 
# glmm_total_effects <- bind_rows(glmm_total_effects_list[!sapply(glmm_total_effects_list, is.null)])
# 
# # ============================================================================
# # PART 5: COMBINE AND PLOT
# # ============================================================================
# 
# # Combine direct effects
# df.direct <- bind_rows(
#   df.true.direct %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`),
#   beta.joint.sim %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`),
#   beta.sim.brm %>% dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)
# )
# 
# # Combine total effects (now including GLMM)
# # First, standardize column names across all dataframes
# true_total_effects_renamed <- true_total_effects %>% 
#   dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)
# 
# scm_total_effects_renamed <- scm_total_effects %>% 
#   dplyr::rename(lower_4.5 = `5.5%`, upper_95.5 = `94.5%`)
# 
# glmm_total_effects_renamed <- glmm_total_effects %>% 
#   dplyr::rename(lower_4.5 = `4.5%`, upper_95.5 = `95.5%`)
# 
# df.total <- bind_rows(
#   true_total_effects_renamed,
#   scm_total_effects_renamed,
#   glmm_total_effects_renamed
# )
# 
# # Relabel direct effects
# df.direct_relabeled <- df.direct %>%
#   mutate(
#     original_variable = variable,
#     variable = if_else(str_detect(variable, "_vec\\[1\\]"), str_remove(variable, "_vec\\[1\\]"), variable),
#     variable = str_replace(variable, "_2$", "²"),
#     variable = case_when(
#       original_variable == "DO_vec[2]" ~ "DO → SB",
#       original_variable == "DO_vec[3]" ~ "DO → SC",
#       original_variable == "DO_vec[4]" ~ "DO → SAV",
#       original_variable == "SAV_vec[2]" ~ "SAV → SB",
#       original_variable == "SAV_vec[3]" ~ "SAV → SC",
#       original_variable == "Temp_vec[2]" ~ "Temp → DO",
#       original_variable == "Temp_vec[3]" ~ "Temp → SAV",
#       original_variable == "BreachDays_vec[2]" ~ "BreachDays → Temp",
#       original_variable == "Wind_vec[2]" ~ "Wind → Temp",
#       original_variable == "Substrate_vec[2]" ~ "Substrate → SC",
#       TRUE ~ variable
#     ),
#     variable = str_replace(variable, "^Goby_lag$", "Goby lag"),
#     variable = str_replace(variable, "^SB_count$", "SB"),
#     variable = str_replace(variable, "^SC_count$", "SC")
#   ) %>%
#   dplyr::select(-original_variable)
# 
# # Create plots
# df_direct_plot <- df.direct_relabeled %>% mutate(variable = factor(variable))
# axis_labels_direct <- levels(fct_rev(df_direct_plot$variable))
# label_colors_direct <- ifelse(str_detect(axis_labels_direct, "Wind|Rain"), "blue", "black")
# 
# p.direct <- df_direct_plot %>%
#   ggplot(aes(x = fct_rev(variable), y = mean, color = model)) +
#   geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), size = 0.7, position = position_dodge(width = 0.7)) +
#   coord_flip() +
#   geom_hline(yintercept = 0, linetype = 2) +
#   scale_color_manual(values = c("True" = "black", "SCM" = "red3", "GLMM" = "blue")) +
#   theme_minimal(base_size = 18) +
#   theme(
#     legend.position = c(0.80, 0.96),
#     legend.title = element_blank(),
#     legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = NA),
#     axis.text.y = element_text(colour = label_colors_direct),
#     plot.title = element_text(face = "bold", size = 20)
#   ) +
#   guides(color = guide_legend(reverse = TRUE)) +
#   labs(title = "A. Direct Effects on Goby", y = "Beta estimate", x = "")
# 
# # Total effects plot
# df_total_plot <- df.total %>% mutate(variable = factor(variable))
# axis_labels_total <- levels(fct_rev(df_total_plot$variable))
# label_colors_total <- ifelse(str_detect(axis_labels_total, "Wind|Rain"), "blue", "black")
# 
# p.total <- df_total_plot %>%
#   ggplot(aes(x = fct_rev(variable), y = mean, color = model)) +
#   geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), size = 0.7, position = position_dodge(width = 0.7)) +
#   coord_flip() +
#   geom_hline(yintercept = 0, linetype = 2) +
#   scale_color_manual(values = c("True" = "black", "SCM" = "red3", "GLMM" = "blue"), 
#                      breaks = c("True", "SCM", "GLMM")) +
#   theme_minimal(base_size = 18) +
#   theme(
#     legend.position = c(0.85, 0.96),
#     legend.title = element_blank(),
#     legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = NA),
#     axis.text.y = element_text(colour = label_colors_total),
#     plot.title = element_text(face = "bold", size = 20)
#   ) +
#   guides(color = guide_legend(reverse = TRUE)) +
#   labs(title = "B. Total Causal Effects on Goby", 
#        subtitle = "Direct + all mediated pathways (Year at mean; GLMM confounds effects)", 
#        y = "Total effect (log scale)", x = "")
# 
# p.combined <- p.direct / p.total + plot_layout(heights = c(2, 1))
# print(p.combined)
# 
# # Save plots
# ggsave("Output/Plots/sim.forest.combined.png", plot = p.combined, width = 25, height = 35, units = "cm", dpi = 300)
# ggsave("Output/Plots/sim.forest.direct.png", plot = p.direct, width = 20, height = 30, units = "cm")
# ggsave("Output/Plots/sim.forest.total.png", plot = p.total, width = 20, height = 15, units = "cm")
# 
# # Calculate coverage
# df.all <- bind_rows(df.direct_relabeled, df.total)
# coverage_summary <- df.all %>%
#   dplyr::filter(model != "True") %>%
#   left_join(df.all %>% dplyr::filter(model == "True") %>% dplyr::select(variable, effect_type, true_mean = mean),
#             by = c("variable", "effect_type")) %>%
#   mutate(covered = (true_mean >= lower_4.5) & (true_mean <= upper_95.5),
#          abs_error = abs(mean - true_mean),
#          rel_error = abs_error / abs(true_mean))
# 
# coverage_by_model <- coverage_summary %>%
#   group_by(model, effect_type) %>%
#   summarise(n_effects = n(), 
#             coverage_rate = mean(covered, na.rm = TRUE),
#             mean_abs_error = mean(abs_error, na.rm = TRUE),
#             mean_rel_error = mean(rel_error, na.rm = TRUE),
#             .groups = "drop")
# 
# cat("\n=== PARAMETER RECOVERY SUMMARY ===\n")
# print(coverage_by_model, n = 20)
# 
# cat("\n=== TOTAL EFFECT COMPARISON ===\n")
# total_comparison <- df.total %>%
#   dplyr::select(variable, model, mean, lower_4.5, upper_95.5) %>%
#   arrange(variable, model)
# print(total_comparison, n = 30)
# 
# cat("\n=== YEAR EFFECT COMPARISON ===\n")
# cat("All models now report the marginal effect of Year at mean Year\n")
# cat("This is beta_Year (the linear coefficient) when Year is centered\n\n")
# 
# year_comparison <- df.total %>%
#   filter(str_detect(variable, "Year")) %>%
#   select(variable, model, mean, lower_4.5, upper_95.5) %>%
#   arrange(model)
# 
# print(year_comparison)
# 
# cat("\n=== EXPLANATION ===\n")
# cat("Year has a quadratic term (Year²), so its effect is non-linear.\n")
# cat("We report the marginal effect at the mean Year (centered at 0).\n")
# cat("At Year = 0: d(response)/d(Year) = beta_Year + 2*beta_Year_2*0 = beta_Year\n")
# cat("\nThis gives a single, comparable number across all three models.\n")
# cat("For the full non-linear relationship, see a marginal effects curve.\n")
# cat("===========================================\n\n")
# 
# cat("\n=== ANALYSIS COMPLETE ===\n")
# cat("Plots saved to Output/Plots/\n\n")
# cat("Key insight: Total effects are more accurate and interpretable than individual paths!\n")
# cat("SCM correctly recovers TOTAL causal effects, which is what decision-makers need.\n")

