# SIMULATION FOR ZONE-LEVEL LATENT ENVIRONMENT MODEL
# Includes measurement error in Temp and DO
# SC and SB as measured (no submodels)
# 2026-01-12

library(MASS)

# --- 1. Set Constants & True Parameter Values ---
set.seed(42)
N <- 314
J <- 3

## Hyperpriors & Goby-specific parameters
a_bar_Goby <- 2.6
sigma_Goby <- 0.5
phi <- 1.0

## Pathway coefficients (simplified - no hierarchical except Wind)
beta_Rain_Breach <- 0.5
beta_Breach_Temp <- 0.4
beta_Temp_DO <- -0.4  # Physics: warm water holds less O2
beta_Wind_Temp <- -0.3
beta_Wind_DO <- 0.4

# Wind effects (the only hierarchical structure)
b_bar_Wind <- c(0.4, -0.3)  # [1]=DO, [2]=Temp
sigma_Wind <- c(0.1, 0.1)
rho_Wind <- 0.3

make_sigma <- function(k, sd_vec, rho) {
  rho_mat <- matrix(rho, nrow = k, ncol = k)
  diag(rho_mat) <- 1
  sigma_mat <- diag(sd_vec[1:k]) %*% rho_mat %*% diag(sd_vec[1:k])
  return(sigma_mat)
}

Sigma_Wind <- make_sigma(2, sigma_Wind, rho_Wind)
beta_Wind_vec <- mvrnorm(1, b_bar_Wind, Sigma_Wind)

## Direct effects on Goby
beta_BreachDays_Goby <- 0.3
beta_Temp_Goby <- 0.3
beta_DO_Goby <- 0.6
beta_SAV_Goby <- 0.4
beta_SC_count <- -0.5
beta_SB_count <- 0.5
beta_Substrate <- -0.3
beta_Micro <- 0.2
beta_Year <- 0.2
beta_Year_2 <- -0.2
beta_Goby_lag <- 0.1

## Quadratic terms
beta_BreachDays_2 <- -0.1
beta_Temp_2 <- -0.4
beta_SAV_2 <- -0.3

## Submodel intercepts
a_BreachDays <- 0
a_Temp <- 0
a_DO <- 0

## Error term
tau <- 0.5

## Latent variable parameters (only U_phys)
k_U_phys <- 0.4

## Wind/Rain correlation parameters
a_windrain <- c(0, 0)
rho_windrain <- 0.4
sigma_wind_weather <- 1.0
sigma_rain <- 1.0

## *** MEASUREMENT ERROR PARAMETERS ***
sigma_Temp_measure <- 0.5  # Substantial measurement error in Temp
sigma_DO_measure <- 0.4    # Substantial measurement error in DO
sigma_Temp_within <- 0.3   # Within-zone true variation
sigma_DO_within <- 0.3     # Within-zone true variation

# --- 2. Simulate Zone-Level Parameters ---
a_Goby_zone <- rnorm(J, mean = a_bar_Goby, sd = sigma_Goby)

## *** Only U_phys (physical confounding) ***
U_phys <- rnorm(J, mean = 0, sd = 1)

# --- 3. Simulate Observation-Level Data ---
sim_data <- data.frame(Zone = sample(1:J, size = N, replace = TRUE))

sim_data$Area <- rnorm(N, mean = 2.4, sd = 0.5)
sim_data$Substrate <- rbinom(N, size = 1, prob = 0.5)

## Generate correlated Wind and Rain
Sigma_windrain <- matrix(c(
  sigma_wind_weather^2, rho_windrain * sigma_wind_weather * sigma_rain,
  rho_windrain * sigma_wind_weather * sigma_rain, sigma_rain^2
), nrow = 2)

windrain_data <- mvrnorm(N, mu = a_windrain, Sigma = Sigma_windrain)
sim_data$Wind <- windrain_data[, 1]
sim_data$Rain <- windrain_data[, 2]

sim_data$Micro <- rnorm(N)
sim_data$Year <- rnorm(N)
sim_data$Goby_lag <- rnorm(N)

# --- 4. Generative Cascade with ZONE-LEVEL TRUE ENVIRONMENT ---

## BreachDays (affected by Rain + U_phys)
mu_breach <- a_BreachDays + 
  beta_Rain_Breach * sim_data$Rain + 
  k_U_phys * U_phys[sim_data$Zone]
sim_data$BreachDays <- rnorm(N, mean = mu_breach, sd = tau)

## *** ZONE-LEVEL TRUE TEMPERATURE ***
Temp_true_zone <- numeric(J)
for (j in 1:J) {
  obs_in_zone <- which(sim_data$Zone == j)
  avg_breach <- mean(sim_data$BreachDays[obs_in_zone])
  avg_wind <- mean(sim_data$Wind[obs_in_zone])
  
  Temp_true_zone[j] <- a_Temp + 
    beta_Breach_Temp * avg_breach + 
    beta_Wind_vec[2] * avg_wind + 
    k_U_phys * U_phys[j]
}

## *** ZONE-LEVEL TRUE DO ***
DO_true_zone <- numeric(J)
for (j in 1:J) {
  obs_in_zone <- which(sim_data$Zone == j)
  avg_wind <- mean(sim_data$Wind[obs_in_zone])
  
  DO_true_zone[j] <- a_DO + 
    beta_Temp_DO * Temp_true_zone[j] + 
    beta_Wind_vec[1] * avg_wind + 
    k_U_phys * U_phys[j]
}

## Assign TRUE zone-level environment to observations
sim_data$Temp_true <- Temp_true_zone[sim_data$Zone]
sim_data$DO_true <- DO_true_zone[sim_data$Zone]

## *** MEASURED VALUES WITH ERROR ***
# Each measurement deviates from zone mean due to:
# 1. True within-zone variation
# 2. Measurement error from limited sampling
total_Temp_sd <- sqrt(sigma_Temp_within^2 + sigma_Temp_measure^2)
total_DO_sd <- sqrt(sigma_DO_within^2 + sigma_DO_measure^2)

sim_data$Temp_measured <- rnorm(N, mean = sim_data$Temp_true, sd = total_Temp_sd)
sim_data$DO_measured <- rnorm(N, mean = sim_data$DO_true, sd = total_DO_sd)

## SAV generated but treated as measured (no submodel in Stan)
mu_sav_true <- 0 + 0.5 * sim_data$DO_true + 0.3 * sim_data$Temp_true + rnorm(N, 0, 0.2)
sim_data$SAV <- rnorm(N, mean = mu_sav_true, sd = tau)

## *** SC and SB as measured (no submodels) ***
# Generate them as if affected by environment, but model treats as measured
logit_sc <- -0.5 + 
  -0.3 * sim_data$Substrate + 
  0.3 * sim_data$DO_true + 
  0.3 * sim_data$SAV
sim_data$SC_count <- rbinom(N, size = 1, prob = plogis(logit_sc))

logit_sb <- -0.2 + 
  0.3 * sim_data$DO_true + 
  0.3 * sim_data$SAV
sim_data$SB_count <- rbinom(N, size = 1, prob = plogis(logit_sb))

# --- 5. Simulate the Final Goby Outcome (responds to TRUE environment) ---
log_mu_goby <- (a_Goby_zone[sim_data$Zone] +
                  sim_data$Area +
                  beta_SAV_Goby * sim_data$SAV +
                  beta_SAV_2 * sim_data$SAV^2 +
                  beta_DO_Goby * sim_data$DO_true +  # TRUE DO
                  beta_BreachDays_Goby * sim_data$BreachDays +
                  beta_BreachDays_2 * sim_data$BreachDays^2 +
                  beta_Temp_Goby * sim_data$Temp_true +  # TRUE Temp
                  beta_Temp_2 * sim_data$Temp_true^2 +
                  beta_SC_count * sim_data$SC_count +
                  beta_SB_count * sim_data$SB_count +
                  beta_Substrate * sim_data$Substrate +
                  beta_Micro * sim_data$Micro +
                  beta_Year * sim_data$Year + 
                  beta_Year_2 * sim_data$Year^2 +
                  beta_Goby_lag * sim_data$Goby_lag)

mu_goby <- exp(log_mu_goby)
sim_data$Goby <- rnbinom(N, mu = mu_goby, size = phi)

# --- 6. Calculate Goby Density ---
sim_data$Goby_Density <- sim_data$Goby / exp(sim_data$Area)

# --- 7. Finalize ---
sim_data$Year_2 <- sim_data$Year^2
sim_data$Temp_2 <- sim_data$Temp_measured^2  # For data preparation
sim_data$BreachDays_2 <- sim_data$BreachDays^2
sim_data$SAV_2 <- sim_data$SAV^2

# Reorder columns - NOTE: We keep both measured and true for comparison
sim_data <- sim_data[, c("Goby", "Goby_Density", "Zone", "Area", "Rain", "Wind", 
                         "BreachDays", 
                         "Temp_measured", "Temp_true",  # Both versions
                         "DO_measured", "DO_true",      # Both versions
                         "SAV", "SC_count", "SB_count", 
                         "Substrate", "Micro", "Year", "Goby_lag", "Year_2", 
                         "Temp_2", "BreachDays_2", "SAV_2")]

# Check the summaries
cat("--- Summary of Simulated Goby Counts ---\n")
print(summary(sim_data$Goby))

cat("\n--- Measurement Error Check ---\n")
cat("True Temp (zone means):", round(Temp_true_zone, 2), "\n")
cat("Measured Temp range:", round(range(sim_data$Temp_measured), 2), "\n")
cat("Correlation (measured vs true Temp):", round(cor(sim_data$Temp_measured, sim_data$Temp_true), 3), "\n")
cat("Correlation (measured vs true DO):", round(cor(sim_data$DO_measured, sim_data$DO_true), 3), "\n")

# --- 8. Store true parameter values ---
true_params <- list(
  # Zone-level parameters
  a_Goby_zone = a_Goby_zone,
  U_phys = U_phys,
  Temp_true_zone = Temp_true_zone,
  DO_true_zone = DO_true_zone,
  
  # Latent variable loading
  k_U_phys = k_U_phys,
  
  # Measurement error
  sigma_Temp_measure = sigma_Temp_measure,
  sigma_DO_measure = sigma_DO_measure,
  sigma_Temp_within = sigma_Temp_within,
  sigma_DO_within = sigma_DO_within,
  
  # Wind hierarchical
  beta_Wind_vec = beta_Wind_vec,
  
  # Pathway coefficients
  beta_Rain_Breach = beta_Rain_Breach,
  beta_Breach_Temp = beta_Breach_Temp,
  beta_Temp_DO = beta_Temp_DO,
  beta_Wind_Temp = beta_Wind_Temp,
  beta_Wind_DO = beta_Wind_DO,
  
  # Direct effects on Goby
  beta_BreachDays_Goby = beta_BreachDays_Goby,
  beta_Temp_Goby = beta_Temp_Goby,
  beta_DO_Goby = beta_DO_Goby,
  beta_SAV_Goby = beta_SAV_Goby,
  beta_SC_count = beta_SC_count,
  beta_SB_count = beta_SB_count,
  beta_Substrate = beta_Substrate,
  beta_Micro = beta_Micro,
  beta_Year = beta_Year,
  beta_Year_2 = beta_Year_2,
  beta_Goby_lag = beta_Goby_lag,
  
  # Quadratic terms
  beta_BreachDays_2 = beta_BreachDays_2,
  beta_Temp_2 = beta_Temp_2,
  beta_SAV_2 = beta_SAV_2,
  
  # Dispersion and error
  phi = phi,
  tau = tau,
  
  # Wind/Rain correlation
  rho_windrain = rho_windrain,
  sigma_wind = sigma_wind_weather,
  sigma_rain = sigma_rain
)

cat("\n=== TRUE PARAMETER VALUES ===\n")
cat("\n--- Zone-Level True Environment ---\n")
cat("Temp_true_zone:", round(Temp_true_zone, 3), "\n")
cat("DO_true_zone:", round(DO_true_zone, 3), "\n")

cat("\n--- Measurement Error ---\n")
cat("sigma_Temp_measure:", sigma_Temp_measure, "\n")
cat("sigma_DO_measure:", sigma_DO_measure, "\n")
cat("sigma_Temp_within:", sigma_Temp_within, "\n")
cat("sigma_DO_within:", sigma_DO_within, "\n")

cat("\n--- Key Pathways ---\n")
cat("beta_Rain_Breach:", beta_Rain_Breach, "\n")
cat("beta_Breach_Temp:", beta_Breach_Temp, "\n")
cat("beta_Temp_DO:", beta_Temp_DO, "\n")
cat("beta_Wind_vec [DO, Temp]:", round(beta_Wind_vec, 3), "\n")

cat("\n--- Direct Effects on Goby ---\n")
cat("beta_Temp_Goby:", beta_Temp_Goby, "\n")
cat("beta_DO_Goby:", beta_DO_Goby, "\n")
cat("beta_BreachDays_Goby:", beta_BreachDays_Goby, "\n")

# ======================================================
# 9. CALCULATE TRUE TOTAL CAUSAL EFFECTS
# ======================================================

cat("\n==============================================\n")
cat(" TRUE TOTAL CAUSAL EFFECTS\n")
cat("==============================================\n")

# Calculate at mean values
mean_Temp <- mean(sim_data$Temp_true)
mean_Breach <- mean(sim_data$BreachDays)
mean_SAV <- mean(sim_data$SAV)

# Marginal slopes
slope_Temp_to_Goby <- beta_Temp_Goby + (2 * beta_Temp_2 * mean_Temp)
slope_Breach_to_Goby <- beta_BreachDays_Goby + (2 * beta_BreachDays_2 * mean_Breach)
slope_SAV_to_Goby <- beta_SAV_Goby + (2 * beta_SAV_2 * mean_SAV)

# RAIN TOTAL EFFECT
rain_total <- (
  beta_Rain_Breach * slope_Breach_to_Goby +
    beta_Rain_Breach * beta_Breach_Temp * slope_Temp_to_Goby +
    beta_Rain_Breach * beta_Breach_Temp * beta_Temp_DO * beta_DO_Goby
)

# WIND TOTAL EFFECT  
wind_total <- (
  beta_Wind_vec[2] * slope_Temp_to_Goby +
    beta_Wind_vec[1] * beta_DO_Goby +
    beta_Wind_vec[2] * beta_Temp_DO * beta_DO_Goby
)

# BREACHDAYS TOTAL EFFECT
breach_total <- (
  slope_Breach_to_Goby +
    beta_Breach_Temp * slope_Temp_to_Goby +
    beta_Breach_Temp * beta_Temp_DO * beta_DO_Goby
)

# TEMP TOTAL EFFECT
temp_total <- (
  slope_Temp_to_Goby +
    beta_Temp_DO * beta_DO_Goby
)

# DO TOTAL EFFECT
do_total <- beta_DO_Goby

# SAV TOTAL EFFECT
sav_total <- slope_SAV_to_Goby

# Others
substrate_total <- beta_Substrate
micro_total <- beta_Micro
year_total <- beta_Year

true_total_effects <- data.frame(
  Variable = c("Rain", "Wind", "BreachDays", "Temp", "DO", "SAV", 
               "Substrate", "Micro", "Year"),
  Total_Effect = c(rain_total, wind_total, breach_total, temp_total, 
                   do_total, sav_total, substrate_total, micro_total, year_total),
  Note = c("Rain → Breach → Temp → DO → Goby",
           "Wind → Temp/DO → Goby",
           "Direct + Temp + DO pathways",
           "Direct + DO pathway",
           "Direct only",
           "Direct only (no submodel)",
           "Direct only",
           "Direct only",
           "At mean Year")
)

print(true_total_effects, digits = 4)

# Save everything
save(sim_data, true_params, true_total_effects, 
     file = "simulated_data_zone_latent.RData")

cat("\n=== DATA SAVED ===\n")
cat("File: simulated_data_zone_latent.RData\n")
cat("Contains: sim_data, true_params, true_total_effects\n")
cat("\nThis matches the zone_latent_environment_SEM.stan model\n")

cat("\n=== MODEL SPECIFICATIONS ===\n")
cat("Submodels: BreachDays, Temp, DO (zone-level latent)\n")
cat("Measured with error: Temp, DO\n")
cat("Measured directly: SAV, SC, SB\n")
cat("Latent variables: U_phys only\n")
cat("Hierarchical: Wind effects (DO, Temp)\n")
cat("Parameters: ~44\n")
cat("Sample size: N =", N, "\n")
cat("Zones: J =", J, "\n")
cat("\nKEY FEATURE: Temp_measured and DO_measured deviate from\n")
cat("zone-level truth (Temp_true_zone, DO_true_zone) due to\n")
cat("measurement error. Gobies respond to TRUE environment.\n")

###----
###old below


# # CORRECTED SIMULATION CODE WITH LATENT VARIABLES
# # 2025-12-28
# library(MASS)
# 
# # --- 1. Set Constants & True Parameter Values ---
# set.seed(42)
# N <- 314
# J <- 3
# 
# ## Hyperpriors & Goby-specific parameters
# a_bar_Goby <- 2.6
# sigma_Goby <- 0.5
# phi <- 1.0
# 
# ## --- Hierarchical Coefficient Hyperpriors ---
# b_bar_DO <- c(0.6, 0.6, 0.6, 0.6)
# b_bar_SAV <- c(0.4, 0.4, 0.4)
# b_bar_Temp <- c(0.4, 0.4, 0.4)
# b_bar_BreachDays <- c(0.3, 0.3)
# b_bar_Wind <- c(-0.5, -0.5)
# b_bar_Substrate <- c(-0.3, -0.3)
# 
# sigma_vec <- c(0.1, 0.1, 0.1, 0.1)
# rho_val <- 0.3
# 
# make_sigma <- function(k, sd_vec, rho) {
#   rho_mat <- matrix(rho, nrow = k, ncol = k)
#   diag(rho_mat) <- 1
#   sigma_mat <- diag(sd_vec[1:k]) %*% rho_mat %*% diag(sd_vec[1:k])
#   return(sigma_mat)
# }
# 
# Sigma_DO <- make_sigma(4, sigma_vec, rho_val)
# Sigma_SAV <- make_sigma(3, sigma_vec, rho_val)
# Sigma_Temp <- make_sigma(3, sigma_vec, rho_val)
# Sigma_BreachDays <- make_sigma(2, sigma_vec, rho_val)
# Sigma_Wind <- make_sigma(2, sigma_vec, rho_val)
# Sigma_Substrate <- make_sigma(2, sigma_vec, rho_val)
# 
# ## Draw coefficients
# beta_DO_vec <- mvrnorm(1, b_bar_DO, Sigma_DO)
# beta_SAV_vec <- mvrnorm(1, b_bar_SAV, Sigma_SAV)
# beta_Temp_vec <- mvrnorm(1, b_bar_Temp, Sigma_Temp)
# beta_BreachDays_vec <- mvrnorm(1, b_bar_BreachDays, Sigma_BreachDays)
# beta_Wind_vec <- mvrnorm(1, b_bar_Wind, Sigma_Wind)
# beta_Substrate_vec <- mvrnorm(1, b_bar_Substrate, Sigma_Substrate)
# 
# ## Other fixed coefficients
# beta_SC_count <- -0.5
# beta_SB_count <- 0.5
# beta_Micro <- 0.2
# beta_Year <- 0.2
# beta_Year_2 <- -0.2
# beta_Goby_lag <- 0.1
# beta_Temp_2 <- -0.4
# beta_BreachDays_2 <- -0.1
# beta_SAV_2 <- -0.3
# 
# ## Submodel intercepts and error terms
# a_BreachDays <- 0
# a_Temp <- 0
# a_DO <- 0
# a_SAV <- 0
# a_SC <- -0.5
# a_SB <- -0.2
# tau <- 0.5
# 
# ## *** NEW: Latent variable parameters ***
# k_U_bio <- 0.5   # Loading for biological latent variable
# k_U_phys <- 0.4  # Loading for physical latent variable
# 
# ## *** NEW: Wind/Rain correlation parameters ***
# a_windrain <- c(0, 0)  # Means for Wind and Rain
# rho_windrain <- 0.4    # Correlation between Wind and Rain
# sigma_wind <- 1.0
# sigma_rain <- 1.0
# 
# # --- 2. Simulate Zone-Level Parameters ---
# a_Goby_zone <- rnorm(J, mean = a_bar_Goby, sd = sigma_Goby)
# 
# ## *** NEW: Simulate zone-level latent variables ***
# U_bio <- rnorm(J, mean = 0, sd = 1)   # Standard normal for each zone
# U_phys <- rnorm(J, mean = 0, sd = 1)  # Standard normal for each zone
# 
# # --- 3. Simulate Observation-Level Data ---
# sim_data <- data.frame(Zone = sample(1:J, size = N, replace = TRUE))
# 
# sim_data$Area <- rnorm(N, mean = 2.4, sd = 0.5)
# sim_data$Substrate <- rbinom(N, size = 1, prob = 0.5)
# 
# ## *** NEW: Generate correlated Wind and Rain ***
# Sigma_windrain <- matrix(c(
#   sigma_wind^2, rho_windrain * sigma_wind * sigma_rain,
#   rho_windrain * sigma_wind * sigma_rain, sigma_rain^2
# ), nrow = 2)
# 
# windrain_data <- mvrnorm(N, mu = a_windrain, Sigma = Sigma_windrain)
# sim_data$Wind <- windrain_data[, 1]
# sim_data$Rain <- windrain_data[, 2]
# 
# sim_data$Micro <- rnorm(N)
# sim_data$Year <- rnorm(N)
# sim_data$Goby_lag <- rnorm(N)
# 
# # --- 4. Generative Cascade (Following the Causal DAG) ---
# 
# ## *** MODIFIED: Add U_phys to BreachDays ***
# mu_breach <- a_BreachDays + 
#   0.5 * sim_data$Rain +  # beta_Rain (not in hierarchical structure)
#   k_U_phys * U_phys[sim_data$Zone]
# sim_data$BreachDays <- rnorm(N, mean = mu_breach, sd = tau)
# 
# ## *** MODIFIED: Add U_phys to Temp ***
# mu_temp <- a_Temp + 
#   beta_BreachDays_vec[2] * sim_data$BreachDays + 
#   beta_Wind_vec[2] * sim_data$Wind + 
#   k_U_phys * U_phys[sim_data$Zone]
# sim_data$Temp <- rnorm(N, mean = mu_temp, sd = tau)
# 
# ## *** MODIFIED: Add U_phys to DO ***
# mu_do <- a_DO + 
#   beta_Temp_vec[2] * sim_data$Temp + 
#   beta_Wind_vec[1] * sim_data$Wind + 
#   k_U_phys * U_phys[sim_data$Zone]
# sim_data$DO <- rnorm(N, mean = mu_do, sd = tau)
# 
# ## *** MODIFIED: Add U_bio to SAV ***
# mu_sav <- a_SAV + 
#   beta_DO_vec[4] * sim_data$DO + 
#   beta_Temp_vec[3] * sim_data$Temp + 
#   k_U_bio * U_bio[sim_data$Zone]
# sim_data$SAV <- rnorm(N, mean = mu_sav, sd = tau)
# 
# ## *** MODIFIED: Add U_bio to SC ***
# logit_sc <- a_SC + 
#   beta_Substrate_vec[2] * sim_data$Substrate + 
#   beta_DO_vec[3] * sim_data$DO + 
#   beta_SAV_vec[3] * sim_data$SAV + 
#   k_U_bio * U_bio[sim_data$Zone]
# sim_data$SC_count <- rbinom(N, size = 1, prob = plogis(logit_sc))
# 
# ## *** MODIFIED: Add U_bio to SB ***
# logit_sb <- a_SB + 
#   beta_DO_vec[2] * sim_data$DO + 
#   beta_SAV_vec[2] * sim_data$SAV + 
#   k_U_bio * U_bio[sim_data$Zone]
# sim_data$SB_count <- rbinom(N, size = 1, prob = plogis(logit_sb))
# 
# # --- 5. Simulate the Final Goby Outcome ---
# log_mu_goby <- (a_Goby_zone[sim_data$Zone] +
#                   sim_data$Area +
#                   beta_SAV_vec[1] * sim_data$SAV +
#                   beta_DO_vec[1] * sim_data$DO +
#                   beta_BreachDays_vec[1] * sim_data$BreachDays +
#                   beta_Temp_vec[1] * sim_data$Temp +
#                   beta_SC_count * sim_data$SC_count +
#                   beta_SB_count * sim_data$SB_count +
#                   beta_Micro * sim_data$Micro +
#                   beta_Substrate_vec[1] * sim_data$Substrate +
#                   beta_Year * sim_data$Year + 
#                   beta_Year_2 * sim_data$Year^2 +
#                   beta_Goby_lag * sim_data$Goby_lag +
#                   beta_BreachDays_2 * sim_data$BreachDays^2 +
#                   beta_Temp_2 * sim_data$Temp^2 +
#                   beta_SAV_2 * sim_data$SAV^2)
# 
# mu_goby <- exp(log_mu_goby)
# sim_data$Goby <- rnbinom(N, mu = mu_goby, size = phi)
# 
# # --- 6. Calculate Goby Density ---
# sim_data$Goby_Density <- sim_data$Goby / exp(sim_data$Area)
# 
# # --- 7. Finalize and Inspect ---
# sim_data$Year_2 <- sim_data$Year^2
# sim_data$Temp_2 <- sim_data$Temp^2
# sim_data$BreachDays_2 <- sim_data$BreachDays^2
# sim_data$SAV_2 <- sim_data$SAV^2
# 
# # Reorder columns
# sim_data <- sim_data[, c("Goby", "Goby_Density", "Zone", "Area", "Rain", "Wind", 
#                          "BreachDays", "Temp", "DO", "SAV", "SC_count", "SB_count", 
#                          "Substrate", "Micro", "Year", "Goby_lag", "Year_2", 
#                          "Temp_2", "BreachDays_2", "SAV_2")]
# 
# # Check the summaries
# cat("--- Summary of Simulated Goby Counts ---\n")
# print(summary(sim_data$Goby))
# cat("\n--- Summary of Simulated Area (log-transformed) ---\n")
# print(summary(sim_data$Area))
# 
# # View the first few rows
# cat("\n--- First 6 Rows of Simulated Dataset ---\n")
# print(head(sim_data))
# 
# # *** NEW: Store true parameter values for comparison ***
# true_params <- list(
#   # Zone-level parameters
#   a_Goby_zone = a_Goby_zone,
#   U_bio = U_bio,
#   U_phys = U_phys,
#   
#   # Latent variable loadings
#   k_U_bio = k_U_bio,
#   k_U_phys = k_U_phys,
#   
#   # Hierarchical coefficient vectors
#   beta_DO_vec = beta_DO_vec,
#   beta_SAV_vec = beta_SAV_vec,
#   beta_Temp_vec = beta_Temp_vec,
#   beta_BreachDays_vec = beta_BreachDays_vec,
#   beta_Wind_vec = beta_Wind_vec,
#   beta_Substrate_vec = beta_Substrate_vec,
#   
#   # Fixed effects
#   beta_Rain = 0.5,  # Not in hierarchical structure
#   beta_SC_count = beta_SC_count,
#   beta_SB_count = beta_SB_count,
#   beta_Micro = beta_Micro,
#   beta_Year = beta_Year,
#   beta_Year_2 = beta_Year_2,
#   beta_Goby_lag = beta_Goby_lag,
#   beta_Temp_2 = beta_Temp_2,
#   beta_BreachDays_2 = beta_BreachDays_2,
#   beta_SAV_2 = beta_SAV_2,
#   
#   # Dispersion and error parameters
#   phi = phi,
#   tau = tau,
#   
#   # Wind/Rain correlation
#   rho_windrain = rho_windrain,
#   sigma_wind = sigma_wind,
#   sigma_rain = sigma_rain
# )
# 
# cat("\n=== TRUE PARAMETER VALUES ===\n")
# cat("\n--- Latent Variables by Zone ---\n")
# cat("U_bio:", round(U_bio, 3), "\n")
# cat("U_phys:", round(U_phys, 3), "\n")
# cat("k_U_bio:", k_U_bio, "\n")
# cat("k_U_phys:", k_U_phys, "\n")
# 
# cat("\n--- Goby Intercepts by Zone ---\n")
# cat("a_Goby:", round(a_Goby_zone, 3), "\n")
# 
# cat("\n--- Hierarchical Coefficients ---\n")
# cat("beta_DO_vec:", round(beta_DO_vec, 3), "\n")
# cat("beta_SAV_vec:", round(beta_SAV_vec, 3), "\n")
# cat("beta_Temp_vec:", round(beta_Temp_vec, 3), "\n")
# 
# cat("\n--- Wind/Rain Correlation ---\n")
# cat("Correlation:", rho_windrain, "\n")
# cat("Actual sample correlation:", round(cor(sim_data$Wind, sim_data$Rain), 3), "\n")
# 
# # ======================================================
# # 8. CALCULATION OF PATH-SPECIFIC CAUSAL EFFECTS
# # ======================================================
# 
# # --- Step A: Calculate Means for Quadratic Derivatives ---
# mean_Temp <- mean(sim_data$Temp)
# mean_Breach <- mean(sim_data$BreachDays)
# mean_SAV <- mean(sim_data$SAV)
# 
# # --- Step B: Calculate "Marginal Slopes" for the Goby Model ---
# slope_Temp_to_Goby <- beta_Temp_vec[1] + (2 * beta_Temp_2 * mean_Temp)
# slope_Breach_to_Goby <- beta_BreachDays_vec[1] + (2 * beta_BreachDays_2 * mean_Breach)  # FIXED
# slope_SAV_to_Goby <- beta_SAV_vec[1] + (2 * beta_SAV_2 * mean_SAV)
# 
# # --- Step C: Calculate Specific Path Effects ---
# path_Breach_Temp_Goby <- beta_BreachDays_vec[2] * slope_Temp_to_Goby
# path_DO_SAV_Goby <- beta_DO_vec[4] * slope_SAV_to_Goby
# path_DO_SB_Goby <- beta_DO_vec[2] * beta_SB_count
# path_DO_SC_Goby <- beta_DO_vec[3] * beta_SC_count
# path_SAV_SB_Goby <- beta_SAV_vec[2] * beta_SB_count
# path_SAV_SC_Goby <- beta_SAV_vec[3] * beta_SC_count
# path_Substrate_SC_Goby <- beta_Substrate_vec[2] * beta_SC_count
# path_Temp_DO_Goby <- beta_Temp_vec[2] * beta_DO_vec[1]
# path_Temp_SAV_Goby <- beta_Temp_vec[3] * slope_SAV_to_Goby
# path_Wind_Temp_Goby <- beta_Wind_vec[2] * slope_Temp_to_Goby
# 
# # --- Step D: Print Results ---
# cat("\n==============================================\n")
# cat(" TRUE PATH-SPECIFIC CAUSAL EFFECTS (Simulated)\n")
# cat(" (Calculated at mean values for quadratic terms)\n")
# cat("==============================================\n")
# 
# results <- data.frame(
#   Path = c("Breach -> Temp -> Goby",
#            "DO -> SAV -> Goby",
#            "DO -> SB -> Goby",
#            "DO -> SC -> Goby",
#            "SAV -> SB -> Goby",
#            "SAV -> SC -> Goby",
#            "Substrate -> SC -> Goby",
#            "Temp -> DO -> Goby",
#            "Temp -> SAV -> Goby",
#            "Wind -> Temp -> Goby"),
#   Effect = c(path_Breach_Temp_Goby,
#              path_DO_SAV_Goby,
#              path_DO_SB_Goby,
#              path_DO_SC_Goby,
#              path_SAV_SB_Goby,
#              path_SAV_SC_Goby,
#              path_Substrate_SC_Goby,
#              path_Temp_DO_Goby,
#              path_Temp_SAV_Goby,
#              path_Wind_Temp_Goby)
# )
# 
# print(results, digits = 4)
# 
# # *** NEW: Function to compare Stan results with true values ***
# compare_with_stan <- function(stan_fit, true_params) {
#   # Extract posterior means from Stan fit
#   posterior_means <- summary(stan_fit)$summary[, "mean"]
#   
#   # Compare latent variables
#   cat("\n=== COMPARISON: Latent Variables ===\n")
#   for(j in 1:J) {
#     cat(sprintf("U_bio[%d]: True = %.3f, Estimated = %.3f\n", 
#                 j, true_params$U_bio[j], posterior_means[paste0("U_bio[", j, "]")]))
#     cat(sprintf("U_phys[%d]: True = %.3f, Estimated = %.3f\n", 
#                 j, true_params$U_phys[j], posterior_means[paste0("U_phys[", j, "]")]))
#   }
#   
#   # Compare loadings
#   cat("\n=== COMPARISON: Loadings ===\n")
#   cat(sprintf("k_U_bio: True = %.3f, Estimated = %.3f\n", 
#               true_params$k_U_bio, posterior_means["k_U_bio"]))
#   cat(sprintf("k_U_phys: True = %.3f, Estimated = %.3f\n", 
#               true_params$k_U_phys, posterior_means["k_U_phys"]))
#   
#   # Add more comparisons as needed...
# }
# 
# # Save everything for later use
# save(sim_data, true_params, file = "simulated_data_with_latents.RData")
# cat("\nSimulated data and true parameters saved to 'simulated_data_with_latents.RData'\n")

## OLD below-


# #adding SAV_2 ----------------------
# library(MASS) # Required for the mvrnorm function
# 
# # --- 1. Set Constants & True Parameter Values ---
# set.seed(42)
# N <- 314 # Total observations
# J <- 3   # Number of zones
# 
# ## Hyperpriors & Goby-specific parameters
# a_bar_Goby <- 2.6
# sigma_Goby <- 0.5
# phi <- 1.0
# 
# ## --- Hierarchical Coefficient Hyperpriors ---
# b_bar_DO <- c(0.6, 0.6, 0.6, 0.6)
# b_bar_SAV <- c(0.4, 0.4, 0.4)
# b_bar_Temp <- c(0.4, 0.4, 0.4)
# b_bar_BreachDays <- c(0.3, 0.3)
# b_bar_Wind <- c(-0.5, -0.5)
# b_bar_Substrate <- c(-0.3, -0.3)
# 
# sigma_vec <- c(0.1, 0.1, 0.1, 0.1)
# rho_val <- 0.3
# 
# make_sigma <- function(k, sd_vec, rho) {
#   rho_mat <- matrix(rho, nrow = k, ncol = k)
#   diag(rho_mat) <- 1
#   sigma_mat <- diag(sd_vec[1:k]) %*% rho_mat %*% diag(sd_vec[1:k])
#   return(sigma_mat)
# }
# 
# Sigma_DO <- make_sigma(4, sigma_vec, rho_val)
# Sigma_SAV <- make_sigma(3, sigma_vec, rho_val)
# Sigma_Temp <- make_sigma(3, sigma_vec, rho_val)
# Sigma_BreachDays <- make_sigma(2, sigma_vec, rho_val)
# Sigma_Wind <- make_sigma(2, sigma_vec, rho_val)
# Sigma_Substrate <- make_sigma(2, sigma_vec, rho_val)
# 
# ## --- Draw the actual, unshared coefficients from these distributions ---
# beta_DO_vec <- mvrnorm(1, b_bar_DO, Sigma_DO)
# beta_SAV_vec <- mvrnorm(1, b_bar_SAV, Sigma_SAV)
# beta_Temp_vec <- mvrnorm(1, b_bar_Temp, Sigma_Temp)
# beta_BreachDays_vec <- mvrnorm(1, b_bar_BreachDays, Sigma_BreachDays)
# beta_Wind_vec <- mvrnorm(1, b_bar_Wind, Sigma_Wind)
# beta_Substrate_vec <- mvrnorm(1, b_bar_Substrate, Sigma_Substrate)
# 
# ## Other fixed coefficients
# beta_Rain <- 0.5
# beta_SC_count <- -0.5
# beta_SB_count <- 0.5
# beta_Micro <- 0.2
# beta_Year <- 0.2
# beta_Year_2 <- -0.2
# beta_Goby_lag <- 0.1
# beta_Temp_2 <- -0.4
# beta_BreachDays_2 <- -0.1
# beta_SAV_2 <- -0.3 # <-- NEW: True effect for SAV squared
# 
# ## Submodel intercepts and error terms
# a_BreachDays <- 0; a_Temp <- 0; a_DO <- 0; a_SAV <- 0; a_SC <- -0.5; a_SB <- -0.2
# tau <- 0.5
# 
# # --- 2. Simulate Zone-Level Parameters ---
# a_Goby_zone <- rnorm(J, mean = a_bar_Goby, sd = sigma_Goby)
# 
# # --- 3. Simulate Observation-Level Data ---
# sim_data <- data.frame(Zone = sample(1:J, size = N, replace = TRUE))
# 
# sim_data$Area <- rnorm(N, mean = 2.4, sd = 0.5) # Already logged
# sim_data$Substrate <- rbinom(N, size = 1, prob = 0.5)
# sim_data$Rain <- rnorm(N)
# sim_data$Wind <- rnorm(N)
# sim_data$Micro <- rnorm(N)
# sim_data$Year <- rnorm(N)
# sim_data$Goby_lag <- rnorm(N)
# 
# # --- 4. Generative Cascade (Following the Causal DAG) --- 
# mu_breach <- a_BreachDays + beta_Rain * sim_data$Rain
# sim_data$BreachDays <- rnorm(N, mean = mu_breach, sd = tau)
# 
# mu_temp <- a_Temp + beta_BreachDays_vec[2] * sim_data$BreachDays + beta_Wind_vec[2] * sim_data$Wind
# sim_data$Temp <- rnorm(N, mean = mu_temp, sd = tau)
# 
# mu_do <- a_DO + beta_Temp_vec[2] * sim_data$Temp + beta_Wind_vec[1] * sim_data$Wind
# sim_data$DO <- rnorm(N, mean = mu_do, sd = tau)
# 
# mu_sav <- a_SAV + beta_DO_vec[4] * sim_data$DO + beta_Temp_vec[3] * sim_data$Temp
# sim_data$SAV <- rnorm(N, mean = mu_sav, sd = tau)
# 
# logit_sc <- a_SC + beta_Substrate_vec[2] * sim_data$Substrate + beta_DO_vec[3] * sim_data$DO + beta_SAV_vec[3] * sim_data$SAV
# sim_data$SC_count <- rbinom(N, size = 1, prob = plogis(logit_sc))
# 
# logit_sb <- a_SB + beta_DO_vec[2] * sim_data$DO + beta_SAV_vec[2] * sim_data$SAV
# sim_data$SB_count <- rbinom(N, size = 1, prob = plogis(logit_sb))
# 
# # --- 5. Simulate the Final Goby Outcome --- 
# log_mu_goby <- (a_Goby_zone[sim_data$Zone] +
#                   sim_data$Area + # Area is the offset
#                   beta_SAV_vec[1] * sim_data$SAV +
#                   beta_DO_vec[1] * sim_data$DO +
#                   beta_BreachDays_vec[1] * sim_data$BreachDays +
#                   beta_Temp_vec[1] * sim_data$Temp +
#                   beta_SC_count * sim_data$SC_count +
#                   beta_SB_count * sim_data$SB_count +
#                   beta_Micro * sim_data$Micro +
#                   beta_Substrate_vec[1] * sim_data$Substrate +
#                   beta_Year * sim_data$Year + beta_Year_2 * sim_data$Year^2 +
#                   beta_Goby_lag * sim_data$Goby_lag +
#                   beta_BreachDays_2 * sim_data$BreachDays^2 +
#                   beta_Temp_2 * sim_data$Temp^2 +
#                   beta_SAV_2 * sim_data$SAV^2) # <-- NEW: SAV^2 term added
# 
# mu_goby <- exp(log_mu_goby)
# sim_data$Goby <- rnbinom(N, mu = mu_goby, size = phi)
# 
# # --- 6. Calculate Goby Density ---
# sim_data$Goby_Density <- sim_data$Goby / exp(sim_data$Area) # exp() to get back to natural scale
# 
# # --- 7. Finalize and Inspect ---
# sim_data$Year_2 <- sim_data$Year^2
# sim_data$Temp_2 <- sim_data$Temp^2
# sim_data$BreachDays_2 <- sim_data$BreachDays^2
# sim_data$SAV_2 <- sim_data$SAV^2 # <-- NEW: Add SAV_2 column to final data
# 
# # Reorder columns
# sim_data <- sim_data[, c("Goby", "Goby_Density", "Zone", "Area", "Rain", "Wind", "BreachDays",
#                          "Temp", "DO", "SAV", "SC_count", "SB_count", "Substrate",
#                          "Micro", "Year", "Goby_lag", "Year_2", "Temp_2", "BreachDays_2", "SAV_2")] # <-- NEW: SAV_2 added to reorder
# 
# # Check the summaries
# cat("--- Summary of Simulated Goby Counts ---\n")
# print(summary(sim_data$Goby))
# cat("\n--- Summary of Simulated Area (log-transformed) ---\n")
# print(summary(sim_data$Area))
# 
# # View the first few rows
# cat("\n--- First 6 Rows of Simulated Dataset ---\n")
# print(head(sim_data))
# 
# 
# # ======================================================
# # 8. CALCULATION OF PATH-SPECIFIC CAUSAL EFFECTS
# # ======================================================
# 
# # --- Step A: Calculate Means for Quadratic Derivatives ---
# # We need the average environmental conditions to evaluate the curves
# mean_Temp <- mean(sim_data$Temp)
# mean_Breach <- mean(sim_data$BreachDays)
# mean_SAV <- mean(sim_data$SAV)
# 
# # --- Step B: Calculate "Marginal Slopes" for the Goby Model ---
# # Since Goby has quadratic inputs (Temp^2, Breach^2, SAV^2), the "Slope" 
# # is the derivative: beta + 2*beta_2*X
# 
# # Slope of Temp --> Goby
# slope_Temp_to_Goby <- beta_Temp_vec[1] + (2 * beta_Temp_2 * mean_Temp)
# 
# # Slope of Breach --> Goby
# slope_Breach_to_Goby <- beta_BreachDays_vec[1] + (2 * beta_BreachDays_2 * mean_Temp)
# 
# # Slope of SAV --> Goby
# slope_SAV_to_Goby <- beta_SAV_vec[1] + (2 * beta_SAV_2 * mean_SAV)
# 
# # --- Step C: Calculate Specific Path Effects ---
# # We multiply the coefficients along the chain.
# 
# # 1. BreachDays --> Temp --> Goby
# # Path: (Breach->Temp) * (Temp->Goby)
# path_Breach_Temp_Goby <- beta_BreachDays_vec[2] * slope_Temp_to_Goby
# 
# # 2. DO --> SAV --> Goby (Assuming 'SAC' in your prompt was a typo for SAV)
# # Path: (DO->SAV) * (SAV->Goby)
# # Note: beta_DO_vec[4] is the effect of DO on SAV in your code
# path_DO_SAV_Goby <- beta_DO_vec[4] * slope_SAV_to_Goby
# 
# # 3. DO --> SB --> Goby
# # Path: (DO->SB_logit) * (SB->Goby)
# # Note: beta_DO_vec[2] is DO effect on SB
# path_DO_SB_Goby <- beta_DO_vec[2] * beta_SB_count
# 
# # 4. DO --> SC --> Goby
# # Path: (DO->SC_logit) * (SC->Goby)
# # Note: beta_DO_vec[3] is DO effect on SC
# path_DO_SC_Goby <- beta_DO_vec[3] * beta_SC_count
# 
# # 5. SAV --> SB --> Goby
# # Path: (SAV->SB_logit) * (SB->Goby)
# path_SAV_SB_Goby <- beta_SAV_vec[2] * beta_SB_count
# 
# # 6. SAV --> SC --> Goby
# # Path: (SAV->SC_logit) * (SC->Goby)
# path_SAV_SC_Goby <- beta_SAV_vec[3] * beta_SC_count
# 
# # 7. Substrate --> SC --> Goby
# # Path: (Substrate->SC_logit) * (SC->Goby)
# path_Substrate_SC_Goby <- beta_Substrate_vec[2] * beta_SC_count
# 
# # 8. Temp --> DO --> Goby
# # Path: (Temp->DO) * (DO->Goby)
# # Note: beta_Temp_vec[2] is Temp on DO; beta_DO_vec[1] is DO on Goby
# path_Temp_DO_Goby <- beta_Temp_vec[2] * beta_DO_vec[1]
# 
# # 9. Temp --> SAV --> Goby
# # Path: (Temp->SAV) * (SAV->Goby)
# # Note: beta_Temp_vec[3] is Temp on SAV
# path_Temp_SAV_Goby <- beta_Temp_vec[3] * slope_SAV_to_Goby
# 
# # 10. Wind --> Temp --> Goby
# # Path: (Wind->Temp) * (Temp->Goby)
# path_Wind_Temp_Goby <- beta_Wind_vec[2] * slope_Temp_to_Goby
# 
# 
# # --- Step D: Print Results ---
# cat("\n==============================================\n")
# cat(" TRUE PATH-SPECIFIC CAUSAL EFFECTS (Simulated)\n")
# cat(" (Calculated at mean values for quadratic terms)\n")
# cat("==============================================\n")
# results <- data.frame(
#   Path = c("Breach -> Temp -> Goby",
#            "DO -> SAV -> Goby",
#            "DO -> SB -> Goby",
#            "DO -> SC -> Goby",
#            "SAV -> SB -> Goby",
#            "SAV -> SC -> Goby",
#            "Substrate -> SC -> Goby",
#            "Temp -> DO -> Goby",
#            "Temp -> SAV -> Goby",
#            "Wind -> Temp -> Goby"),
#   Effect = c(path_Breach_Temp_Goby,
#              path_DO_SAV_Goby,
#              path_DO_SB_Goby,
#              path_DO_SC_Goby,
#              path_SAV_SB_Goby,
#              path_SAV_SC_Goby,
#              path_Substrate_SC_Goby,
#              path_Temp_DO_Goby,
#              path_Temp_SAV_Goby,
#              path_Wind_Temp_Goby)
# )
# 
# print(results, digits = 4)


