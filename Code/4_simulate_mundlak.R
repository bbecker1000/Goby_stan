# Help from gemini
# Title: R Script for Simulation with Final Area Range
# Description: This is the final script where Area is a uniform value
# between 1 and 3, and all model dependencies are correctly specified.

# --- 1. Simulation Setup ---
N <- 314 # Number of observations

# Define an inverse logit function
inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}


# --- 2. Set Known Coefficients ---
# --- CHANGE 1: Adjust intercept to compensate for smaller Area values ---
mu_Zone <- 4.0 # Increased from 3.0 to keep Goby counts in the target range
beta_Area <- 0.01

# Other parameters
tau_Zone <- 0.5
a_BreachDays <- 0
a_Temp <- 0
a_DO <- 0
a_SAV <- 0
a_SC <- -0.5
a_SB <- -0.2
beta_Rain <- 0.5
beta_Wind <- -0.5
beta_Substrate <- -0.3
beta_Temp <- 0.4
beta_Temp_2 <- -0.4
beta_DO <- 0.6
beta_SAV <- 0.4
beta_SAV_2 <- -0.2
beta_SC_count <- -0.5
beta_SB_count <-  0.5
beta_BreachDays <- 0.3
beta_BreachDays_2 <- -0.1
beta_Micro <- 0.2
beta_Year <- 0.2
beta_Year_2 <- -0.2
beta_Goby_lag <- 0.1
tau <- 0.5
phi <- 1.0


# --- 3. Generate Base Predictor Variables ---
predictors <- data.frame(
  # --- Unscaled Predictors ---
  Zone = sample(1:3, N, replace = TRUE),
  Substrate = rbinom(N, 1, 0.5),
  # --- CHANGE 2: Generate Area from a uniform 1-3 range ---
  Area = runif(N, min = 1, max = 3),
  
  # --- Scaled Predictors ---
  Rain = scale(rnorm(N, mean = 10, sd = 3)),
  Wind = scale(rnorm(N, mean = 5, sd = 2)),
  Micro = scale(rnorm(N, mean = 1, sd = 0.5)),
  Year = scale(sample(2000:2020, N, replace = TRUE)),
  Goby_lag = scale(rpois(N, lambda = 50))
)


# --- 4. Simulate Outcome Variables Sequentially ---
sim_outcomes <- data.frame(matrix(NA, nrow = N, ncol = 0))

# 4a. Simulate base continuous variables
breach_nu <- a_BreachDays + beta_Rain * predictors$Rain
sim_outcomes$BreachDays <- scale(rnorm(N, mean = breach_nu, sd = tau))

temp_nu <- a_Temp + beta_BreachDays * sim_outcomes$BreachDays + beta_Wind * predictors$Wind
sim_outcomes$Temp <- scale(rnorm(N, mean = temp_nu, sd = tau))

do_nu <- a_DO + beta_Temp * sim_outcomes$Temp + beta_Wind * predictors$Wind
sim_outcomes$DO <- scale(rnorm(N, mean = do_nu, sd = tau))

sav_nu <- a_SAV + beta_DO * sim_outcomes$DO + beta_Temp * sim_outcomes$Temp
sim_outcomes$SAV <- scale(rnorm(N, mean = sav_nu, sd = tau))

# 4b. Simulate binary outcomes
sc_mu_logit <- a_SC + beta_Substrate * predictors$Substrate + beta_DO * sim_outcomes$DO + beta_SAV * sim_outcomes$SAV
sc_prob <- inv_logit(sc_mu_logit)
sim_outcomes$SC_count <- rbinom(N, size = 1, prob = sc_prob)

sb_mu_logit <- a_SB + beta_DO * sim_outcomes$DO + beta_SAV * sim_outcomes$SAV
sb_prob <- inv_logit(sb_mu_logit)
sim_outcomes$SB_count <- rbinom(N, size = 1, prob = sb_prob)

# 4c. Create squared terms from the simulated and base variables
predictors$Year_2 <- predictors$Year^2
sim_outcomes$BreachDays_2 <- sim_outcomes$BreachDays^2
sim_outcomes$Temp_2 <- sim_outcomes$Temp^2
sim_outcomes$SAV_2 <- sim_outcomes$SAV^2

# 4d. Simulate final Goby count with all predictors
a_Goby <- rnorm(3, mean = mu_Zone, sd = tau_Zone)
goby_log_mu <- (
  a_Goby[predictors$Zone]
  + beta_Year * predictors$Year
  + beta_SC_count * sim_outcomes$SC_count
  + beta_SAV * sim_outcomes$SAV
  + beta_SB_count * sim_outcomes$SB_count
  + beta_DO * sim_outcomes$DO
  + beta_Micro * predictors$Micro
  + beta_BreachDays * sim_outcomes$BreachDays
  + beta_Substrate * predictors$Substrate
 # + beta_Wind * predictors$Wind
  + beta_Temp * sim_outcomes$Temp
  + beta_Goby_lag * predictors$Goby_lag
  + beta_Area * predictors$Area
  + beta_Year_2 * predictors$Year_2
  + beta_SAV_2 * sim_outcomes$SAV_2
  + beta_BreachDays_2 * sim_outcomes$BreachDays_2
  + beta_Temp_2 * sim_outcomes$Temp_2
)
goby_mu <- exp(goby_log_mu)
goby_raw_counts <- rnbinom(N, mu = goby_mu, size = exp(phi))

# Enforce density cap
goby_max_allowed <- floor(75 * predictors$Area)
sim_outcomes$Goby <- pmin(goby_raw_counts, goby_max_allowed)


# --- 5. Finalize and Verify ---
simulated_data <- cbind(predictors, sim_outcomes)

print("--- First 6 rows of the final dataset ---")
head(simulated_data)


print("Summary of final Goby counts:")
summary(simulated_data$Goby)

write.csv(simulated_data, "Output/Data/simulated_data.csv", row.names = FALSE)

