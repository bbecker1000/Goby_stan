# --- 1. Set Constants & True Parameter Values ---
set.seed(42)
N <- 314 # Total observations
J <- 3   # Number of zones

## Hyperpriors & Goby-specific parameters
mu_Zone <- 3.5 ## CHANGED: Recalibrated to target a Goby mean of ~100
tau_Zone <- 0.5
phi <- 0.5 # A smaller phi increases variance, allowing for max counts ~1000

## Coefficients provided by the user
beta_Rain <- 0.5
beta_Wind <- -0.5
beta_Substrate <- -0.3
beta_Temp <- 0.4
beta_Temp_2 <- -0.4
beta_DO <- 0.6
beta_SAV <- 0.4
beta_SAV_2 <- -0.2
beta_SC_count <- -0.5
beta_SB_count <- 0.5
beta_BreachDays <- 0.3
beta_BreachDays_2 <- -0.1
beta_Micro <- 0.2
beta_Year <- 0.2
beta_Year_2 <- -0.2
beta_Goby_lag <- 0.1

## Submodel intercepts and error terms
a_BreachDays <- 0
a_Temp <- 0
a_DO <- 0
a_SAV <- 0
a_SC <- -0.5
a_SB <- -0.2
tau <- 0.5

# --- 2. Simulate Zone-Level Parameters ---
a_Goby_zone <- rnorm(J, mean = mu_Zone, sd = tau_Zone)

# --- 3. Simulate Observation-Level Data ---
sim_data <- data.frame(
  Zone = sample(1:J, size = N, replace = TRUE)
)

# Simulate exogenous predictors
sim_data$Area <- rlnorm(N, meanlog = 1.109, sdlog = 1) # Area with mean ~5
sim_data$Substrate <- rbinom(N, size = 1, prob = 0.5)
sim_data$Rain <- rnorm(N)
sim_data$Wind <- rnorm(N)
sim_data$Micro <- rnorm(N)
sim_data$Year <- rnorm(N)
sim_data$Goby_lag <- rnorm(N)

# --- 4. Generative Cascade (Following the Causal DAG) --- 🌊
mu_breach <- a_BreachDays + beta_Rain * sim_data$Rain
sim_data$BreachDays <- rnorm(N, mean = mu_breach, sd = tau)

mu_temp <- a_Temp + beta_BreachDays * sim_data$BreachDays + beta_Wind * sim_data$Wind
sim_data$Temp <- rnorm(N, mean = mu_temp, sd = tau)

mu_do <- a_DO + beta_Temp * sim_data$Temp + beta_Temp_2 * sim_data$Temp^2
sim_data$DO <- rnorm(N, mean = mu_do, sd = tau)

mu_sav <- a_SAV + beta_DO * sim_data$DO + beta_Temp * sim_data$Temp
sim_data$SAV <- rnorm(N, mean = mu_sav, sd = tau)

logit_sc <- a_SC + beta_Substrate * sim_data$Substrate + beta_DO * sim_data$DO + beta_SAV * sim_data$SAV
prob_sc <- plogis(logit_sc)
sim_data$SC_count <- rbinom(N, size = 1, prob = prob_sc)

logit_sb <- a_SB + beta_DO * sim_data$DO + beta_SAV * sim_data$SAV
prob_sb <- plogis(logit_sb)
sim_data$SB_count <- rbinom(N, size = 1, prob = prob_sb)

# --- 5. Simulate the Final Goby Outcome --- 🐟
log_mu_goby <- (a_Goby_zone[sim_data$Zone] +
                  log(sim_data$Area) +
                  beta_SAV * sim_data$SAV + beta_SAV_2 * sim_data$SAV^2 +
                  beta_DO * sim_data$DO +
                  beta_BreachDays * sim_data$BreachDays + beta_BreachDays_2 * sim_data$BreachDays^2 +
                  beta_Temp * sim_data$Temp + beta_Temp_2 * sim_data$Temp^2 +
                  beta_SC_count * sim_data$SC_count +
                  beta_SB_count * sim_data$SB_count +
                  beta_Micro * sim_data$Micro +
                  beta_Substrate * sim_data$Substrate +
                  beta_Year * sim_data$Year + beta_Year_2 * sim_data$Year^2 +
                  beta_Goby_lag * sim_data$Goby_lag)

mu_goby <- exp(log_mu_goby)
sim_data$Goby <- rnbinom(N, mu = mu_goby, size = phi)

# --- 6. Calculate Goby Density ---
sim_data$Goby_Density <- sim_data$Goby / sim_data$Area

# --- 7. Finalize and Inspect ---
sim_data$Year_2 <- sim_data$Year^2
sim_data$Temp_2 <- sim_data$Temp^2
sim_data$BreachDays_2 <- sim_data$BreachDays^2
sim_data$SAV_2 <- sim_data$SAV^2

sim_data <- sim_data[, c("Goby", "Goby_Density", "Zone", "Rain", "Wind", "BreachDays", "Temp", "DO", "SAV",
                         "SC_count", "SB_count", "Substrate", "Micro", "Year",
                         "Goby_lag", "Area", "Year_2", "Temp_2", "BreachDays_2", "SAV_2")]

cat("--- Summary of Simulated Goby Counts ---\n")
print(summary(sim_data$Goby))

cat("\n--- Summary of Simulated Area ---\n")
print(summary(sim_data$Area))

cat("\n--- First 6 Rows of Simulated Dataset ---\n")
print(head(sim_data))


hist(sim_data$Goby_Density)

simulated_data <- sim_data


##model##modelGoby_Density
library(cmdstanr)
library(posterior)

# Compile
#write model
f.mod <- write_stan_file(goby_model.stan)

# compile model 
# about 3 minutes
mod <- cmdstan_model(f.mod)


# Simulate or use your data frame
data_list <- list(
  N = nrow(simulated_data),
  J = length(unique(simulated_data$Zone)),
  Zone = as.integer(as.factor(simulated_data$Zone)),
  Goby = simulated_data$Goby,
  BreachDays = simulated_data$BreachDays,
  DO = simulated_data$DO,
  Temp = simulated_data$Temp,
  SAV = simulated_data$SAV,
  SB_count = simulated_data$SB_count,
  SC_count = simulated_data$SC_count,
  Wind = simulated_data$Wind,
  Rain = simulated_data$Rain,
  Area = simulated_data$Area,
  Substrate = simulated_data$Substrate,
  Micro = simulated_data$Micro,
  Year = simulated_data$Year,
  Year_2 = simulated_data$Year_2,
  Temp_2 = simulated_data$Temp_2,
  Goby_lag = simulated_data$Goby_lag,
  BreachDays_2 = simulated_data$BreachDays_2
)

# Fit
fit <- mod$sample(
  data = data_list,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 3000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)

fit$summary()
fit$diagnostic_summary()

# for saving and seeing
print(fit$summary(), n = 90)
fit$metadata()$model_params
# print only
fit$cmdstan_summary()



