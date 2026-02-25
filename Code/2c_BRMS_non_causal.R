# Create a simulated dataset with known coefficients
# compare glmer coefficients vs shared estimate coefficients
# then compare a do calculus estimate door estimate to shared model.





#fit a brms model to compare to the causal model





library(brms)
library(lme4)
#build a full model per the DAG and latest 2025 ulam model

m1.glmer <- glmer(Goby ~
  Year +
  Year_2 +
  SC_count + 
  Rain + 
  SAV + 
  SAV_2 +
  SB_count + 
  DO +
  Micro + 
  BreachDays + 
  BreachDays_2 + 
  Substrate +
  Wind +
  Temp +   
  Temp_2 +  
  Goby_lag + 
  offset(Area) + 
    (1 | Zone),
  family = negative.binomial(2),
  data = dat) #from 2_Full_Luxury_Gobies.R
  
summary(m1.glmer)
plot_model(m1.glmer, sort.est = TRUE)
  

m1.brm <- brm(Goby ~
                    Year +
                    Year_2 +
                    SC_count + 
                    #Rain + 
                    SAV + 
                    SAV_2 +
                    SB_count + 
                    DO +
                    Micro + 
                    BreachDays + 
                    BreachDays_2 + 
                    Substrate +
                    Wind +
                    Temp +   
                    Temp_2 +  
                    Goby_lag + 
                    offset(Area) + 
                    (1 | Zone),
                  family = negbinomial(link = "log", link_shape = "log"),
                  data = dat,
                   backend = getOption("brms.backend", "rstan")
                   ) #from 2_Full_Luxury_Gobies.R

summary(m1.brm) # glmer coefficients  plot vs network coefficients

library(sjPlot)

plot_model(m1.brm)






#Simulate dataset
# Load necessary library for the negative binomial distribution
library(MASS)
library(dplyr)

# --- 1. Set Simulation Parameters ---

# Set a seed for reproducibility
set.seed(42)

# Number of observations (rows) in the dataset
n <- 300

# Define the random intercepts for the three zones
zone_intercepts <- c("Zone1" = -0.16, "Zone2" = 0.65, "Zone3" = 0.79)
zone_names <- names(zone_intercepts)

# Define the overdispersion parameter for the negative binomial distribution.
# A lower value means more overdispersion.
theta <- 2


# --- 2. Simulate Independent Variables (Covariates) ---

# Zone (random effect grouping factor)
Zone <- sample(zone_names, n, replace = TRUE)

# Year (scaled to mean=0, sd=1)
Year <- rnorm(n, mean = 0, sd = 1)
Year_2 <- Year^2

# SC_count (presence/absence)
SC_count <- rbinom(n, size = 1, prob = 0.45) # 45% presence

# SAV (ordinal, 5 classes from 0-4)
SAV <- sample(0:4, n, replace = TRUE)

# SB_count (presence/absence)
SB_count <- rbinom(n, size = 1, prob = 0.55) # 55% presence

# DO (Dissolved Oxygen, Gaussian)
DO <- rnorm(n, mean = 0, sd = 1)

# Micro (Microbiota count, Poisson)
Micro <- rpois(n, lambda = 2)

# BreachDays (Poisson)
BreachDays <- rnorm(n, mean = 0, sd = 1)
BreachDays_2 <- BreachDays^2

# Substrate (binary, e.g., sandy vs. rocky)
Substrate <- rbinom(n, size = 1, prob = 0.6) # 60% one type

# Wind (continuous)
Wind <- rnorm(n, mean = 0, sd = 1)

# Temp (continuous)
Temp <- rnorm(n, mean = 0, sd = 1)
Temp_2 <- Temp^2

# Area (for the offset term, must be > 0)
Area <- runif(n, min = 1, max = 10)


# --- 3. Construct the Linear Predictor (eta) ---

# Get the random intercept value for each observation based on its Zone
random_effect <- zone_intercepts[Zone]

# Calculate the linear predictor by multiplying each covariate by its coefficient
# The offset is added directly on the log scale
eta <- (
  -1.52 * Year +
    1.52 * Year_2 +
    -0.82 * SC_count +
    0.19 * SAV +
    0.67 * SB_count +
    0.06 * DO +
    0.33 * Micro +
    0.08 * BreachDays +
    0.01 * BreachDays_2 +
    0.26 * Substrate +
    -0.02 * Wind +
    0.18 * Temp +
    -0.37 * Temp_2 +
    offset(Area) +        # Offset term
    random_effect      # Random intercept
)

# --- 4. Simulate the Dependent Variable (Goby count) ---

# Convert the linear predictor to the mean (mu) scale using the inverse link function (exp)
mu <- exp(eta)

# Generate the Goby counts from a negative binomial distribution
Goby <- rnbinom(n = n, size = theta, mu = mu)


# --- 5. Assemble the Final Dataset ---
simulated_data <- data.frame(
  Goby,
  Year,
  Year_2,
  SC_count,
  SAV,
  SB_count,
  DO,
  Micro,
  BreachDays,
  BreachDays_2,
  Substrate,
  Wind,
  Temp,
  Temp_2,
  Area,
  Zone
)

# --- 6. View the Result ---
# Display the first few rows of the dataset
print(head(simulated_data))

# run glmm

m1.simulated.glmer <- glmer(Goby ~
                    Year +
                    Year_2 +
                    SC_count + 
                    #Rain + 
                    SAV + 
                    #SAV_2 +
                    SB_count + 
                    DO +
                    Micro + 
                    BreachDays + 
                    BreachDays_2 + 
                    Substrate +
                    Wind +
                    Temp +   
                    Temp_2 +  
                    #Goby_lag + 
                    offset(Area) + 
                    (1 | Zone),
                  family = negative.binomial(2),
                  data = simulated_data) #from 2_Full_Luxury_Gobies.R

summary(m1.simulated.glmer)


## glm simulator
library(stats)
library(MASS)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(GlmSimulatoR)

set.seed(1)
simdata <- simulate_negative_binomial(
  N = 314, weights = c(
   -1.52,# * Year +
    1.52,# * Year_2 +
   -0.82,# * SC_count +
    0.19,#* SAV +
    0.67,# * SB_count +
    0.06,# * DO +
    0.33,# * Micro +
    0.08,#,* BreachDays +
    0.01,# * BreachDays_2 +
    0.26,# * Substrate +
   -0.02,# * Wind +
    0.18,# * Temp +
   -0.37),# * Temp_2 +
    #ignore offset
  
  ancillary = 2, # ancillary is theta.
  link = "log"
)

# Response looks like a negative binomial distribution.
ggplot(simdata, aes(Y)) +
  geom_histogram(bins = 100)

glm_nb <- glm.nb(Y ~ X1 + X2 + X3 + X4 +X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13, data = simdata, link = "log")
summary



## another try
# Load libraries for data manipulation and the negative binomial distribution
library(dplyr)
library(MASS)

# For reproducibility of the simulation
set.seed(42)

## --------------------------------------
## 1. DRAW "TRUE" PARAMETERS FROM PRIORS
## --------------------------------------
# This mimics the generative process by picking one possible reality from the priors.
n_obs <- 314
n_zones <- 3

# Draw global parameters
phi <- rnorm(1, 1, 5)
tau <- rexp(1, 1)

# Draw hyperparameters for the Zone random effect
mu_Zone <- rnorm(1, 0, 0.5)
tau_Zone <- rexp(1, 1)

# Draw Zone-specific random intercepts (a_Goby)
a_Goby <- rnorm(n_zones, mu_Zone, tau_Zone)

# Draw intercepts for submodels
a_Temp <- rnorm(1, 0, 0.5)
a_SAV <- rnorm(1, 0, 0.5)
a_SC <- rnorm(1, 0, 0.5)
a_SB <- rnorm(1, 0, 0.5)
a_DO <- rnorm(1, 0, 0.5)
a_BreachDays <- rnorm(1, 0, 0.5)

# Draw coefficients (betas)
beta_Goby_lag <- rnorm(1, 0.25, 0.25)
beta_Substrate <- rnorm(1, 0.25, 0.25)
beta_BreachDays_2 <- rnorm(1, -0.1, 0.25)
beta_BreachDays <- rnorm(1, 0.25, 0.25)
beta_DO <- rnorm(1, 0.25, 0.25)
beta_SAV_2 <- rnorm(1, -0.1, 0.25) # Note: Not used in the model block, but included for completeness
beta_SAV <- rnorm(1, 0, 0.25)
beta_SC_count <- rnorm(1, -0.1, 0.25)
beta_Rain <- rnorm(1, 0.25, 0.25)
beta_Micro <- rnorm(1, 0.25, 0.25)
beta_Temp_2 <- rnorm(1, -0.1, 0.25)
beta_SB_count <- rnorm(1, 0, 0.5)
beta_Year <- rnorm(1, 0, 0.5)
beta_Year_2 <- rnorm(1, 0, 0.5)
beta_Temp <- rnorm(1, 0, 0.5)
beta_Wind <- rnorm(1, 0, 0.5)


## --------------------------------------
## 2. SIMULATE DATA IN ORDER OF DEPENDENCY
## --------------------------------------

# Create a data frame and simulate exogenous predictors (no dependencies)
sim_data <- tibble(
  Zone = sample(1:n_zones, n_obs, replace = TRUE),
  Year = scale(2000:(2000 + n_obs - 1))[,1],
  Micro = rpois(n_obs, lambda = 15),
  Substrate = rbinom(n_obs, 1, 0.5),
  Rain = rnorm(n_obs, 10, 5),
  Wind = rnorm(n_obs, 15, 4),
  Area = log(runif(n_obs, 10, 50)) # Area is the offset, so simulate on log scale
) %>%
  mutate(Year_2 = Year^2)

# --- Simulate variables from submodels in sequence ---

# BreachDays depends on Rain
Breach_nu <- a_BreachDays + beta_Rain * sim_data$Rain
sim_data$BreachDays <- rnorm(n_obs, mean = Breach_nu, sd = tau)
sim_data$BreachDays_2 <- sim_data$BreachDays^2

# Temp depends on BreachDays and Wind
Temp_nu <- a_Temp + beta_BreachDays * sim_data$BreachDays + beta_Wind * sim_data$Wind
sim_data$Temp <- rnorm(n_obs, mean = Temp_nu, sd = tau)
sim_data$Temp_2 <- sim_data$Temp^2

# DO depends on Temp and Wind
DO_nu <- a_DO + beta_Temp * sim_data$Temp + beta_Wind * sim_data$Wind
sim_data$DO <- rnorm(n_obs, mean = DO_nu, sd = tau)

# SAV depends on DO and Temp
SAV_nu <- a_SAV + beta_DO * sim_data$DO + beta_Temp * sim_data$Temp
sim_data$SAV <- rnorm(n_obs, mean = SAV_nu, sd = tau)
# The Stan code has SAV_2 in data but not model; we'll omit it from the final data.

# SC_count depends on Substrate, DO, SAV
SC_mu_linear <- a_SC + beta_Substrate * sim_data$Substrate + beta_DO * sim_data$DO + beta_SAV * sim_data$SAV
SC_prob <- 1 / (1 + exp(-SC_mu_linear)) # Inverse logit
sim_data$SC_count <- rbinom(n_obs, 1, SC_prob)

# SB_count depends on DO, SAV
SB_mu_linear <- a_SB + beta_DO * sim_data$DO + beta_SAV * sim_data$SAV
SB_prob <- 1 / (1 + exp(-SB_mu_linear)) # Inverse logit
sim_data$SB_count <- rbinom(n_obs, 1, SB_prob)


## --------------------------------------
## 3. SEQUENTIALLY SIMULATE THE FINAL OUTCOME (Goby)
## --------------------------------------
# This must be done in a loop because of the Goby_lag term

# Initialize empty columns
sim_data$Goby <- NA_integer_
sim_data$Goby_lag <- NA_real_

for (i in 1:n_obs) {
  # The first observation has no lag, so we'll set its effect to 0
  lag_effect <- if (i == 1) 0 else beta_Goby_lag * sim_data$Goby[i - 1]
  
  # Set the lag value in the data frame (for row 2 onwards)
  if (i > 1) sim_data$Goby_lag[i] <- sim_data$Goby[i - 1]
  
  # Calculate the linear predictor for mu
  mu_linear <- with(sim_data[i, ], {
    a_Goby[Zone] +
      beta_Year * Year +
      beta_Year_2 * Year_2 + # Corrected from Stan code which had a bug
      beta_SC_count * SC_count +
      beta_SAV * SAV +
      beta_SB_count * SB_count +
      beta_DO * DO +
      beta_Micro * Micro +
      beta_BreachDays * BreachDays +
      beta_BreachDays_2 * BreachDays_2 +
      beta_Substrate * Substrate +
      beta_Wind * Wind +
      beta_Temp * Temp +
      beta_Temp_2 * Temp_2 +
      lag_effect +
      Area # The offset is added directly
  })
  
  # Apply inverse link to get the mean
  mu <- exp(mu_linear)
  
  # Simulate the Goby count from the negative binomial distribution
  # The second parameter of neg_binomial_2 is the 'precision' or 'shape'
  shape_param <- exp(phi)
  sim_data$Goby[i] <- rnbinom(1, mu = mu, size = shape_param)
}

## --------------------------------------
## 4. FINALIZE AND VIEW DATA
## --------------------------------------

# Select the columns that would be passed to Stan
final_simulated_data <- sim_data %>%
  select(
    Goby, Zone, Year, Year_2, SC_count, SB_count, SAV, DO, Micro,
    BreachDays, BreachDays_2, Substrate, Wind, Temp, Temp_2,
    Goby_lag, Rain, Area
  )

print(head(final_simulated_data))

