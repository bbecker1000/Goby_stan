# ZONE-LEVEL LATENT ENVIRONMENT MODEL (FEASIBLE)
# Models zone-average "true" Temp and DO
# Individual observations deviate from zone means
# Strong informative priors prevent attenuation

 

# ZONE-LEVEL LATENT ENVIRONMENT MODEL (FEASIBLE)
# Models zone-average "true" Temp and DO
# Individual observations deviate from zone means
# Strong informative priors prevent attenuation

stan_SEM_code <- "data {
  int<lower=1> N;
  int<lower=1> J;
  
  // Responses and predictors
  array[N] int Goby;
  vector[N] Area;
  vector[N] BreachDays_2;
  vector[N] Goby_lag;
  vector[N] Micro;
  vector[N] Year_2;
  vector[N] Year;
  vector[N] Rain;
  vector[N] Wind;
  vector[N] BreachDays;
  array[N] int Substrate;
  vector[N] SAV;
  vector[N] SAV_2;
  array[N] int Zone;
  
  // MEASURED (noisy) environmental variables
  vector[N] Temp_measured;
  vector[N] DO_measured;
  
  // MEASURED biological variables (no submodels)
  array[N] int SC_count;
  array[N] int SB_count;
}
parameters {
  // Goby varying intercepts
  vector[J] a_Goby;
  real a_bar_Goby;
  real<lower=0> sigma_Goby;
  
  // Latent variable for unmeasured physical confounding
  vector[J] U_phys;
  
  // *** ZONE-LEVEL LATENT ENVIRONMENT (FEASIBLE!) ***
  vector[J] Temp_true_zone;  // True average temp per zone
  vector[J] DO_true_zone;    // True average DO per zone
  
  // Within-zone variation
  real<lower=0> sigma_Temp_within;
  real<lower=0> sigma_DO_within;
  
  // Measurement error (spatial/temporal sampling)
  real<lower=0> sigma_Temp_measure;
  real<lower=0> sigma_DO_measure;
  
  // Correlated weather errors
  corr_matrix[2] Rho_weather;
  vector<lower=0>[2] sigma_weather;
  vector[2] a_weather;
  
  // Pathway coefficients
  real beta_Rain_Breach;
  real beta_Breach_Temp;
  real beta_Wind_Temp;
  real beta_Wind_DO;
  
  // *** KEY PARAMETERS - STRONG PRIORS ***
  real beta_Temp_DO;
  real beta_Temp_Goby;
  real beta_DO_Goby;
  
  // Other effects
  real beta_BreachDays_Goby;
  real beta_SAV_Goby;
  real beta_SC_count;
  real beta_SB_count;
  real beta_Substrate_Goby;
  real beta_Micro;
  real beta_Year;
  real beta_Year_2;
  real beta_Goby_lag;
  
  // Quadratic terms
  real beta_BreachDays_2;
  real beta_Temp_2;
  real beta_SAV_2;
  
  // Submodel intercepts
  real a_BreachDays;
  real a_Temp;
  real a_DO;
  
  // Latent variable loading
  real<lower=0> k_U_phys;
  
  // Error terms
  real<lower=0> tau;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] Temp_true;  // Observation-level true temperature
  vector[N] DO_true;    // Observation-level true DO
  vector[N] Temp_true_2;
  
  // Each observation gets zone mean
  for (i in 1:N) {
    Temp_true[i] = Temp_true_zone[Zone[i]];
    DO_true[i] = DO_true_zone[Zone[i]];
    Temp_true_2[i] = Temp_true[i]^2;
  }
}
model {
  // Define linear predictors
  vector[N] mu_Goby;
  vector[N] Breach_nu;
  vector[J] Temp_zone_nu;
  vector[J] DO_zone_nu;
  
  // --- Priors ---
  phi ~ exponential(1);
  tau ~ exponential(1);
  k_U_phys ~ exponential(1);
  
  {a_DO, a_Temp, a_BreachDays} ~ normal(0, 0.5);
  
  // *** MEASUREMENT ERROR PRIORS ***
  sigma_Temp_measure ~ exponential(0.5);
  sigma_DO_measure ~ exponential(0.5);
  sigma_Temp_within ~ exponential(1);
  sigma_DO_within ~ exponential(1);
  
  // *** STRONG INFORMATIVE PRIORS ***
  beta_Temp_DO ~ normal(-0.4, 0.15);    // Physics: warm water holds less O2
  beta_DO_Goby ~ normal(0.6, 0.2);      // Biology: fish need oxygen
  beta_Temp_Goby ~ normal(0.3, 0.3);    // Biology: temperature effects
  
  // Other priors
  beta_Wind_Temp ~ normal(0, 0.5);
  beta_Wind_DO ~ normal(0, 0.5);
  beta_Substrate_Goby ~ normal(0, 0.5);
  beta_Micro ~ normal(0.25, 0.25);
  beta_Year ~ normal(0, 0.5);
  beta_Year_2 ~ normal(0, 0.5);
  beta_Goby_lag ~ normal(0.25, 0.25);
  beta_BreachDays_2 ~ normal(-0.10, 0.25);
  beta_Temp_2 ~ normal(-0.1, 0.25);
  beta_SAV_2 ~ normal(-0.1, 0.25);
  beta_SAV_Goby ~ normal(0, 0.5);
  beta_BreachDays_Goby ~ normal(0, 0.5);
  beta_Rain_Breach ~ normal(0.5, 0.3);
  beta_Breach_Temp ~ normal(0, 0.5);
  beta_SC_count ~ normal(0, 0.5);
  beta_SB_count ~ normal(0, 0.5);
  
  // Weather correlation
  a_weather ~ normal(0, 0.5);
  sigma_weather ~ exponential(1);
  Rho_weather ~ lkj_corr(2);
  
  // Latent variable
  U_phys ~ std_normal();
  
  // Goby varying intercepts
  a_bar_Goby ~ normal(0, 0.5);
  sigma_Goby ~ exponential(1);
  a_Goby ~ normal(a_bar_Goby, sigma_Goby);
  
  // --- Likelihoods ---
  
  // *** ZONE-LEVEL TRUE ENVIRONMENTAL PROCESSES ***
  for (i in 1:N) {
    Breach_nu[i] = a_BreachDays + 
                   beta_Rain_Breach * Rain[i] + 
                   k_U_phys * U_phys[Zone[i]];
  }
  
  // Zone-level environmental conditions
  for (j in 1:J) {
    // Average BreachDays per zone
    vector[N] breach_in_zone;
    int n_in_zone = 0;
    real avg_breach = 0;
    real avg_wind = 0;
    
    for (i in 1:N) {
      if (Zone[i] == j) {
        n_in_zone += 1;
        avg_breach += BreachDays[i];
        avg_wind += Wind[i];
      }
    }
    
    if (n_in_zone > 0) {
      avg_breach /= n_in_zone;
      avg_wind /= n_in_zone;
    }
    
    // TRUE zone-level temperature
    Temp_zone_nu[j] = a_Temp + 
                      beta_Breach_Temp * avg_breach + 
                      beta_Wind_Temp * avg_wind + 
                      k_U_phys * U_phys[j];
    
    // TRUE zone-level DO
    DO_zone_nu[j] = a_DO + 
                    beta_Temp_DO * Temp_true_zone[j] + 
                    beta_Wind_DO * avg_wind + 
                    k_U_phys * U_phys[j];
  }
  
  // Process likelihoods
  BreachDays ~ normal(Breach_nu, tau);
  Temp_true_zone ~ normal(Temp_zone_nu, tau);
  DO_true_zone ~ normal(DO_zone_nu, tau);
  
  // *** MEASUREMENT MODEL ***
  // Measured values deviate from zone mean due to:
  // 1. True within-zone spatial/temporal variation
  // 2. Measurement error from limited sampling
  for (i in 1:N) {
    Temp_measured[i] ~ normal(Temp_true_zone[Zone[i]], 
                              sqrt(square(sigma_Temp_within) + square(sigma_Temp_measure)));
    DO_measured[i] ~ normal(DO_true_zone[Zone[i]], 
                            sqrt(square(sigma_DO_within) + square(sigma_DO_measure)));
  }
  
  // SC and SB are measured (not modeled)
  
  // Wind/Rain likelihood
  {
    array[N] vector[2] YY;
    for (j in 1:N) YY[j] = [Wind[j], Rain[j]]';
    YY ~ multi_normal(a_weather, quad_form_diag(Rho_weather, sigma_weather));
  }
  
  // *** GOBY RESPONDS TO TRUE ZONE ENVIRONMENT ***
  for (i in 1:N) {
    mu_Goby[i] = a_Goby[Zone[i]] + 
                 beta_SAV_Goby * SAV[i] + 
                 beta_SAV_2 * SAV_2[i] +
                 beta_DO_Goby * DO_true[i] +
                 beta_BreachDays_Goby * BreachDays[i] + 
                 beta_BreachDays_2 * BreachDays_2[i] +
                 beta_Temp_Goby * Temp_true[i] +
                 beta_Temp_2 * Temp_true_2[i] +
                 beta_SC_count * SC_count[i] +
                 beta_SB_count * SB_count[i] +
                 beta_Substrate_Goby * Substrate[i] + 
                 beta_Micro * Micro[i] + 
                 beta_Year * Year[i] + 
                 beta_Year_2 * Year_2[i] + 
                 beta_Goby_lag * Goby_lag[i] + 
                 Area[i];
  }

  Goby ~ neg_binomial_2_log(mu_Goby, phi);
}
generated quantities {
  vector[N] log_lik;
  array[N] int Goby_pred;
  
  // Measurement error estimates
  real measurement_error_Temp = sigma_Temp_measure;
  real measurement_error_DO = sigma_DO_measure;
  real within_zone_var_Temp = sigma_Temp_within;
  real within_zone_var_DO = sigma_DO_within;
  
  for (n in 1:N) {
    real mu_n = a_Goby[Zone[n]] + 
                beta_SAV_Goby * SAV[n] + 
                beta_SAV_2 * SAV_2[n] +
                beta_DO_Goby * DO_true[n] + 
                beta_BreachDays_Goby * BreachDays[n] + 
                beta_BreachDays_2 * BreachDays_2[n] +
                beta_Temp_Goby * Temp_true[n] + 
                beta_Temp_2 * Temp_true_2[n] +
                beta_SC_count * SC_count[n] +
                beta_SB_count * SB_count[n] +
                beta_Substrate_Goby * Substrate[n] + 
                beta_Micro * Micro[n] + 
                beta_Year * Year[n] + 
                beta_Year_2 * Year_2[n] + 
                beta_Goby_lag * Goby_lag[n] + 
                Area[n];

    log_lik[n] = neg_binomial_2_log_lpmf(Goby[n] | mu_n, phi);
    Goby_pred[n] = neg_binomial_2_log_rng(mu_n, phi);
  }
}
"

# Save the model
cat(stan_SEM_code, file = "zone_latent_environment_SEM.stan")

f <- write_stan_file(stan_SEM_code)

# compile model 
mod.SEM <- cmdstan_model(f)   # about 3 minutes
# look
mod.SEM$print()

# Save the model
cat(zone_latent_environment_SEM_code, file = "zone_latent_environment_SEM.stan")



# ============================================================================
# PARAMETER COUNT - FEASIBLE VERSION
# ============================================================================

cat("
=== STREAMLINED ZONE-LEVEL LATENT ENVIRONMENT MODEL ===

PARAMETER COUNT:

ZONE-LEVEL LATENT ENVIRONMENT:
- Temp_true_zone: 3 (one per zone)
- DO_true_zone: 3 (one per zone)
- Within-zone variation SDs: 2 
- Measurement error SDs: 2
Subtotal: 10 parameters

GOBY MODEL: 15 parameters (removed SC/SB pathways)
SUBMODELS: 8 parameters (BreachDays, Temp, DO only)
LATENT CONFOUNDING: 4 parameters (U_phys only)
WEATHER CORRELATION: 5 parameters
ERROR: 2 parameters (tau, phi)

TOTAL: ~44 parameters
Obs per parameter: 314/44 = 7.1

COMPARISON:
┌──────────────────────────────────────────────────────────┐
│ Original model:              ~55-60 params,  5.2-5.7 o/p │
│ Moderate simplified:         ~53 params,     5.9 o/p     │
│ With SC/SB latent env:       ~64 params,     4.9 o/p     │
│ STREAMLINED (no SC/SB):      ~44 params,     7.1 o/p     │
│ Maximal simple (broken):     ~43 params,     7.3 o/p     │
└──────────────────────────────────────────────────────────┘

WHAT WAS REMOVED:
✗ SC and SB submodels (also poorly measured)
✗ U_bio (only needed for biological submodels)
✗ ~20 parameters related to SC/SB pathways

WHAT WAS KEPT (Essential):
✓ Zone-level latent Temp/DO (measurement error correction)
✓ Strong informative priors (prevent attenuation)
✓ U_phys (physical confounding)
✓ Full environmental cascade: Rain → Breach → Temp → DO → Goby
✓ SAV → Goby (SAV measured, direct effect)
✓ SC/SB → Goby (measured, direct effects only)
✓ Wind/Rain correlation
✓ All quadratics

KEY PATHWAYS:
Rain → BreachDays → Temp_true → DO_true → Goby
                              ↓
                            SAV → Goby
                            
SC, SB as measured predictors (no submodels)

ADVANTAGES:
✓ Addresses Temp/DO measurement error problem
✓ Best obs/param ratio while keeping causal structure
✓ Strong priors prevent attenuation
✓ Simpler than original, more correct than maximal simplification
✓ Focus on physical pathway (your main interest)

This should give you the best balance of:
- Statistical power (7.1 obs/param)
- Causal validity (environmental cascade intact)
- Measurement error correction (zone-level latents)
- Interpretability (clear Rain → Breaching story)
")



# ============================================================================
# PARAMETER COUNT - FEASIBLE VERSION
# ============================================================================

cat("
=== STREAMLINED ZONE-LEVEL LATENT ENVIRONMENT MODEL ===

PARAMETER COUNT:

ZONE-LEVEL LATENT ENVIRONMENT:
- Temp_true_zone: 3 (one per zone)
- DO_true_zone: 3 (one per zone)
- Within-zone variation SDs: 2 
- Measurement error SDs: 2
Subtotal: 10 parameters

GOBY MODEL: 15 parameters (removed SC/SB pathways)
SUBMODELS: 8 parameters (BreachDays, Temp, DO only)
LATENT CONFOUNDING: 4 parameters (U_phys only)
WEATHER CORRELATION: 5 parameters
ERROR: 2 parameters (tau, phi)

TOTAL: ~44 parameters
Obs per parameter: 314/44 = 7.1

COMPARISON:
┌──────────────────────────────────────────────────────────┐
│ Original model:              ~55-60 params,  5.2-5.7 o/p │
│ Moderate simplified:         ~53 params,     5.9 o/p     │
│ With SC/SB latent env:       ~64 params,     4.9 o/p     │
│ STREAMLINED (no SC/SB):      ~44 params,     7.1 o/p     │
│ Maximal simple (broken):     ~43 params,     7.3 o/p     │
└──────────────────────────────────────────────────────────┘

WHAT WAS REMOVED:
✗ SC and SB submodels (also poorly measured)
✗ U_bio (only needed for biological submodels)
✗ ~20 parameters related to SC/SB pathways

WHAT WAS KEPT (Essential):
✓ Zone-level latent Temp/DO (measurement error correction)
✓ Strong informative priors (prevent attenuation)
✓ U_phys (physical confounding)
✓ Full environmental cascade: Rain → Breach → Temp → DO → Goby
✓ SAV → Goby (SAV measured, direct effect)
✓ SC/SB → Goby (measured, direct effects only)
✓ Wind/Rain correlation
✓ All quadratics

KEY PATHWAYS:
Rain → BreachDays → Temp_true → DO_true → Goby
                              ↓
                            SAV → Goby
                            
SC, SB as measured predictors (no submodels)

ADVANTAGES:
✓ Addresses Temp/DO measurement error problem
✓ Best obs/param ratio while keeping causal structure
✓ Strong priors prevent attenuation
✓ Simpler than original, more correct than maximal simplification
✓ Focus on physical pathway (your main interest)

This should give you the best balance of:
- Statistical power (7.1 obs/param)
- Causal validity (environmental cascade intact)
- Measurement error correction (zone-level latents)
- Interpretability (clear Rain → Breaching story)
")


#old below-----original full model
