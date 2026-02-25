# 4c flbi.goby
# with SAV, Rain, and DO as multivariation prior by zone
# thanks Gemini
# with geneartive submodels and shared coefficients
# remember to chec if using dat or simulated_data

Goby.Simple_DAG.sim <- ulam(
  alist(
    ## --- Goby Model ---
    Goby ~ dgampois(mu, exp(phi)),
    log(mu) <-
      a_Goby[Zone] +
      beta_SAV * SAV +
      beta_DO * DO +
      beta_BreachDays * BreachDays +
      beta_Temp * Temp +
      beta_SC_count * SC_count +
      beta_Year * Year +
      beta_Year_2 * Year_2 +
      beta_SB_count * SB_count +
      beta_Micro * Micro +
      beta_Substrate * Substrate +
      beta_Temp_2 * Temp_2 +
      beta_Goby_lag * Goby_lag +
      beta_BreachDays_2 * BreachDays_2 +
      Area,
    
    ## Varying intercept model
    a_Goby[Zone] ~ dnorm(a_bar_Goby, sigma_Goby),
    a_bar_Goby ~ dnorm(0, 0.5),
    sigma_Goby ~ dexp(1),
    
    ## --- Submodels with Independent Errors ---
    # Wind and Rain are now separate models
    Wind ~ normal(wind_mu, tau_wind),
    wind_mu <- a_Wind,
    
    Rain ~ normal(rain_mu, tau_rain),
    rain_mu <- a_Rain,
    
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <- a_BreachDays + beta_Rain * Rain,
    
    DO ~ normal(DO_nu, tau),
    DO_nu <- a_DO + beta_Temp * Temp + beta_Wind * Wind,
    
    Temp ~ normal(Temp_nu, tau),
    Temp_nu <- a_Temp + beta_BreachDays * BreachDays + beta_Wind * Wind,
    
    SB_count ~ dbinom(1, SB_mu),
    logit(SB_mu) <- a_SB + beta_DO * DO + beta_SAV * SAV,
    
    SC_count ~ dbinom(1, SC_mu),
    logit(SC_mu) <- a_SC + beta_Substrate * Substrate + beta_DO * DO + beta_SAV * SAV,
    
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- a_SAV + beta_DO * DO + beta_Temp * Temp,
    
    ## --- Priors ---
    c(
      a_DO, a_SB, a_SC, a_SAV, a_Temp, a_BreachDays, a_Wind, a_Rain,
      beta_SAV, beta_DO, beta_BreachDays, beta_Temp, beta_SC_count,
      beta_SB_count, beta_Year, beta_Year_2, beta_Wind, beta_Rain
    ) ~ normal(0, 0.5),
    
    beta_Temp_2 ~ normal(-0.10, 0.25),
    beta_Micro ~ normal(0.25, 0.25),
    beta_SAV_2 ~ normal(-0.10, 0.25),
    beta_Substrate ~ normal(0.25, 0.25),
    beta_Goby_lag ~ normal(0.25, 0.25),
    beta_BreachDays_2 ~ dnorm(-0.1, 0.25),
    
    # Priors for all error terms
    c(tau, tau_wind, tau_rain) ~ dexp(1),
    phi ~ dnorm(1, 5)
  ),
  data = simulated_data, chains = 1, cores = parallel::detectCores(), iter = 4000, warmup = 3000,
  cmdstan = TRUE
)

precis(Goby.Simple_DAG.sim, depth = 1)


###

#makes most biological sense - using this one...

Goby.MVnorm_U_shared_bio.sim <- ulam(
  alist(
    # --- Goby Model (Unchanged) ---
    Goby ~ dgampois(mu, exp(phi)),
    log(mu) <-
      a_Goby[Zone] +
      beta_SAV * SAV +
      beta_DO * DO +
      beta_BreachDays * BreachDays +
      beta_Temp * Temp +
      beta_SC_count * SC_count +
      beta_Year * Year +
      beta_Year_2 * Year_2 +
      beta_SB_count * SB_count +
      beta_Micro * Micro +
      beta_Substrate * Substrate +
      beta_Temp_2 * Temp_2 +
      beta_Goby_lag * Goby_lag +
      beta_BreachDays_2 * BreachDays_2 +
      Area,
    
    # Varying intercept model (Unchanged)
    a_Goby[Zone] ~ dnorm(a_bar_Goby, sigma_Goby),
    a_bar_Goby ~ dnorm(0, 0.5),
    sigma_Goby ~ dexp(1),
    
    # Latent Variable Definition (Unchanged)
    U[Zone] ~ dnorm(0, 1),
    
    # --- Submodels with U as a common cause ---
    # Wind and Rain still have correlated errors
    c(Wind, Rain) ~ multi_normal( c(wind_mu, rain_mu), Rho_errors, sigma_errors ),
    wind_mu <- a_Wind,
    rain_mu <- a_Rain,
    Rho_errors ~ dlkjcorr(2),
    sigma_errors ~ dexp(1),
    
    # BreachDays, Temp, and DO influenced by a SHARED U effect
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <- a_BreachDays + beta_Rain * Rain + k_U * U[Zone],
    
    DO ~ normal(DO_nu, tau),
    DO_nu <- a_DO + beta_Temp * Temp + beta_Wind * Wind + k_U * U[Zone],
    
    Temp ~ normal(Temp_nu, tau),
    Temp_nu <- a_Temp + beta_BreachDays * BreachDays + beta_Wind * Wind + k_U * U[Zone],
    
    # --- Other submodels now ALSO influenced by a NEW shared U effect ---
    SB_count ~ dbinom(1, SB_mu),
    logit(SB_mu) <- a_SB + beta_DO * DO + beta_SAV * SAV + k_U_bio * U[Zone], ## CHANGED
    
    SC_count ~ dbinom(1, SC_mu),
    logit(SC_mu) <- a_SC + beta_Substrate * Substrate + beta_DO * DO + beta_SAV * SAV + k_U_bio * U[Zone], ## CHANGED
    
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- a_SAV + beta_DO * DO + beta_Temp * Temp + k_U_bio * U[Zone], ## CHANGED
    
    # --- Priors ---
    c(
      a_DO, a_SB, a_SC, a_SAV, a_Temp, a_BreachDays, a_Wind, a_Rain,
      beta_SAV, beta_DO, beta_BreachDays, beta_Temp, beta_SC_count,
      beta_SB_count, beta_Year, beta_Year_2, beta_Wind, beta_Rain
    ) ~ normal(0, 0.5),
    
    beta_Temp_2 ~ normal(-0.10, 0.25),
    beta_Micro ~ normal(0.25, 0.25),
    beta_SAV_2 ~ normal(-0.10, 0.25),
    beta_Substrate ~ normal(0.25, 0.25),
    beta_Goby_lag ~ normal(0.25, 0.25),
    beta_BreachDays_2 ~ dnorm(-0.1, 0.25),
    
    # Priors for the TWO shared latent variable coefficients
    c(k_U, k_U_bio) ~ dexp(1), ## CHANGED
    
    tau ~ dexp(1),
    phi ~ dnorm(1, 5)
  ),
  data = simulated_data, chains = 1, cores = parallel::detectCores(), iter = 4000, warmup = 3000,
  cmdstan = TRUE
)

precis(Goby.MVnorm_U_shared_bio.sim, depth = 1)

stanmodelcode <- rethinking::stancode(Goby.MVnorm_U_shared_bio.sim)

saveRDS(Goby.MVnorm_U_shared_bio.sim, file = "Output/Models/Goby.MVnorm_U_shared_bio.sim.rds")


# and add breach as Random slope
Goby.MVnorm_U_shared_bio_RS.sim <- ulam(
  alist(
    # --- Goby Model ---
    Goby ~ dgampois(mu, exp(phi)),
    log(mu) <-
      a_Goby[Zone] +                 # Varying intercept (now part of MVN)
      b_BreachDays[Zone] * BreachDays +  ## CHANGED: Slope now varies by Zone
      beta_SAV * SAV +
      beta_DO * DO +
      beta_Temp * Temp +
      beta_SC_count * SC_count +
      beta_Year * Year +
      beta_Year_2 * Year_2 +
      beta_SB_count * SB_count +
      beta_Micro * Micro +
      beta_Substrate * Substrate +
      beta_Temp_2 * Temp_2 +
      beta_Goby_lag * Goby_lag +
      beta_BreachDays_2 * BreachDays_2 +
      Area,
    
    ## CHANGED: Upgraded to a multivariate prior for intercepts and slopes
    c(a_Goby, b_BreachDays)[Zone] ~ multi_normal( c(a_bar_Goby, b_bar_Breach), Rho_Zone, sigma_Zone ),
    a_bar_Goby ~ dnorm(0, 0.5),
    b_bar_Breach ~ dnorm(0, 0.5), # Hyperprior for the average BreachDays slope
    Rho_Zone ~ dlkjcorr(2),      # Prior for the correlation matrix
    sigma_Zone ~ dexp(1),        # Prior for the standard deviations
    
    # Latent Variable Definition (Unchanged)
    U[Zone] ~ dnorm(0, 1),
    
    # --- Submodels ---
    # Wind and Rain still have correlated errors
    c(Wind, Rain) ~ multi_normal( c(wind_mu, rain_mu), Rho_errors, sigma_errors ),
    wind_mu <- a_Wind,
    rain_mu <- a_Rain,
    Rho_errors ~ dlkjcorr(2),
    sigma_errors ~ dexp(1),
    
    # Physical submodels influenced by shared U effect
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <- a_BreachDays + beta_Rain * Rain + k_U * U[Zone],
    
    DO ~ normal(DO_nu, tau),
    DO_nu <- a_DO + beta_Temp * Temp + beta_Wind * Wind + k_U * U[Zone],
    
    Temp ~ normal(Temp_nu, tau),
    Temp_nu <- a_Temp + beta_BreachDays_sub * BreachDays + beta_Wind * Wind + k_U * U[Zone], ## CHANGED
    
    # Biological submodels influenced by a second shared U effect
    SB_count ~ dbinom(1, SB_mu),
    logit(SB_mu) <- a_SB + beta_DO * DO + beta_SAV * SAV + k_U_bio * U[Zone],
    
    SC_count ~ dbinom(1, SC_mu),
    logit(SC_mu) <- a_SC + beta_Substrate * Substrate + beta_DO * DO + beta_SAV * SAV + k_U_bio * U[Zone],
    
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- a_SAV + beta_DO * DO + beta_Temp * Temp + k_U_bio * U[Zone],
    
    # --- Priors ---
    c(
      a_DO, a_SB, a_SC, a_SAV, a_Temp, a_BreachDays, a_Wind, a_Rain,
      beta_SAV, beta_DO, beta_Temp, beta_SC_count,
      beta_SB_count, beta_Year, beta_Year_2, beta_Wind, beta_Rain,
      beta_BreachDays_sub ## ADDED: Renamed parameter for Temp submodel
    ) ~ normal(0, 0.5),
    
    beta_Temp_2 ~ normal(-0.10, 0.25),
    beta_Micro ~ normal(0.25, 0.25),
    beta_SAV_2 ~ normal(-0.10, 0.25),
    beta_Substrate ~ normal(0.25, 0.25),
    beta_Goby_lag ~ normal(0.25, 0.25),
    beta_BreachDays_2 ~ dnorm(-0.1, 0.25),
    
    # Priors for the two shared latent variable coefficients
    c(k_U, k_U_bio) ~ dexp(1),
    
    tau ~ dexp(1),
    phi ~ dnorm(1, 5)
  ),
  data = simulated_data, chains = 1, cores = parallel::detectCores(), iter = 4000, warmup = 3000,
  cmdstan = TRUE
)

precis(Goby.MVnorm_U_shared_bio_RS.sim, depth = 1)

stanmodelcode <- rethinking::stancode(Goby.MVnorm_U_shared_bio_RS.sim)

saveRDS(Goby.MVnorm_U_shared_bio.sim, file = "Output/Models/Goby.MVnorm_U_shared_bio_RS.sim.rds")



