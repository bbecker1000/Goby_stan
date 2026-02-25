library(rethinking)

Goby.Hierarchical_SEM <- ulam(
  alist(
    # --- Goby Model ---
    # Coefficients that were shared are now indexed from their respective vectors
    Goby ~ dgampois(mu, exp(phi)),
    log(mu) <-
      a_Goby[Zone] +
      beta_SAV_vec[1] * SAV +         ## CHANGED
      beta_DO_vec[1] * DO +           ## CHANGED
      beta_BreachDays_vec[1] * BreachDays + ## CHANGED
      beta_Temp_vec[1] * Temp +       ## CHANGED
      beta_SC_count * SC_count +
      beta_Year * Year +
      beta_Year_2 * Year_2 +
      beta_SB_count * SB_count +
      beta_Micro * Micro +
      beta_Substrate_vec[1] * Substrate + ## CHANGED
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
    
    # --- Submodels with Unshared but Correlated Parameters ---
    # Wind and Rain still have correlated errors (Unchanged)
    c(Wind, Rain) ~ multi_normal( c(wind_mu, rain_mu), Rho_errors, sigma_errors ),
    wind_mu <- a_Wind,
    rain_mu <- a_Rain,
    Rho_errors ~ dlkjcorr(2),
    sigma_errors ~ dexp(1),
    
    # Physical submodels
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <- a_BreachDays + beta_Rain * Rain + k_U * U[Zone],
    
    DO ~ normal(DO_nu, tau),
    DO_nu <- a_DO + beta_Temp_vec[2] * Temp + beta_Wind_vec[1] * Wind + k_U * U[Zone], ## CHANGED
    
    Temp ~ normal(Temp_nu, tau),
    Temp_nu <- a_Temp + beta_BreachDays_vec[2] * BreachDays + beta_Wind_vec[2] * Wind + k_U * U[Zone], ## CHANGED
    
    # Biological submodels
    SB_count ~ dbinom(1, SB_mu),
    logit(SB_mu) <- a_SB + beta_DO_vec[2] * DO + beta_SAV_vec[2] * SAV + k_U_bio * U[Zone], ## CHANGED
    
    SC_count ~ dbinom(1, SC_mu),
    logit(SC_mu) <- a_SC + beta_Substrate_vec[2] * Substrate + beta_DO_vec[3] * DO + beta_SAV_vec[3] * SAV + k_U_bio * U[Zone], ## CHANGED
    
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- a_SAV + beta_DO_vec[4] * DO + beta_Temp_vec[3] * Temp + k_U_bio * U[Zone], ## CHANGED
    
    # --- Priors ---
    ## ADDED: Hierarchical priors for the now-unshared coefficients
    # Each vector of related betas gets its own multivariate normal prior
    
    # 1. DO effects (4 parameters)
    vector[4]:beta_DO_vec ~ multi_normal(b_bar_DO, Rho_DO, sigma_DO),
    vector[4]:b_bar_DO ~ normal(0, 0.5),
    Rho_DO ~ dlkjcorr(2),
    sigma_DO ~ dexp(1),
    
    # 2. SAV effects (3 parameters)
    vector[3]:beta_SAV_vec ~ multi_normal(b_bar_SAV, Rho_SAV, sigma_SAV),
    vector[3]:b_bar_SAV ~ normal(0, 0.5),
    Rho_SAV ~ dlkjcorr(2),
    sigma_SAV ~ dexp(1),
    
    # 3. Temp effects (3 parameters)
    vector[3]:beta_Temp_vec ~ multi_normal(b_bar_Temp, Rho_Temp, sigma_Temp),
    vector[3]:b_bar_Temp ~ normal(0, 0.5),
    Rho_Temp ~ dlkjcorr(2),
    sigma_Temp ~ dexp(1),
    
    # 4. BreachDays effects (2 parameters)
    vector[2]:beta_BreachDays_vec ~ multi_normal(b_bar_Breach, Rho_Breach, sigma_Breach),
    vector[2]:b_bar_Breach ~ normal(0, 0.5),
    Rho_Breach ~ dlkjcorr(2),
    sigma_Breach ~ dexp(1),
    
    # 5. Wind effects (2 parameters)
    vector[2]:beta_Wind_vec ~ multi_normal(b_bar_Wind, Rho_Wind, sigma_Wind),
    vector[2]:b_bar_Wind ~ normal(0, 0.5),
    Rho_Wind ~ dlkjcorr(2),
    sigma_Wind ~ dexp(1),
    
    # 6. Substrate effects (2 parameters)
    vector[2]:beta_Substrate_vec ~ multi_normal(b_bar_Substrate, Rho_Substrate, sigma_Substrate),
    vector[2]:b_bar_Substrate ~ normal(0, 0.5),
    Rho_Substrate ~ dlkjcorr(2),
    sigma_Substrate ~ dexp(1),
    
    ## Priors for remaining fixed effects
    c(
      a_DO, a_SB, a_SC, a_SAV, a_Temp, a_BreachDays, a_Wind, a_Rain,
      beta_SC_count, beta_SB_count, beta_Year, beta_Year_2, beta_Rain,
      beta_Temp_2, beta_Micro, beta_Goby_lag, beta_BreachDays_2
    ) ~ normal(0, 0.5),
    
    ## Priors for other parameters (Unchanged)
    c(k_U, k_U_bio) ~ dexp(1),
    tau ~ dexp(1),
    phi ~ dnorm(1, 5)
  ),
  data = sim_data, chains = 1, cores = parallel::detectCores(), iter = 4000, warmup = 3000,
  cmdstan = TRUE
)

precis(Goby.Hierarchical_SEM, depth = 2)

stanmodelcode <- rethinking::stancode(Goby.Hierarchical_SEM)


