goby_model.stan <- "// Corrected and vectorized Stan model for the Goby Bayesian SEM
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> Zone;

  // Responses and predictors
  array[N] int<lower=0> Goby;
  vector[N] BreachDays;
  vector[N] DO;
  vector[N] Temp;
  vector[N] SAV;
  array[N] int<lower=0, upper=1> SB_count;
  array[N] int<lower=0, upper=1> SC_count;
  vector[N] Wind;
  vector[N] Rain;
  vector[N] Area; // Must be log-transformed in the data list
  array[N] int<lower=0, upper=1> Substrate;
  vector[N] Micro;
  vector[N] Year;
  vector[N] Year_2;
  vector[N] Temp_2;
  vector[N] Goby_lag;
  vector[N] BreachDays_2;
}

parameters {
  // Hyperpriors for Goby varying effects
  real a_bar_Goby;
  real b_bar_Breach;
  vector<lower=0>[2] sigma_Zone;
  cholesky_factor_corr[2] L_Rho_Zone;
  matrix[2, J] z_ab_Goby; // Non-centered raw effects

  // Hyperpriors for shared coefficients
  vector[5] beta_raw;
  real<lower=0> tau_beta;
  matrix[5, J] z_beta_group; // Non-centered raw effects

  // Latent variable
  vector[J] U_raw; // Non-centered

  // Submodel intercepts
  real a_Wind;
  real a_Rain;
  real a_BreachDays;
  real a_DO;
  real a_Temp;
  real a_SB;
  real a_SC;
  real a_SAV;

  // Error terms
  real<lower=0> tau;
  real<lower=0> phi;

  // Fixed coefficients
  real beta_Temp_2;
  real beta_Micro;
  real beta_Substrate;
  real beta_Goby_lag;
  real beta_BreachDays_2;
  real<lower=0> k_U;
  real<lower=0> k_U_bio;

  // Correlated residuals for Wind & Rain
  vector[2] mu_windrain;
  vector<lower=0>[2] sigma_windrain;
  cholesky_factor_corr[2] L_Rho_windrain;
}

transformed parameters {
  // Non-centered parameterization for all varying effects
  vector[J] a_Goby;
  vector[J] b_BreachDays;
  array[J] vector[5] beta_group;
  vector[J] U = U_raw; // Centered U

  // Reconstruct Goby varying intercepts and slopes
  {
    matrix[J, 2] ab = (diag_pre_multiply(sigma_Zone, L_Rho_Zone) * z_ab_Goby)';
    a_Goby = ab[, 1] + a_bar_Goby;
    b_BreachDays = ab[, 2] + b_bar_Breach;
  }

  // Reconstruct hierarchical coefficients for each zone
  for (j in 1:J) {
    beta_group[j] = beta_raw + z_beta_group[, j] * tau_beta;
  }
}

model {
  // --- Priors ---
  a_bar_Goby ~ normal(0, 0.5);
  b_bar_Breach ~ normal(0, 0.5);
  L_Rho_Zone ~ lkj_corr_cholesky(2);
  sigma_Zone ~ exponential(1);
  to_vector(z_ab_Goby) ~ std_normal();

  U_raw ~ std_normal();

  beta_raw ~ normal(0, 0.5);
  tau_beta ~ exponential(1);
  to_vector(z_beta_group) ~ std_normal();

  // Priors for submodel intercepts
  {a_Wind, a_Rain, a_BreachDays, a_DO, a_Temp, a_SB, a_SC, a_SAV} ~ normal(0, 0.5);
  
  // Priors for fixed coefficients
  {beta_Temp_2, beta_Micro, beta_Substrate, beta_Goby_lag, beta_BreachDays_2} ~ normal(0, 0.5);
  {k_U, k_U_bio} ~ exponential(1);

  // Priors for correlated residuals
  mu_windrain ~ normal(0, 0.5);
  sigma_windrain ~ exponential(1);
  L_Rho_windrain ~ lkj_corr_cholesky(2);
  
  // Error terms
  tau ~ exponential(1);
  phi ~ normal(1, 5);

  // --- Likelihoods (Vectorized) ---
  // Joint model for Wind and Rain
  {
    array[N] vector[2] Y_windrain;
    for (n in 1:N) Y_windrain[n] = [Wind[n], Rain[n]]';
    Y_windrain ~ multi_normal_cholesky(mu_windrain, diag_pre_multiply(sigma_windrain, L_Rho_windrain));
  }

  // Submodels
  vector[N] breach_nu;
  vector[N] do_nu;
  vector[N] temp_nu;
  vector[N] sav_nu;
  vector[N] sb_mu;
  vector[N] sc_mu;

  for (n in 1:N) {
    breach_nu[n] = a_BreachDays + beta_group[Zone[n]][1] * Rain[n] + k_U * U[Zone[n]];
    temp_nu[n]   = a_Temp + beta_group[Zone[n]][4] * BreachDays[n] + beta_group[Zone[n]][5] * Wind[n] + k_U * U[Zone[n]];
    do_nu[n]     = a_DO + beta_group[Zone[n]][2] * Temp[n] + beta_group[Zone[n]][3] * Wind[n] + k_U * U[Zone[n]];
    sav_nu[n]    = a_SAV + beta_group[Zone[n]][3] * DO[n] + beta_group[Zone[n]][2] * Temp[n] + k_U_bio * U[Zone[n]];
    sb_mu[n]     = a_SB + beta_group[Zone[n]][1] * DO[n] + beta_group[Zone[n]][2] * SAV[n] + k_U_bio * U[Zone[n]];
    sc_mu[n]     = a_SC + beta_group[Zone[n]][3] * Substrate[n] + beta_group[Zone[n]][4] * DO[n] + k_U_bio * U[Zone[n]];
  }
  
  BreachDays ~ normal(breach_nu, tau);
  DO ~ normal(do_nu, tau);
  Temp ~ normal(temp_nu, tau);
  SAV ~ normal(sav_nu, tau);

  SB_count ~ bernoulli_logit(sb_mu);
  SC_count ~ bernoulli_logit(sc_mu);

  // Goby model
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = exp(
      a_Goby[Zone[n]] +
      b_BreachDays[Zone[n]] * BreachDays[n] +
      beta_group[Zone[n]][1] * SAV[n] +
      beta_group[Zone[n]][2] * DO[n] +
      beta_group[Zone[n]][3] * Temp[n] +
      beta_group[Zone[n]][4] * SC_count[n] +
      beta_group[Zone[n]][5] * SB_count[n] +
      beta_Micro * Micro[n] +
      beta_Substrate * Substrate[n] +
      beta_Goby_lag * Goby_lag[n] +
      beta_BreachDays_2 * BreachDays_2[n] +
      beta_Temp_2 * Temp_2[n] +
      Area[n] // Offset term (coefficient is fixed to 1)
    );
  }
  Goby ~ neg_binomial_2(mu, exp(phi));
}
"
