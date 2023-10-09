

stan_program <- "

data{
  vector[299] Breach;
  vector[299] SC;
  vector[299] SB;
  array[299] int Year_int;
  array[299] int Goby;
  array[299] int Substrate;
  vector[299] Year;
  array[299] int Zone;
  array[299] int SB_count;
  array[299] int SC_count;
  vector[299] DO;
  array[299] int Breach_count;
  vector[299] Area;
  vector[299] Wind;
  vector[299] SAV;
  vector[299] Temp;
  vector[299] Rain;
}
parameters{
  real mu_Zone;
  real<lower=0> tau_Zone;
  real beta_Substrate;
  real beta_Wind;
  real beta_Breach_count;
  real beta_Temp;
  real beta_DO;
  real beta_SB_count;
  real beta_SAV;
  real beta_SC_count;
  real beta_Rain;
  real beta_Year;
  real a_Temp;
  real a_SAV;
  real a_SC;
  real a_SB;
  real a_DO;
  real a_Breach;
  vector[3] a_Goby;
  real<lower=0> tau;
  real phi;
  real SB_phi;
  real SC_phi;
}
model{
  vector[299] mu;
  vector[299] DO_nu;
  vector[299] Temp_nu;
  vector[299] SB_mu;
  vector[299] SC_mu;
  vector[299] Breach_nu;
  vector[299] SAV_nu;
  SC_phi ~ normal( 1 , 10 );
  SB_phi ~ normal( 1 , 10 );
  phi ~ normal( 1 , 10 );
  tau ~ exponential( 1 );
  a_Goby ~ normal( 0 , 0.5 );
  a_Breach ~ normal( 0 , 0.5 );
  a_DO ~ normal( 0 , 0.5 );
  a_SB ~ normal( 0 , 0.5 );
  a_SC ~ normal( 0 , 0.5 );
  a_SAV ~ normal( 0 , 0.5 );
  a_Temp ~ normal( 0 , 0.5 );
  beta_Year ~ normal( 0 , 0.5 );
  beta_Rain ~ normal( 0 , 0.5 );
  beta_SC_count ~ normal( 0 , 0.5 );
  beta_SAV ~ normal( 0 , 0.5 );
  beta_SB_count ~ normal( 0 , 0.5 );
  beta_DO ~ normal( 0 , 0.5 );
  beta_Temp ~ normal( 0 , 0.5 );
  beta_Breach_count ~ normal( 0 , 0.5 );
  beta_Wind ~ normal( 0 , 0.5 );
  beta_Substrate ~ normal( 0 , 0.5 );
  for ( i in 1:299 ) {
    SAV_nu[i] = a_SAV + beta_Rain * Rain[i] + beta_Temp * Temp[i];
  }
  SAV ~ normal( SAV_nu , tau );
  for ( i in 1:299 ) {
    Breach_nu[i] = a_Breach + beta_Rain * Rain[i] + beta_Wind * Wind[i] + Area[i];
    Breach_nu[i] = exp(Breach_nu[i]);
  }
  Breach_count ~ poisson( Breach_nu );
  for ( i in 1:299 ) {
    SC_mu[i] = a_SC + beta_DO * DO[i] + beta_SAV * SAV[i] + Area[i];
    SC_mu[i] = exp(SC_mu[i]);
  }
  SC_count ~ neg_binomial_2( SC_mu , exp(SC_phi) );
  for ( i in 1:299 ) {
    SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i] + Area[i];
    SB_mu[i] = exp(SB_mu[i]);
  }
  SB_count ~ neg_binomial_2( SB_mu , exp(SB_phi) );
  for ( i in 1:299 ) {
    Temp_nu[i] = a_Temp + beta_Breach_count * Breach_count[i] + beta_Wind * Wind[i];
  }
  Temp ~ normal( Temp_nu , tau );
  for ( i in 1:299 ) {
    DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_Breach_count * Breach_count[i] + beta_Wind * Wind[i] + beta_SAV * SAV[i];
  }
  DO ~ normal( DO_nu , tau );
  tau_Zone ~ exponential( 1 );
  mu_Zone ~ normal( 0 , 5 );
  a_Goby ~ normal( mu_Zone , tau_Zone );
  for ( i in 1:299 ) {
    mu[i] = a_Goby[Zone[i]] + beta_Year * Year[i] + beta_SC_count * SC_count[i] + beta_SAV * SAV[i] + beta_SB_count * SB_count[i] + beta_DO * DO[i] + beta_Substrate * Substrate[i] + Area[i];
    mu[i] = exp(mu[i]);
  }
  Goby ~ neg_binomial_2( mu , exp(phi) );
}
"


library(rstan)
fit1 <- stan(
  model_code = stan_program,  # Stan program
  data = dat.list,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 250,          # number of warmup iterations per chain
  iter = 500,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = 1            # no progress shown
)











