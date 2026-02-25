#stan for sims
# has priors similar to the simulations

stan_program_sim <- "
data{
    array[314] int Breach;
    array[314] int Year_int;
    array[314] int SB;
    array[314] int SC;
    array[314] int Breach_count;
    array[314] int BreachDays_Count_2;
    array[314] int BreachDays_Count;
     vector[314] SAV_2;
    array[314] int Goby;
     vector[314] Area;
     vector[314] BreachDays_2;
     vector[314] Goby_lag;
     vector[314] Temp_2;
     vector[314] Micro;
     vector[314] Year_2;
     vector[314] Year;
     vector[314] Rain;
     vector[314] Wind;
     vector[314] BreachDays;
    array[314] int SB_count;
    array[314] int SC_count;
    array[314] int Substrate;
     vector[314] SAV;
    array[314] int Zone;
     vector[314] Temp;
     vector[314] DO;
}
parameters{
     vector[3] a_Goby;
     real a_bar_Goby;
     real<lower=0> sigma_Goby;
     vector[3] U;
     corr_matrix[2] Rho_errors;
     vector<lower=0>[2] sigma_errors;
     real beta_Rain;
     real beta_Wind;
     real beta_Year_2;
     real beta_Year;
     real beta_SB_count;
     real beta_SC_count;
     real beta_Temp;
     real beta_BreachDays;
     real beta_DO;
     real beta_SAV;
     real a_Rain;
     real a_Wind;
     real a_BreachDays;
     real a_Temp;
     real a_SAV;
     real a_SC;
     real a_SB;
     real a_DO;
     real beta_Temp_2;
     real beta_Micro;
     real beta_SAV_2;
     real beta_Substrate;
     real beta_Goby_lag;
     real beta_BreachDays_2;
     real<lower=0> k_U_bio;
     real<lower=0> k_U;
     real<lower=0> tau;
     real phi;
}
model{
     vector[314] mu;
     real wind_mu;
     real rain_mu;
     vector[314] Breach_nu;
     vector[314] DO_nu;
     vector[314] Temp_nu;
     vector[314] SB_mu;
     vector[314] SC_mu;
     vector[314] SAV_nu;
    phi ~ normal( 1 , 5 );
    tau ~ exponential( 1 );
    k_U ~ exponential( 1 );
    k_U_bio ~ exponential( 1 );
    beta_BreachDays_2 ~ normal( -0.1 , 0.25 );
    beta_Goby_lag ~ normal( 0.25 , 0.25 );
    beta_Substrate ~ normal( 0.25 , 0.25 );
    beta_SAV_2 ~ normal( -0.1 , 0.25 );
    beta_Micro ~ normal( 0.25 , 0.25 );
    beta_Temp_2 ~ normal( -0.1 , 0.25 );
    a_DO ~ normal( 0 , 0.5 );
    a_SB ~ normal( 0 , 0.5 );
    a_SC ~ normal( 0 , 0.5 );
    a_SAV ~ normal( 0 , 0.5 );
    a_Temp ~ normal( 0 , 0.5 );
    a_BreachDays ~ normal( 0 , 0.5 );
    a_Wind ~ normal( 0 , 0.5 );
    a_Rain ~ normal( 0 , 0.5 );
    beta_SAV ~ normal( 0 , 0.5 );
    beta_DO ~ normal( 0 , 0.5 );
    beta_BreachDays ~ normal( 0 , 0.5 );
    beta_Temp ~ normal( 0 , 0.5 );
    beta_SC_count ~ normal( 0 , 0.5 );
    beta_SB_count ~ normal( 0 , 0.5 );
    beta_Year ~ normal( 0 , 0.5 );
    beta_Year_2 ~ normal( 0 , 0.5 );
    beta_Wind ~ normal( 0 , 0.5 );
    beta_Rain ~ normal( 0 , 0.5 );
    for ( i in 1:314 ) {
        SAV_nu[i] = a_SAV + beta_DO * DO[i] + beta_Temp * Temp[i] + k_U_bio * U[Zone[i]];
    }
    SAV ~ normal( SAV_nu , tau );
    for ( i in 1:314 ) {
        SC_mu[i] = a_SC + beta_Substrate * Substrate[i] + beta_DO * DO[i] + beta_SAV * SAV[i] + k_U_bio * U[Zone[i]];
        SC_mu[i] = inv_logit(SC_mu[i]);
    }
    SC_count ~ binomial( 1 , SC_mu );
    for ( i in 1:314 ) {
        SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i] + k_U_bio * U[Zone[i]];
        SB_mu[i] = inv_logit(SB_mu[i]);
    }
    SB_count ~ binomial( 1 , SB_mu );
    for ( i in 1:314 ) {
        Temp_nu[i] = a_Temp + beta_BreachDays * BreachDays[i] + beta_Wind * Wind[i] + k_U * U[Zone[i]];
    }
    Temp ~ normal( Temp_nu , tau );
    for ( i in 1:314 ) {
        DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_Wind * Wind[i] + k_U * U[Zone[i]];
    }
    DO ~ normal( DO_nu , tau );
    for ( i in 1:314 ) {
        Breach_nu[i] = a_BreachDays + beta_Rain * Rain[i] + k_U * U[Zone[i]];
    }
    BreachDays ~ normal( Breach_nu , tau );
    sigma_errors ~ exponential( 1 );
    Rho_errors ~ lkj_corr( 2 );
    rain_mu = a_Rain;
    wind_mu = a_Wind;
    {
    array[314] vector[2] YY;
    vector[2] MU;
    MU = [ wind_mu , rain_mu ]';
    for ( j in 1:314 ) YY[j] = [ Wind[j] , Rain[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho_errors , sigma_errors) );
    }
    U ~ normal( 0 , 1 );
    sigma_Goby ~ exponential( 1 );
    a_bar_Goby ~ normal( 0 , 0.5 );
    a_Goby ~ normal( a_bar_Goby , sigma_Goby );
    for ( i in 1:314 ) {
        mu[i] = a_Goby[Zone[i]] + beta_SAV * SAV[i] + beta_DO * DO[i] + beta_BreachDays * BreachDays[i] + beta_Temp * Temp[i] + beta_SC_count * SC_count[i] + beta_Year * Year[i] + beta_Year_2 * Year_2[i] + beta_SB_count * SB_count[i] + beta_Micro * Micro[i] + beta_Substrate * Substrate[i] + beta_Temp_2 * Temp_2[i] + beta_Goby_lag * Goby_lag[i] + beta_BreachDays_2 * BreachDays_2[i] + Area[i];
        mu[i] = exp(mu[i]);
    }
    Goby ~ neg_binomial_2( mu , exp(phi) );
}
"

# write to stan file
f.stan_program_sim <- write_stan_file(stan_program_sim)

# compile model 
# about 3 minutes
mod.sim <- cmdstan_model(f.stan_program_sim)


## or for RS for Breach


stan_program.RS_Breach.sim <- "
data{
  array[314] int Breach;
  array[314] int Year_int;
  array[314] int SB;
  array[314] int SC;
  array[314] int Breach_count;
  array[314] int BreachDays_Count_2;
  array[314] int BreachDays_Count;
  vector[314] SAV_2;
  array[314] int Goby;
  vector[314] Area;
  vector[314] BreachDays_2;
  vector[314] Goby_lag;
  vector[314] Temp_2;
  vector[314] Micro;
  vector[314] Year_2;
  vector[314] Year;
  vector[314] Rain;
  vector[314] Wind;
  vector[314] BreachDays;
  array[314] int SB_count;
  array[314] int SC_count;
  array[314] int Substrate;
  vector[314] SAV;
  array[314] int Zone;
  vector[314] Temp;
  vector[314] DO;
}
parameters{
  vector[3] b_BreachDays;
  vector[3] a_Goby;
  real a_bar_Goby;
  real b_bar_Breach;
  corr_matrix[2] Rho_Zone;
  vector<lower=0>[2] sigma_Zone;
  vector[3] U;
  corr_matrix[2] Rho_errors;
  vector<lower=0>[2] sigma_errors;
  real beta_BreachDays_sub;
  real beta_Rain;
  real beta_Wind;
  real beta_Year_2;
  real beta_Year;
  real beta_SB_count;
  real beta_SC_count;
  real beta_Temp;
  real beta_DO;
  real beta_SAV;
  real a_Rain;
  real a_Wind;
  real a_BreachDays;
  real a_Temp;
  real a_SAV;
  real a_SC;
  real a_SB;
  real a_DO;
  real beta_Temp_2;
  real beta_Micro;
  real beta_SAV_2;
  real beta_Substrate;
  real beta_Goby_lag;
  real beta_BreachDays_2;
  real<lower=0> k_U_bio;
  real<lower=0> k_U;
  real<lower=0> tau;
  real phi;
}
model{
  vector[314] mu;
  real wind_mu;
  real rain_mu;
  vector[314] Breach_nu;
  vector[314] DO_nu;
  vector[314] Temp_nu;
  vector[314] SB_mu;
  vector[314] SC_mu;
  vector[314] SAV_nu;
  phi ~ normal( 1 , 5 );
  tau ~ exponential( 1 );
  k_U ~ exponential( 1 );
  k_U_bio ~ exponential( 1 );
  beta_BreachDays_2 ~ normal( -0.1 , 0.25 );
  beta_Goby_lag ~ normal( 0.25 , 0.25 );
  beta_Substrate ~ normal( 0.25 , 0.25 );
  beta_SAV_2 ~ normal( -0.1 , 0.25 );
  beta_Micro ~ normal( 0.25 , 0.25 );
  beta_Temp_2 ~ normal( -0.1 , 0.25 );
  a_DO ~ normal( 0 , 0.5 );
  a_SB ~ normal( 0 , 0.5 );
  a_SC ~ normal( 0 , 0.5 );
  a_SAV ~ normal( 0 , 0.5 );
  a_Temp ~ normal( 0 , 0.5 );
  a_BreachDays ~ normal( 0 , 0.5 );
  a_Wind ~ normal( 0 , 0.5 );
  a_Rain ~ normal( 0 , 0.5 );
  beta_SAV ~ normal( 0 , 0.5 );
  beta_DO ~ normal( 0 , 0.5 );
  beta_Temp ~ normal( 0 , 0.5 );
  beta_SC_count ~ normal( 0 , 0.5 );
  beta_SB_count ~ normal( 0 , 0.5 );
  beta_Year ~ normal( 0 , 0.5 );
  beta_Year_2 ~ normal( 0 , 0.5 );
  beta_Wind ~ normal( 0 , 0.5 );
  beta_Rain ~ normal( 0 , 0.5 );
  beta_BreachDays_sub ~ normal( 0 , 0.5 );
  for ( i in 1:314 ) {
    SAV_nu[i] = a_SAV + beta_DO * DO[i] + beta_Temp * Temp[i] + k_U_bio * U[Zone[i]];
  }
  SAV ~ normal( SAV_nu , tau );
  for ( i in 1:314 ) {
    SC_mu[i] = a_SC + beta_Substrate * Substrate[i] + beta_DO * DO[i] + beta_SAV * SAV[i] + k_U_bio * U[Zone[i]];
    SC_mu[i] = inv_logit(SC_mu[i]);
  }
  SC_count ~ binomial( 1 , SC_mu );
  for ( i in 1:314 ) {
    SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i] + k_U_bio * U[Zone[i]];
    SB_mu[i] = inv_logit(SB_mu[i]);
  }
  SB_count ~ binomial( 1 , SB_mu );
  for ( i in 1:314 ) {
    Temp_nu[i] = a_Temp + beta_BreachDays_sub * BreachDays[i] + beta_Wind * Wind[i] + k_U * U[Zone[i]];
  }
  Temp ~ normal( Temp_nu , tau );
  for ( i in 1:314 ) {
    DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_Wind * Wind[i] + k_U * U[Zone[i]];
  }
  DO ~ normal( DO_nu , tau );
  for ( i in 1:314 ) {
    Breach_nu[i] = a_BreachDays + beta_Rain * Rain[i] + k_U * U[Zone[i]];
  }
  BreachDays ~ normal( Breach_nu , tau );
  sigma_errors ~ exponential( 1 );
  Rho_errors ~ lkj_corr( 2 );
  rain_mu = a_Rain;
  wind_mu = a_Wind;
  {
    array[314] vector[2] YY;
    vector[2] MU;
    MU = [ wind_mu , rain_mu ]';
    for ( j in 1:314 ) YY[j] = [ Wind[j] , Rain[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho_errors , sigma_errors) );
  }
  U ~ normal( 0 , 1 );
  sigma_Zone ~ exponential( 1 );
  Rho_Zone ~ lkj_corr( 2 );
  b_bar_Breach ~ normal( 0 , 0.5 );
  a_bar_Goby ~ normal( 0 , 0.5 );
  {
    array[3] vector[2] YY;
    vector[2] MU;
    MU = [ a_bar_Goby , b_bar_Breach ]';
    for ( j in 1:3 ) YY[j] = [ a_Goby[j] , b_BreachDays[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho_Zone , sigma_Zone) );
  }
  for ( i in 1:314 ) {
    mu[i] = a_Goby[Zone[i]] + b_BreachDays[Zone[i]] * BreachDays[i] + beta_SAV * SAV[i] + beta_DO * DO[i] + beta_Temp * Temp[i] + beta_SC_count * SC_count[i] + beta_Year * Year[i] + beta_Year_2 * Year_2[i] + beta_SB_count * SB_count[i] + beta_Micro * Micro[i] + beta_Substrate * Substrate[i] + beta_Temp_2 * Temp_2[i] + beta_Goby_lag * Goby_lag[i] + beta_BreachDays_2 * BreachDays_2[i] + Area[i];
    mu[i] = exp(mu[i]);
  }
  Goby ~ neg_binomial_2( mu , exp(phi) );
}
"

# write to stan file
f.RS_Breach.sim <- write_stan_file(stan_program.RS_Breach.sim)

# compile model 
# about 3 minutes
mod.sim.RS_Breach <- cmdstan_model(f.RS_Breach.sim)

