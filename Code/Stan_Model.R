library(rstan)

stan_program <- '
  data{
  array[299] int Breach_count;
  vector[299] Breach;
  vector[299] SC;
  vector[299] SB;
  array[299] int Year_int;
  array[299] int Goby;
  vector[299] Micro;
  vector[299] Year;
  array[299] int Zone;
  vector[299] Wind;
  array[299] int SB_count;
  array[299] int SC_count;
  vector[299] Area;
  array[299] int Substrate;
  vector[299] BreachDays;
  vector[299] Rain;
  vector[299] SAV;
  vector[299] Temp;
  vector[299] DO;
}
parameters{
  real mu_Zone;
  real<lower=0> tau_Zone;
  real beta_Substrate;
  real beta_Wind;
  real beta_BreachDays;
  real beta_Micro;
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
  real a_BreachDays;
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
  a_BreachDays ~ normal( 0 , 0.5 );
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
  beta_Micro ~ normal( 0 , 0.5 );
  beta_BreachDays ~ normal( 0 , 0.5 );
  beta_Wind ~ normal( 0 , 0.5 );
  beta_Substrate ~ normal( 0 , 0.5 );
  for ( i in 1:299 ) {
    SAV_nu[i] = a_SAV + beta_DO * DO[i] + beta_Temp * Temp[i];
  }
  SAV ~ normal( SAV_nu , tau );
  for ( i in 1:299 ) {
    Breach_nu[i] = a_BreachDays + beta_Rain * Rain[i];
  }
  BreachDays ~ normal( Breach_nu , tau );
  for ( i in 1:299 ) {
    SC_mu[i] = a_SC + beta_Substrate * Substrate[i] + beta_DO * DO[i] + beta_SAV * SAV[i] + Area[i];
    SC_mu[i] = exp(SC_mu[i]);
  }
  SC_count ~ neg_binomial_2( SC_mu , exp(SC_phi) );
  for ( i in 1:299 ) {
    SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i] + Area[i];
    SB_mu[i] = exp(SB_mu[i]);
  }
  SB_count ~ neg_binomial_2( SB_mu , exp(SB_phi) );
  for ( i in 1:299 ) {
    Temp_nu[i] = a_Temp + beta_BreachDays * BreachDays[i] + beta_Wind * Wind[i];
  }
  Temp ~ normal( Temp_nu , tau );
  for ( i in 1:299 ) {
    DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_BreachDays * BreachDays[i] + beta_Wind * Wind[i];
  }
  DO ~ normal( DO_nu , tau );
  tau_Zone ~ exponential( 1 );
  mu_Zone ~ normal( 0 , 5 );
  a_Goby ~ normal( mu_Zone , tau_Zone );
  for ( i in 1:299 ) {
    mu[i] = a_Goby[Zone[i]] + beta_Year * Year[i] + beta_SC_count * SC_count[i] + beta_SAV * SAV[i] + beta_SB_count * SB_count[i] + beta_DO * DO[i] + beta_Micro * Micro[i] + beta_BreachDays * BreachDays[i] + beta_Substrate * Substrate[i] + Area[i];
    mu[i] = exp(mu[i]);
  }
  Goby ~ neg_binomial_2( mu , exp(phi) );
}
'

dat$Zone <- as.numeric(dat$Zone)
dat$Substrate <- as.numeric(dat$Substrate)
stan_data <- dat

#run model

m1 <- stan(model_code = stan_program, data = stan_data, 
           cores = 4)

#runs in 3.5 minutes after compilation

summary(m1)
plot(m1, pars = c("beta_Year", 
                          "beta_DO", 
                          "beta_Substrate", 
                          "beta_Micro",
                          "beta_BreachDays", 
                          "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain",
                          "beta_Temp",
                          "beta_Wind"))
plot_model(m1)

df_of_draws <- as.data.frame(m1)

ext_fit <- as.data.frame(m1)
beta_BreachDays <- ext_fit$beta_BreachDays

gen_quant_r <- function(x) {
  lin_comb <- sample()
}





library(bayesplot) 

ppc_ribbon(
  y=dat$Goby/dat$Area,
        yrep=posterior_predict(m1),  ##does not work with stanfit
        x = dat$BreachDays,
  prob=0.5
)

library(shinystan)

launch_shinystan(m1)



dat %>%
  data_grid(BreachDays = seq_range(BreachDays, n = 101)) %>%
  add_predicted_draws(m1) %>%
  ggplot(aes(x = BreachDays, y = Goby/Area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "#08519C") +
  geom_point(data = dat, size = 2) +
  scale_fill_brewer()


