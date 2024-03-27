library(rstan)

#updated 2024-03-19
stan_program <- '
data{
     vector[301] Temp_2;
     vector[301] BreachDays_2;
     vector[301] Micro;
     vector[301] SAV_2;
     vector[301] Year_2;
     vector[301] Year;
    array[301] int Zone;
     vector[301] Wind;
    array[301] int SB_count;
    array[301] int SC_count;
     vector[301] Area;
    array[301] int Substrate;
     vector[301] BreachDays;
     vector[301] Rain;
     vector[301] SAV;
     vector[301] Temp;
     vector[301] DO;
}
parameters{
     real mu_Zone;
     real<lower=0> tau_Zone;
     real beta_Wind;
     real beta_Micro;
     real beta_Temp_2;
     real beta_Temp;
     real beta_Year_2;
     real beta_Year;
     real beta_SB_count;
     real a_Temp;
     real a_SAV;
     real a_SC;
     real a_SB;
     real a_DO;
     real a_BreachDays;
     vector[3] a_Goby;
     real beta_Rain;
     real beta_SC_count;
     real beta_SAV;
     real beta_SAV_2;
     real beta_DO;
     real beta_BreachDays;
     real beta_BreachDays_2;
     real beta_Substrate;
     real<lower=0> tau;
     real phi;
     real SB_phi;
     real SC_phi;
}
model{
     vector[301] mu;
     vector[301] DO_nu;
     vector[301] Temp_nu;
     vector[301] SB_mu;
     vector[301] SC_mu;
     vector[301] Breach_nu;
     vector[301] SAV_nu;
    SC_phi ~ normal( 1 , 5 );
    SB_phi ~ normal( 1 , 5 );
    phi ~ normal( 1 , 5 );
    tau ~ exponential( 1 );
    beta_Substrate ~ normal( 0.25 , 0.25 );
    beta_BreachDays_2 ~ normal( 0.25 , 0.25 );
    beta_BreachDays ~ normal( 0.25 , 0.25 );
    beta_DO ~ normal( 0.25 , 0.25 );
    beta_SAV_2 ~ normal( 0 , 0.25 );
    beta_SAV ~ normal( 0 , 0.25 );
    beta_SC_count ~ normal( -0.1 , 0.25 );
    beta_Rain ~ normal( 0.25 , 0.25 );
    a_Goby ~ normal( 0 , 0.75 );
    a_BreachDays ~ normal( 0 , 0.75 );
    a_DO ~ normal( 0 , 0.75 );
    a_SB ~ normal( 0 , 0.75 );
    a_SC ~ normal( 0 , 0.75 );
    a_SAV ~ normal( 0 , 0.75 );
    a_Temp ~ normal( 0 , 0.75 );
    beta_SB_count ~ normal( 0 , 0.75 );
    beta_Year ~ normal( 0 , 0.75 );
    beta_Year_2 ~ normal( 0 , 0.75 );
    beta_Temp ~ normal( 0 , 0.75 );
    beta_Temp_2 ~ normal( 0 , 0.75 );
    beta_Micro ~ normal( 0 , 0.75 );
    beta_Wind ~ normal( 0 , 0.75 );
    for ( i in 1:301 ) {
        SAV_nu[i] = a_SAV + beta_DO * DO[i] + beta_Temp * Temp[i];
    }
    SAV ~ normal( SAV_nu , tau );
    for ( i in 1:301 ) {
        Breach_nu[i] = a_BreachDays + beta_Rain * Rain[i];
    }
    BreachDays ~ normal( Breach_nu , tau );
    for ( i in 1:301 ) {
        SC_mu[i] = a_SC + beta_Substrate * Substrate[i] + beta_DO * DO[i] + beta_SAV * SAV[i] + Area[i];
        SC_mu[i] = exp(SC_mu[i]);
    }
    SC_count ~ neg_binomial_2( SC_mu , exp(SC_phi) );
    for ( i in 1:301 ) {
        SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i] + beta_SC_count * SC_count[i] + Area[i];
        SB_mu[i] = exp(SB_mu[i]);
    }
    SB_count ~ neg_binomial_2( SB_mu , exp(SB_phi) );
    for ( i in 1:301 ) {
        Temp_nu[i] = a_Temp + beta_BreachDays * BreachDays[i] + beta_Wind * Wind[i];
    }
    Temp ~ normal( Temp_nu , tau );
    for ( i in 1:301 ) {
        DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_Wind * Wind[i];
    }
    DO ~ normal( DO_nu , tau );
    tau_Zone ~ exponential( 1 );
    mu_Zone ~ normal( 0 , 0.5 );
    a_Goby ~ normal( mu_Zone , tau_Zone );
    for ( i in 1:301 ) {
        mu[i] = a_Goby[Zone[i]] + beta_Year * Year[i] + beta_Year_2 + Year_2[i] + beta_SC_count * SC_count[i] + beta_SAV * SAV[i] + beta_SAV_2 * SAV_2[i] + beta_SB_count * SB_count[i] + beta_DO * DO[i] + beta_Micro * Micro[i] + beta_BreachDays * BreachDays[i] + beta_BreachDays_2 * BreachDays_2[i] + beta_Substrate * Substrate[i] + beta_Wind * Wind[i] + beta_Temp * Temp[i] + beta_Temp_2 * Temp_2[i] + Area[i];
        mu[i] = exp(mu[i]);
    }
    Goby ~ neg_binomial_2( mu , exp(phi) );
}
'

stan_data <- as.data.frame(dat)

#run model

#remove rethinking packages to make sure rstan uses an rstan stanfit
m1 <- stan(model_code = stan_program, data = list(stan_data), 
           chains=3 , cores=parallel::detectCores() , iter=3000)
beepr::beep()
#runs in 3.5 minutes after compilation

summary(m1, digits = 2)
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


