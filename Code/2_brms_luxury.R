
library(brms)

# 2025-08-28

#some logistic values for other fish
names(dat)
hist(dat$SC_count, n=100)
dat$SC_count <- ifelse(dat$SC_count == 0, 0, 1)
hist(dat$SC_count)

hist(dat$SB_count, n=100)
dat$SB_count <- ifelse(dat$SB_count == 0, 0, 1)
hist(dat$SB_count)




#working test model with 2 formulas and distributions
bf_Goby_test <- bf(Goby ~ Year + BreachDays + (1 | Zone), 
              family = negbinomial(link = "log", link_shape = "log"))

bf_BreachDays_test <- bf(BreachDays ~ Rain, family = gaussian)

multi_formula_test <- bf_Goby + bf_BreachDays

goby.brm.test <- brm(
  formula = multi_formula_test,
  data=dat, 
  chains=3, 
  iter=3000,
  backend = "cmdstanr"
)




#the model from 2_full_luxiury_gobies_2025...

# Goby model
bf_Goby <- bf(Goby ~ 
                beta_Year * Year + 
                beta_Year_2 * Year_2 +
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
                beta_Temp_2 * Temp_2+
                beta_Goby_lag * Goby_lag +
                offset(Area) + #area already logged
                (1 | Zone)# allows to assign shared parameters


get_prior(multi_formula_joint, data = dat)

PRIORS <- c(
  prior(normal(0, 2), nlpar = "beta_Year"),
  prior(normal(0, 2), nlpar = "beta_Year_2"),
  prior(normal(0, 2), nlpar = "beta_SC_count"),
        prior(normal(0, 2), nlpar = "beta_SAV"),
        prior(normal(0, 2), nlpar = "beta_SB_count"),
              prior(normal(0, 2), nlpar = "beta_DO"),
              prior(normal(0, 2), nlpar = "beta_Micro"),
              prior(normal(0, 2), nlpar = "beta_BreachDays",
                    prior(normal(0, 2), nlpar = "beta_BreachDays_2"),
                    prior(normal(0, 2), nlpar = "beta_Substrate"),
                    prior(normal(0, 2), nlpar = "beta_Wind",
                          prior(normal(0, 2), nlpar = "beta_Temp"),
                          prior(normal(0, 2), nlpar = "beta_Temp_2"),
                          prior(normal(0, 2), nlpar = "beta_Goby_lag")
                    )))

              
              
goby.brm <- brm(
  formula = multi_formula_joint,
  family = FAMILY,
  prior = PRIORS,
  data=dat, 
  chains=3, 
  cores=3,
  iter=4000,
  backend = "cmdstanr"
)

summary(goby.brm)


tidybayes::tidy_draws(goby.brm)

# #fixed effects priors
# c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
#   beta_SB_count,
#   beta_Year,
#   beta_Year_2,
#   beta_Temp, 
#   beta_Wind
#   # beta_ZW
# )                 ~ normal( 0 , 0.5 ),      # regularizing
# beta_Temp_2       ~ normal(-0.10 , 0.25),   # goldilocks (corrected to negative 2025-06-07)
# beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
# beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
# beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
# beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
# beta_SAV_2        ~ normal(-0.10 , 0.25 ),  # goldilocks (corrected to negative 2025-06-07)
# beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
# beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
# beta_BreachDays_2 ~ normal(-0.10 , 0.25 ),  # goldilocks (corrected to negative 2025-06-07)
# beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
# beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  # added 2024-03-24
# 
# tau ~ exponential(1),
# #phi ~ dexp(1),
# phi ~ dnorm(1,5) #was ("dexp(1)")
# #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
# #c(SC_phi, SB_phi) ~ dnorm(1,5)
# #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
# ), 


#

#gemini and kurtz chap 14

# Load the libraries and data
library(brms)
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# Standardize the variables for the model
d$D_obs <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A_obs <- standardize(d$MedianAgeMarriage)

# Define the two submodels
bf_divorce_unconstrained <- brms::bf(
  D_obs ~ b_int_d + b_A * A_obs + b_M_D * M,
  nl = TRUE
)

bf_marriage_age_unconstrained <- brms::bf(
  A_obs | mi() ~ b_int_a + b_M_A * M,
  nl = TRUE
)

# Set the priors
priors_unconstrained <-
  prior(normal(0, 0.2), nlpar = "b_int_d") +
  prior(normal(0, 0.5), nlpar = "b_A") +
  prior(normal(0, 0.5), nlpar = "b_M_D") + # <-- Prior for effect on Divorce
  prior(normal(0, 0.2), nlpar = "b_int_a") +
  prior(normal(0, 0.5), nlpar = "b_M_A")   # <-- Prior for effect on Age

# This model would estimate b_M_D and b_M_A separately

# Define the two submodels with a SHARED parameter name
bf_divorce_shared <- brms::bf(
  D_obs ~ b_int_d + b_A * A_obs + b_M_shared * M, # <-- Using b_M_shared
  nl = TRUE
)

bf_marriage_age_shared <- brms::bf(
  A_obs | mi() ~ b_int_a + b_M_shared * M, # <-- Using b_M_shared
  nl = TRUE
)

# Set the priors, now with only one prior for the shared parameter
priors_shared <-
  prior(normal(0, 0.2), nlpar = "b_int_d") +
  prior(normal(0, 0.5), nlpar = "b_A") +
  prior(normal(0, 0.2), nlpar = "b_int_a") +
  prior(normal(0, 0.5), nlpar = "b_M_shared") # <-- SINGLE prior for the shared effect

# Fit the model
fit_shared_parameter <- brm(
  bf_divorce_shared + bf_marriage_age_shared,
  data = d,
  prior = priors_shared,
  chains = 4, cores = 4, seed = 123
)


