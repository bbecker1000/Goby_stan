
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

summary(goby.brm.test)


#the model from 2_full_luxiury_gobies_2025...

# Goby model
bf_Goby <- bf(Goby ~ 
                Year + 
                Year_2 +
                SC_count +
                SAV +
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
                offset(Area) + #area already logged
                (1 | Zone), 
                family = negbinomial(link = "log", link_shape = "log"))

#DO model
bf_DO <- bf(DO ~ Temp + Wind + (1 | Zone))

#Temp model
bf_Temp <- bf(Temp ~ BreachDays + Wind + (1 | Zone))

#SB model as logistic
bf_SB_count <- bf(SB_count ~ DO + SAV + (1 | Zone), family = bernoulli(link = "logit"))

#SC model as logistic
bf_SC_count  <- bf(SC_count ~ Substrate + DO + SAV + (1 | Zone), family = bernoulli(link = "logit"))

#BreachDays as normal
bf_BreachDays <- bf(BreachDays ~ Rain)

#SAV model
bf_SAV <- bf(SAV ~ DO + Temp + (1 | Zone))

multi_formula <- bf_Goby + bf_DO + bf_Temp + bf_SB_count + bf_SC_count + bf_BreachDays + bf_SAV

get_prior(multi_formula, data = dat)


goby.brm <- brm(
  formula = multi_formula,
  data=dat, 
  chains=3, 
  cores=3,
  iter=4000,
  backend = "cmdstanr"
)

summary(goby.brm)


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






