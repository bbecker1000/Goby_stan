# Final Models
source("code/1_DataPrep.R")


library(rethinking)
library(cmdstanr)


##################

## updated 2024-01-19 with BreachDays
## updated 2024-01-22 with updated DAG
## updated 2024-01-23 with Micro
## updated 2024-01-30 with Zone*Wind interaction
## updated 2024-02-09 with breachdays and temp squared
## 2024-03-04 added 2023 data from darren
#### add year for trend.

#### -------------------------

## SC and SB dist from nb to logistic 2024-03-26
# the long tails on the SC and SB counts are hard for dgampois and pois to model (mass too high)
# so cap SC counts at 80 and SB counts at 500
# still has mean too high around 20 for neg bin and 15 for pois
# so use categories of y/n for each species

names(dat)
hist(dat$SC_count, n = 100)
dat$SC_count <- ifelse(dat$SC_count == 0, 0, 1)
hist(dat$SC_count)

hist(dat$SB_count, n = 100)
dat$SB_count <- ifelse(dat$SB_count == 0, 0, 1)
hist(dat$SB_count)



#### lag model 2025-06-07
#### make sure run lines 341-348 above to change to binary
#### fixing the polynomial priors to change sign!!
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07 <- ulam(
  alist(
    # Goby model
    Goby ~ dgampois(mu, exp(phi)), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <-
      a_Goby[Zone] + # random slope and random intercept
      beta_Year * Year +
      beta_Year_2 * Year_2 + # found "+" typo !! 2025-09-13
      beta_SC_count * SC_count + # updated
      beta_SAV * SAV +
      # beta_SAV_2*SAV_2 +
      beta_SB_count * SB_count + # competition
      beta_DO * DO +
      beta_Micro * Micro + # added 2024-01-23
      beta_BreachDays * BreachDays +
      beta_BreachDays_2 * BreachDays_2 + # include? as direct effect on goby?
      beta_Substrate * Substrate +
      beta_Wind * Wind +
      beta_Temp * Temp + # added 2024-02-09
      beta_Temp_2 * Temp_2 + # added 2024-02-09
      beta_Goby_lag * Goby_lag + # added 2024-03-24
      Area, # offset already logged

    # Zone RE model with regularizing priors
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), # the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5), # Hyperprior / hyperparameter
    tau_Zone ~ dexp(1), # Hyperprior / hyperparameter

    # DO model
    DO ~ normal(DO_nu, tau),
    DO_nu <-
      a_DO +
      beta_Temp * Temp +
      beta_Wind * Wind,

    # Temp model
    Temp ~ normal(Temp_nu, tau),
    Temp_nu <-
      a_Temp +
      beta_BreachDays * BreachDays +
      beta_Wind * Wind,

    # SB model as logistic
    SB_count ~ dbinom(1, SB_mu),
    logit(SB_mu) <-
      a_SB +
      beta_DO * DO +
      beta_SAV * SAV,

    # SC model as logistic
    SC_count ~ dbinom(1, SC_mu),
    logit(SC_mu) <-
      a_SC +
      beta_Substrate * Substrate +
      beta_DO * DO +
      beta_SAV * SAV,

    # BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain * Rain,

    # SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <-
      a_SAV +
      beta_DO * DO +
      beta_Temp * Temp,

    # fixed effects priors
    c(
      a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp,
      beta_Wind
    ) ~ normal(0, 0.5), # regularizing
    beta_Temp_2 ~ normal(-0.10, 0.25), # goldilocks (corrected to negative 2025-06-07)
    beta_Micro ~ normal(0.25, 0.25), # more fish = more micro
    beta_Rain ~ normal(0.25, 0.25), # more rain = more goby
    beta_SC_count ~ normal(-0.10, 0.25), # sculpins eat goby larvae
    beta_SAV ~ normal(0.00, 0.25), # goldilocks
    beta_SAV_2 ~ normal(-0.10, 0.25), # goldilocks (corrected to negative 2025-06-07)
    beta_DO ~ normal(0.25, 0.25), # more DO good
    beta_BreachDays ~ normal(0.25, 0.25), # more breach good
    beta_BreachDays_2 ~ normal(-0.10, 0.25), # goldilocks (corrected to negative 2025-06-07)
    beta_Substrate ~ normal(0.25, 0.25), # coarser habitat better
    beta_Goby_lag ~ normal(0.25, 0.25), # added 2024-03-24

    tau ~ exponential(1),
    phi ~ dnorm(1, 5) # was ("dexp(1)")
  ),
  data = dat, chains = 3, cores = parallel::detectCores(), iter = 4000, # high r-hat with 3k iter
  cmdstan = TRUE # FALSE to get stanfit object, but has a problem solution if false
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1 - t0
runtime


precis(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07,
  pars = c(
    "beta_Year",
    "beta_Year_2",
    "beta_Goby_lag",
    "beta_DO",
    "beta_Substrate",
    "beta_Micro",
    "beta_BreachDays",
    "beta_BreachDays_2",
    "beta_SB_count",
    "beta_SC_count",
    "beta_SAV",
    "beta_Rain",
    "beta_Temp",
    "beta_Temp_2",
    "beta_Wind" # ,
    #  "beta_Zone,"#
    #  "beta_ZW"
  ),
  xlab = "Beta Coefficient",
  main = "network model-Breach Direct"
)


save(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07,
  file = "Output/Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07.RData"
)

# load(file = "Output/Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07.RData")

# explort for use with cmdstan
stanmodelcode <- rethinking::stancode(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07)


# for use in 3_plotting
fit <- Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07


rethinking::traceplot(fit)

plot(fit, depth = 2)

#### -------------
