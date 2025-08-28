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
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      #beta_ZW*Wind*Zone + # interaction between wind and Zone
                           # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # was 0,5
    tau_Zone ~ dexp(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,

    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as neg.bin
    SB_count ~ dgampois( SB_mu, exp(SB_phi) ), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # logged offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, exp(SC_phi) ), 
    log(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind
      # beta_ZW
      )                 ~ normal( 0 , 0.5 ),      # uninformed
      beta_Temp_2       ~ normal( 0.10 , 0.25),   # goldilocks 
      beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
      beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
      beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
      beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
      beta_SAV_2        ~ normal( 0.10 , 0.25 ),  # goldilocks 
      beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
      beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
      beta_BreachDays_2 ~ normal( 0.10 , 0.25 ),  # goldilocks
      beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5), #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1,5)
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=3000 , 
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS,
     pars = c("beta_Year", 
              "beta_Year_2", 
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
              "beta_Wind"#,
            #  "beta_Zone,"#
            #  "beta_ZW"
            ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS, file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS)

fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS

plot(fit, depth = 2)

#rethinking::pairs(fit)

#causal path calculations



## SC and SB dist from nb to logistic 2024-03-26
# the long tails on the SC and SB counts are hard for dgampois and pois to model (mass too high)
# so cap SC counts at 80 and SB counts at 500  
# still has mean too high around 20 for neg bin and 15 for pois
# so use categories of y/n for each species

names(dat)
hist(dat$SC_count, n=100)
dat$SC_count <- ifelse(dat$SC_count == 0, 0, 1)
hist(dat$SC_count)

hist(dat$SB_count, n=100)
dat$SB_count <- ifelse(dat$SB_count == 0, 0, 1)
hist(dat$SB_count)


#### lag model _Pre 2025
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      beta_Goby_lag*Goby_lag + #added 2024-03-24
      #beta_ZW*Wind*Zone + # interaction between wind and Zone
      # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model with regularizing priors
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # Hyperprior / hyperparameter
    tau_Zone ~ dexp(1),       # Hyperprior / hyperparameter
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as logistic
    SB_count ~ dbinom( 1, SB_mu ), 
    logit(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #SC model as logistic
    SC_count ~ dbinom( 1, SC_mu ),   
    logit(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind
      # beta_ZW
    )                 ~ normal( 0 , 0.5 ),      # regularizing
    beta_Temp_2       ~ normal( 0.10 , 0.25),   # goldilocks 
    beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
    beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
    beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
    beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
    beta_SAV_2        ~ normal( 0.10 , 0.25 ),  # goldilocks 
    beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
    beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
    beta_BreachDays_2 ~ normal( 0.10 , 0.25 ),  # goldilocks
    beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  # added 2024-03-24
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5) #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    #c(SC_phi, SB_phi) ~ dnorm(1,5)
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=4000 , #high r-hat with 3k iter
  cmdstan=TRUE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime


precis(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag,
     pars = c("beta_Year", 
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
              "beta_Wind"#,
              #  "beta_Zone,"#
              #  "beta_ZW"
     ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag, 
     file = "Output/Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag)

fit <- Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag

plot(fit, depth = 2)

####-------------------------

## SC and SB dist from nb to logistic 2024-03-26
# the long tails on the SC and SB counts are hard for dgampois and pois to model (mass too high)
# so cap SC counts at 80 and SB counts at 500  
# still has mean too high around 20 for neg bin and 15 for pois
# so use categories of y/n for each species

names(dat)
hist(dat$SC_count, n=100)
dat$SC_count <- ifelse(dat$SC_count == 0, 0, 1)
hist(dat$SC_count)

hist(dat$SB_count, n=100)
dat$SB_count <- ifelse(dat$SB_count == 0, 0, 1)
hist(dat$SB_count)



#### lag model 2025-06-07
#### make sure run lines 341-348 above to change to binary
#### fixing the polynomial priors to change sign!!
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07 <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      beta_Goby_lag*Goby_lag + #added 2024-03-24
      #beta_ZW*Wind*Zone + # interaction between wind and Zone
      # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model with regularizing priors
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # Hyperprior / hyperparameter
    tau_Zone ~ dexp(1),       # Hyperprior / hyperparameter
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as logistic
    SB_count ~ dbinom( 1, SB_mu ), 
    logit(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #SC model as logistic
    SC_count ~ dbinom( 1, SC_mu ),   
    logit(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind
      # beta_ZW
    )                 ~ normal( 0 , 0.5 ),      # regularizing
    beta_Temp_2       ~ normal(-0.10 , 0.25),   # goldilocks (corrected to negative 2025-06-07)
    beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
    beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
    beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
    beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
    beta_SAV_2        ~ normal(-0.10 , 0.25 ),  # goldilocks (corrected to negative 2025-06-07)
    beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
    beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
    beta_BreachDays_2 ~ normal(-0.10 , 0.25 ),  # goldilocks (corrected to negative 2025-06-07)
    beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  # added 2024-03-24
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5) #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    #c(SC_phi, SB_phi) ~ dnorm(1,5)
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=4000 , #high r-hat with 3k iter
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime


precis(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07,
     pars = c("beta_Year", 
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
              "beta_Wind"#,
              #  "beta_Zone,"#
              #  "beta_ZW"
     ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07, 
     file = "Output/Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07.RData")
#load(file = "Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07)

fit <- Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07

plot(fit, depth = 2)

####-------------





#### lag model
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      beta_Goby_lag*Goby_lag + #added 2024-03-24
      #beta_ZW*Wind*Zone + # interaction between wind and Zone
      # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # was 0,5
    tau_Zone ~ dexp(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as neg.bin
    SB_count ~ dgampois( SB_mu, exp(SB_phi) ), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # logged offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, exp(SC_phi) ),   
    log(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind
      # beta_ZW
    )                 ~ normal( 0 , 0.5 ),      # uninformed
    beta_Temp_2       ~ normal(-0.10 , 0.25),   # goldilocks 
    beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
    beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
    beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
    beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
    beta_SAV_2        ~ normal(-0.10 , 0.25 ),  # goldilocks 
    beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
    beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
    beta_BreachDays_2 ~ normal( 0.10 , 0.25 ),  # goldilocks
    beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  #added 2024-03-24
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5), #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1,5)
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=4000 , #high r-hat with 3k iter
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime


precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag,
     pars = c("beta_Year", 
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
              "beta_Wind"#,
              #  "beta_Zone,"#
              #  "beta_ZW"
     ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag, 
     file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag)

fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag

plot(fit, depth = 2)


#### lag model
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      beta_Goby_lag*Goby_lag + #added 2024-03-24
      #beta_ZW*Wind*Zone + # interaction between wind and Zone
      # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # was 0,5
    tau_Zone ~ dexp(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as neg.bin
    SB_count ~ dpois( SB_mu ), # dgampois( SB_mu, exp(SB_phi) ), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # logged offset
    
    #SC model as neg.bin
    SC_count ~ dpois( SC_mu ),   #dgampois( SC_mu, exp(SC_phi) ),  
    log(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind
      # beta_ZW
    )                 ~ normal( 0 , 0.5 ),      # uninformed
    beta_Temp_2       ~ normal(-0.10 , 0.25),   # goldilocks 
    beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
    beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
    beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
    beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
    beta_SAV_2        ~ normal(-0.10 , 0.25 ),  # goldilocks 
    beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
    beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
    beta_BreachDays_2 ~ normal( 0.10 , 0.25 ),  # goldilocks
    beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  #added 2024-03-24
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5) #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    #c(SC_phi, SB_phi) ~ dnorm(1,5) # removed if using poisson. 2024-03-25
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=4000 , #high r-hat with 3k iter
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime


precis(Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag,
     pars = c("beta_Year", 
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
              "beta_Wind"#,
              #  "beta_Zone,"#
              #  "beta_ZW"
     ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag, 
     file = "Output/Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag)

fit <- Goby.m2.year.DAG.SC.SB_POISSON.BreachDays.direct.RS.lag

plot(fit, depth = 2)

#####-end



####---
#Add wind*zone interaction

t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp( phi )), # was "phi" exp(phi) to constrain to positive reals
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      #beta_SAV_2*SAV_2 + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +
      beta_Wind*Wind +  #add a wind slope conditional on [Zone]
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  #  added 2024-02-09
      beta_Goby_lag*Goby_lag + #added 2024-03-24
      beta_ZW*Wind*Zone +  # interaction between wind and Zone
      beta_Zone*Zone +   # need for WZ interaction
      # slope of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ dnorm(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ dnorm(0, 0.5),  # was 0,5
    tau_Zone ~ dexp(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as neg.bin
    SB_count ~ dgampois( SB_mu, exp(SB_phi) ), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # logged offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, exp(SC_phi) ), 
    log(SC_mu) <- 
      a_SC + 
      beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_SB_count,
      beta_Year,
      beta_Year_2,
      beta_Temp, 
      beta_Wind,
      beta_ZW,
      beta_Zone
    )                 ~ normal( 0 , 0.5 ),      # uninformed
    beta_Temp_2       ~ normal(-0.10 , 0.25),   # goldilocks 
    beta_Micro        ~ normal( 0.25 , 0.25 ),  # more fish = more micro
    beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
    beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
    beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
    beta_SAV_2        ~ normal(-0.10 , 0.25 ),  # goldilocks 
    beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
    beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
    beta_BreachDays_2 ~ normal( 0.10 , 0.25 ),  # goldilocks
    beta_Substrate    ~ normal( 0.25 , 0.25 ),  # coarser habitat better
    beta_Goby_lag     ~ normal( 0.25 , 0.25 ),  #added 2024-03-24
    
    tau ~ exponential(1),
    #phi ~ dexp(1),
    phi ~ dnorm(1,5), #was ("dexp(1)")
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1,5)
    #c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=3 , cores=parallel::detectCores() , iter=4000 , #high r-hat with 3k iter
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime


precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag,
     pars = c("beta_Year", 
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
              "beta_Wind"
     ),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


save(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag, 
     file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.lag.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag.RData")


rethinking::stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag)

fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS.WZ.lag

plot(fit, depth = 2)





############-----------------


####----------
#### try with imputed data
### need priors for data and the n = 349 data list 

### add year for trend.
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      #beta_Temp*Temp +
      beta_BreachDays*BreachDays +   # include? as direct effect on goby?
      #beta_Substrate*Substrate +    # removed since now RE
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #Substrate RE model
    a_Goby[Substrate] ~ normal(mu_Substrate, tau_Substrate),
    mu_Substrate ~ normal(0, 5),
    tau_Substrate ~ exponential(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_BreachDays*BreachDays + 
      beta_Wind*Wind,
    
    #SB model as neg.bin
    SB_count ~ dgampois( SB_mu, exp(SB_phi)), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, exp(SC_phi)), 
    log(SC_mu) <- 
      a_SC + 
      #beta_Substrate*Substrate +
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #BreachDays as normal
    BreachDays ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_BreachDays +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO +
      beta_Temp*Temp,
    
    #missing data priors
    c(Temp, DO) ~ normal (0,1),
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
      beta_Micro,
      #beta_Breach,
      beta_BreachDays,
      beta_Wind#,
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat.missing , chains=4 , cores=4 , iter=3000 , cmdstan=TRUE
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing, 
     pars = c("beta_Year", "beta_Wind", "beta_BreachDays", "beta_Temp",
              "beta_Micro",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient")

#plot(Goby.m1, depth = 2)
summary(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing)

save(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing, 
     file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.missing.RData")
####









