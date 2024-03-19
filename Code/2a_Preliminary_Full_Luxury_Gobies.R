# Slow Building of Ulam Stan Models

source("code/1_DataPrep.R")


#test lmer
m1.lmer <- glmer(Goby ~ Year + 
                   Rain + I(SC_count/Area) + SAV + I(SB_count/Area) + DO + Temp + 
                   BreachDays +
                   Micro +
                   #Breach_count + 
                   Wind + 
                   Substrate +
                   (1|Zone), #think about random slopes
                 data = dat, family = negative.binomial(0.6), offset = Area)
summary(m1.lmer)
plot(m1.lmer)
plot_model(m1.lmer, type = "pred")


#test gamm
library(mgcv)
m1.gamm <- gamm(Goby ~ Year + 
                  s(Rain, bs="cr") + 
                  I(SC_count/Area) + SAV + I(SB_count/Area) + DO + Temp + 
                  s(BreachDays, bs="cr") +
                  Micro +
                  #Breach_count + 
                  Wind + #maybe smooth
                  Substrate,
                random=(list(Zone=~1)), 
                data = dat, family = negative.binomial(0.6),
                niterPQL=50)
summary(m1.gamm$gam)
plot(m1.gamm$gam)

plot_model(m1.gamm)


library(rstanarm)
m1.stan_glm <- stan_glmer(Goby ~ Year +
                            Rain + SC_count + SAV + SB_count + DO + Temp + BreachDays + 
                            Substrate + (1|Zone), data = dat, family = neg_binomial_2(),
                          chains = 4, cores = 4, iter = 1000)
beepr::beep(0)
summary(m1.stan_glm, digits = 3)
plot_model(m1.stan_glm, type = "pred")


### model with no network
set.seed(321)
Goby.m1.no.network <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  
      beta_Year*Year +
      beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      beta_Temp*Temp +
      beta_BreachDays*BreachDays + 
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #fixed effects priors
    c(a_Goby, a_BreachDays, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC,
      beta_SAV, 
      beta_SB, 
      beta_DO,
      beta_Temp, 
      beta_BreachDays,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ normal( 0 , 10 ),
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE)

beepr::beep(0)

precis(Goby.m1.no.network, depth = 2)
plot(Goby.m1.no.network,
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_Breach", "beta_Temp",
              "beta_DO", "beta_SB", "beta_SC", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient",
     main = "no network model")
#plot(Goby.m1, depth = 2)
summary(Goby.m1.no.network)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m1.no.network)

# "no network" ulam coefficients are very close to stan_glmer coefficients - Yay!
# the Zone random effect intercepts are different, but have similar intervals 

save(Goby.m1.no.network, file = "Output/Goby.m1.no.network.RData")

par(mfrow = c(1,1))


#####------------try for breach = 0,1,2,
##              SB and SC = neg.bin, but made coefficients really tiny

### refining to get closer to Dag
### No year to just get causes

##### - SB and SC neg-bin, breach = poisson

### refining to get closer to Dag
### No year to just get causes

t0 <- Sys.time()

Goby.m2.DAG.SC.SB_counts.Breach.pois <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      #a_Goby[Substrate] + #random effect 
      beta_Year*Year +
      #beta_Rain*Rain + 
      beta_Micro*Micro +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      #beta_Temp*Temp +
      beta_Breach_count*Breach_count + 
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #Substrate RE model
    # a_Goby[Substrate] ~ normal(mu_Substrate, tau_Substrate),
    # mu_Substrate ~ normal(0, 5),
    # tau_Substrate ~ exponential(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Breach_count*Breach_count + 
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_Breach_count*Breach_count + 
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
      beta_DO*DO + 
      beta_SAV*SAV + 
      Area,  # logged offset
    
    #Breach model poisson
    Breach_count ~ poisson(Breach_nu),
    log(Breach_nu) <-
      a_Breach +
      beta_Rain*Rain,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_DO*DO,
    
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
      beta_Micro,
      #beta_Breach,
      beta_Breach_count,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=4 , cores=4 , iter=2000 , cmdstan=TRUE
)

beepr::beep(0)

t1 <- Sys.time()
runtime <- t1-t0
runtime


par(mfrow = c(1,1))

precis(Goby.m2.DAG.SC.SB_counts.Breach.pois, depth = 2)
plot(Goby.m2.DAG.SC.SB_counts.Breach.pois, 
     pars = c("beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta",
     main = "Full model")
#plot(Goby.m1, depth = 2)
summary(Goby.m2.DAG.SC.SB_counts.Breach.pois)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.DAG.SC.SB_counts.Breach.pois)
par(mfrow = c(1,2))

### Simulate Goby with intervention of Wind * Temp --> Goby------------------
# set Wind to -0.25

post <- extract.samples(Goby.m2.DAG.SC.SB_counts.Breach.pois)

#Effect of Rain on Breach
quantile( with(post, beta_Rain*beta_Breach_count), probs = c(0.1, 0.5, 0.9)) ## ns
#effect of breach on DO
quantile( with(post, beta_DO*beta_Breach_count), probs = c(0.1, 0.5, 0.9)) ## ns
#effect of SAV on SC
quantile( with(post, beta_SAV*beta_SC_count), probs = c(0.1, 0.5, 0.9))  ## neg

quantile( with(post, beta_Temp*beta_SAV), probs = c(0.1, 0.5, 0.9))  ## positive

quantile( with(post, beta_DO*beta_Temp), probs = c(0.1, 0.5, 0.9)) ## neg

quantile( with(post, beta_Wind, beta_DO), probs = c(0.1, 0.5, 0.9)) ## positive




## updated 2024-01-19 with BreachDays
## updated 2024-01-22 with updated DAG 

#### add year for trend.
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      #a_Goby[Substrate] + #random effect 
      beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      #beta_Temp*Temp +
      # beta_Breach*BreachDays +   # include? as direct effect on goby?
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #Substrate RE model
    # a_Goby[Substrate] ~ normal(mu_Substrate, tau_Substrate),
    # mu_Substrate ~ normal(0, 5),
    # tau_Substrate ~ exponential(1),
    
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
      beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
      #beta_Breach,
      beta_BreachDays,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=4 , cores=4 , iter=1000 , cmdstan=TRUE
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_BreachDays", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient", 
     main = "network model")

#plot(Goby.m1, depth = 2)
summary(Goby.m2.year.DAG.SC.SB_counts.Breach.pois)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.year.DAG.SC.SB_counts.Breach.pois)


save(Goby.m2.year.DAG.SC.SB_counts.BreachDays, file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.RData")


## updated 2024-01-19 with BreachDays
## updated 2024-01-22 with updated DAG 
## updated 2024-01-23 with Micro
#### add year for trend.
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      #a_Goby[Substrate] + #random effect 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      #beta_Rain*Rain + 
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      beta_Micro*Micro + #added 2024-01-23
      #beta_Temp*Temp +
      beta_BreachDays*BreachDays + 
      beta_BreachDays_2*BreachDays_2 +  # include? as direct effect on goby?
      beta_Substrate*Substrate +    # remove if RE
      beta_Wind*Wind +
      beta_Zone*Zone + # interaction between wind and Zone
      #Zone + #interaction of wind piling up algae
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone), #the "1" adds varying slope.
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #Substrate RE model
    # a_Goby[Substrate] ~ normal(mu_Substrate, tau_Substrate),
    #  mu_Substrate ~ normal(0, 5),
    #  tau_Substrate ~ exponential(1),
    
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
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, exp(SC_phi)), 
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
      beta_BreachDays_2,
      beta_Wind,#,
      beta_Zone,
      beta_Substrate
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=4 , cores=parallel::detectCores() , iter=5000 , 
  cmdstan=FALSE # to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

precis(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct,
     pars = c("beta_Year", 
              "beta_DO", 
              "beta_Substrate", 
              "beta_Micro",
              "beta_BreachDays", 
              "beta_BreachDays_2",
              "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain",
              "beta_Temp",
              "beta_Wind"),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct)

