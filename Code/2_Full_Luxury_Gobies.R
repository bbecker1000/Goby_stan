# Models
source("1_DataPrep.R")

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
    Goby ~ dgampois( mu, phi ), #, 
    log(mu) <- 
      a_Goby[Zone] + #random slope and random intercept 
      beta_Year*Year +
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
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
    SB_count ~ dgampois( SB_mu, SB_phi), 
    log(SB_mu) <- 
      a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV + 
      beta_SC_count*SC_count +  #added per DF 2024-01-25
      Area,  # logged offset
    
    #SC model as neg.bin
    SC_count ~ dgampois( SC_mu, SC_phi), 
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
      beta_Temp, 
      beta_Temp_2,
      beta_Micro,
      beta_Wind
      # beta_ZW
      ) ~ 
      normal( 0 , 0.75 ),

      beta_Rain         ~ normal( 0.25 , 0.25 ),
      beta_SC_count     ~ normal(-0.10 , 0.25 ),
      beta_SAV          ~ normal( 0.25 , 0.25 ),
      beta_DO           ~ normal( 0.25 , 0.25 ),
      beta_BreachDays   ~ normal( 0.25 , 0.25 ),
      beta_BreachDays_2 ~ normal( 0.25 , 0.25 ),
      beta_Substrate    ~ normal( 0.25 , 0.25 ),
    
    tau ~ exponential(1),
    phi ~ dexp(1), 
    #phi ~ dnorm( 1, 3 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(1)  # use dexp(100) if not neg.bin
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
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RData")


stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS)




#causal path calculations

post <- extract.samples(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS)
PROBS = c(0.11, 0.5, 0.9) ##89 % CIs
#Effect of Rain --> Breach --> DO --> Goby
Rain_Breach_DO<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays*beta_DO), probs = PROBS)) ## NS
#Effect of Breach --> DO --> Goby
Breach_DO<-as_tibble(quantile( with(post, beta_BreachDays*beta_DO), probs = PROBS)) ##  NS
#Effect of Wind --> Breach --> DO --> Goby
Wind_DO<-as_tibble(quantile( with(post, beta_Wind*beta_DO), probs = PROBS)) ## ns
#Effect of Wind --> Temp --> DO --> Goby
Wind_Temp_DO<-as_tibble(quantile( with(post, beta_Wind*beta_Temp*beta_DO), probs = PROBS)) ## ns  
#Effect of Wind --> Temp --> SAV --> Goby
Wind_DO_SAV<-as_tibble(quantile( with(post, beta_Wind*beta_SAV*beta_DO), probs = PROBS)) ## ns  
#Effect of Wind --> DO --> SB --> Goby
Wind_DO_SC<-as_tibble(quantile( with(post, beta_Wind*beta_DO*beta_SC_count), probs = PROBS)) ## ns 
#Effect of DO --> SB --> Goby
DO_SB<-as_tibble(quantile( with(post, beta_DO*beta_SB_count), probs = PROBS)) ## ns 

Rain_Breach_Temp<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays*beta_Temp), probs = PROBS)) ## NS
#Effect of Breach --> Temp --> Goby
Rain_Breach_2<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays_2), probs = PROBS)) ## NS
#Effect of Rain --> Breach --> Goby

Breach_Temp<-as_tibble(quantile( with(post, beta_BreachDays*beta_Temp), probs = PROBS)) ##  NS
#Effect of Wind --> Breach --> Temp --> Goby
Wind_Temp<-as_tibble(quantile( with(post, beta_Wind*beta_Temp), probs = PROBS)) ## ns

#Effect of Wind --> Temp --> SAV --> Goby
Wind_Temp_SAV<-as_tibble(quantile( with(post, beta_Wind*beta_SAV*beta_Temp), probs = PROBS)) ## ns  
#Effect of Wind --> Temp --> SB --> Goby
Wind_Temp_SC<-as_tibble(quantile( with(post, beta_Wind*beta_Temp*beta_SC_count), probs = PROBS)) ## ns 
#Effect of Temp --> SB --> Goby
Temp_SB<-as_tibble(quantile( with(post, beta_Temp*beta_SB_count), probs = PROBS)) ## ns 
Rain_Temp_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_Temp*beta_SAV*beta_SB_count), probs = PROBS))

Temp<-as_tibble(quantile( with(post, beta_Temp), probs = PROBS)) ## ns 
Temp_2<-as_tibble(quantile( with(post, beta_Temp_2), probs = PROBS)) ## ns 

Breach_2<-as_tibble(quantile( with(post, beta_BreachDays_2), probs = PROBS)) ## ns 

Year<-as_tibble(quantile( with(post, beta_Year), probs = PROBS)) ## ns 
Substrate_SC<-as_tibble(quantile( with(post, beta_Substrate*beta_SC_count), probs = PROBS)) ## negative
Substrate<-as_tibble(quantile( with(post, beta_Substrate), probs = PROBS))
Rain_DO_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_DO*beta_SAV*beta_SB_count), probs = PROBS))

Breach<-as_tibble(quantile( with(post, beta_BreachDays), probs = PROBS)) #positive
Rain_Micro<-as_tibble(quantile( with(post, beta_Rain*beta_Micro), probs = PROBS))
Micro<-as_tibble(quantile( with(post, beta_Micro), probs = PROBS))
SAV_SB<-as_tibble(quantile( with(post, beta_SAV*beta_SB_count), probs = PROBS))
SAV_SC<-as_tibble(quantile( with(post, beta_SAV*beta_SC_count), probs = PROBS))
SC<-as_tibble(quantile( with(post, beta_SC_count), probs = PROBS))
SB<-as_tibble(quantile( with(post, beta_SB_count), probs = PROBS))
Breach_DO_SC<-as_tibble(quantile( with(post, beta_BreachDays*beta_DO*beta_SC_count), probs = PROBS)) ##  NS

# 2024-02-01
# ADD RAIN --> BREACH



#names for tibble
names<- c("Rain_Breach_DO", "Rain_Breach_2", "Breach_DO", "Wind_DO", "Wind_Temp_DO", "Wind_DO_SAV", "Wind_DO_SC",
          "DO_SB", "Year", "Substrate_SC", "Substrate", "Rain_DO_SAV_SB", "Breach", 
          "DO_SB", "Year", "Substrate_SC", "Rain_DO_SAV_SB", "Breach", 
          "Micro", "SAV_SB", "SAV_SC", "SC", "SB","Breach_DO_SC", "Rain_Breach_Temp",
          #"Zone", "Zone_Wind_int", 
          "Breach_Temp",
          "Wind_Temp",
          "Wind_Temp_SAV",
          "Wind_Temp_SC",
          "Temp_SB",
          "Rain_Temp_SAV_SB",
          "Temp",
          "Temp_2",
          "Breach_2")
#add probabilities
plot.posteriors<-rbind(Rain_Breach_DO, Rain_Breach_2, Breach_DO, Wind_DO, Wind_Temp_DO, Wind_DO_SAV, Wind_DO_SC,
        DO_SB, Year, Substrate_SC, Substrate, Rain_DO_SAV_SB, Breach, 
        DO_SB, Year, Substrate_SC, Rain_DO_SAV_SB, Breach, 
        Micro, SAV_SB, SAV_SC, SC, SB, Breach_DO_SC, Rain_Breach_Temp,
        #Zone, Zone_Wind_int,
        Breach_Temp,
        Wind_Temp,
        Wind_Temp_SAV,
        Wind_Temp_SC,
        Temp_SB,
        Rain_Temp_SAV_SB,
        Temp,
        Temp_2,
        Breach_2)
#add names
plot.posteriors$names <- rep(names, each=3)
#add probabilities names
plot.posteriors$probability <- rep(c("lower", "median", "upper"), times = length(names))

plot.posteriors.wide <- plot.posteriors %>% 
  group_by(probability) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = probability, values_from = value) %>%
  select(-row)

#add codes for positive or negative coefficients
plot.posteriors.wide$effect <- ifelse(
                                  plot.posteriors.wide$lower<0 & 
                                    plot.posteriors.wide$median <0 & 
                                      plot.posteriors.wide$upper <0, 
                                        "negative", ifelse(plot.posteriors.wide$lower>0 
                                                           & plot.posteriors.wide$median >0 & 
                                                             plot.posteriors.wide$upper >0, 
                                                                "positive",
                                                                  "neutral"))


print(plot.posteriors.wide, n = 31)

ggplot(plot.posteriors.wide, aes(x = names, y = median, color = effect)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Covariate") +
  ylab("Causal Effect on Goby Density") +
  coord_flip() 



str(rethinking::extract.samples(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS))

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS %>%
  spread_draws(a_Goby[Zone]) %>%
  head(10)






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









