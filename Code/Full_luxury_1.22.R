

## end examples ------------

## note code offset using simply "+ log(Area)"

#Goby

library(readr)
library(tidyverse)
library(lme4)
library(rethinking)
library(rstan)
library(rstanarm)
library(sjPlot)
library(marginaleffects)
goby_master <- read_csv("C:/projects/Goby/data/goby_master.csv")



#goby_master$Since_Breach[is.na(goby_master$Since_Breach)] <- 0
goby_master$Since_Breach <- ifelse(goby_master$Since_Breach == 0.5, 1, 
                                          goby_master$Since_Breach)

#make since breach a 0/1
#wonky categorical compared to rstanarm, so need to fix later
#goby_master$Since_Breach <- ifelse(goby_master$Since_Breach ==  0, 0, 1)
                       
goby_master$Since_Breach


dat.temp <- goby_master




## Prep Data
#IVs
Goby   <- dat.temp$Sum_TW  # must be non-negative
Year   <- scale(dat.temp$Year)
Year_int <- (dat.temp$Year-1995)
SAV    <- scale(dat.temp$SAV)
SB     <- scale(dat.temp$Sum_SB)
SB_count     <- dat.temp$Sum_SB # if using counts
SC     <- scale(dat.temp$Sum_SC)
SC_count     <- dat.temp$Sum_SC
Rain   <- scale(dat.temp$Rain_Sum)
Temp   <- scale(dat.temp$temp_mean)
DO     <- scale(dat.temp$min_DO)
Breach <- scale(dat.temp$Since_Breach)  
Breach_count <- dat.temp$Since_Breach  
#Breach <- as.factor(dat.temp$Since_Breach)  # need to fix this, but categorical is wonky
Wind   <- scale(dat.temp$u_mean)

Zone   <- as.factor(dat.temp$Zone)
Substrate   <- as.factor(dat.temp$Dom_substrate)
#Random Effects
unique(dat.temp$Zone)
# zone to integer, pool NW and W
Zone   <- as.integer(ifelse(dat.temp$Zone == "E", 0,
                  ifelse(dat.temp$Zone == "NW", 1,
                         ifelse(dat.temp$Zone == "W", 2, dat.temp$Zone))))
unique(dat.temp$Dom_substrate)

plot(dat.temp$Sum_TW, Zone)

#pool substrates
Substrate <- as.integer(ifelse(dat.temp$Dom_substrate == "corophium_tubes", "0",
                ifelse(dat.temp$Dom_substrate == "mud", "0",
                       ifelse(dat.temp$Dom_substrate == "muck", "0",
                       ifelse(dat.temp$Dom_substrate == "gravel", "1",
                       ifelse(dat.temp$Dom_substrate == "sand", "1",
                       ifelse(dat.temp$Dom_substrate == "cobble", "1",
                       ifelse(dat.temp$Dom_substrate == "riprap", "1",
                            dat.temp$Dom_substrate))))))))

#Offset
Area <- log(dat.temp$Area)


dat <- data.frame(Goby=Goby, Year=Year, Year_int = Year_int, SAV=SAV, SB=SB, SB_count=SB_count, SC=SC, 
                  SC_count=SC_count, Rain=Rain, Temp=Temp, 
            DO=DO, Breach=Breach,  Breach_count=Breach_count, Wind=Wind, 
            Zone=Zone, Substrate=Substrate, Area=Area)


#View(dat)
nrow(dat) #n = 349

#save a file with all cases including missing cases
dat.missing <- tibble(Goby=dat$Goby, 
              Year=dat$Year, 
              Year_int=dat$Year_int, 
              SAV=dat$SAV, 
              SB=dat$SB,
              SB_count=dat$SB_count,
              SC=dat$SC, 
              SC_count=dat$SC_count, 
              Rain=dat$Rain, 
              Temp=dat$Temp, 
              DO=dat$DO, 
              Breach=dat$Breach, 
              Breach_count=dat$Breach_count, 
              Wind=dat$Wind, 
              Zone=as.factor(dat$Zone), 
              Substrate=as.factor(dat$Substrate), 
              Area=dat$Area)
dat.missing

#remove variables with NAs except Temp, DO, and SAV ok
dat.missing <- dat.missing %>%
  filter_at(vars(Year, SAV, SB_count, Goby, SC_count, SB, SC, Substrate, Wind, Zone, Area), 
            all_vars(!is.na(.)))

##remove cases with NAs
dat <- na.omit(dat) #removes 50 cases !  new n = 299, most missingness due to DO
nrow(dat)




dat.list <- as.list(dat)  # for stan


#need to convert substrate and zone to numbers
#for now leave out


dat <- tibble(Goby=dat$Goby, 
            Year=dat$Year, 
            Year_int=dat$Year_int, 
            SAV=dat$SAV, 
            SB=dat$SB,
            SB_count=dat$SB_count,
            SC=dat$SC, 
            SC_count=dat$SC_count, 
            Rain=dat$Rain, 
            Temp=dat$Temp, 
            DO=dat$DO, 
            Breach=dat$Breach, 
            Breach_count=dat$Breach_count, 
            Wind=dat$Wind, 
            Zone=as.factor(dat$Zone), 
            Substrate=as.factor(dat$Substrate), 
            Area=dat$Area)
dat
summary(dat$Area)
sd(dat$Area)

#test lmer
m1.lmer <- glmer(Goby ~ Year + 
                   Rain + I(SC_count/Area) + SAV + I(SB_count/Area) + DO + Temp + 
                Breach_count + Wind + 
                Substrate +(1|Zone), data = dat, family = negative.binomial(0.6), offset = Area)
summary(m1.lmer)

library(rstanarm)
m1.stan_glm <- stan_glmer(Goby ~ Year 
                          Rain + SC_count + SAV + SB + DO + Temp + Breach + 
                   Substrate + (1|Zone), data = dat, family = neg_binomial_2(),
                   chains = 4, cores = 4, iter = 1000)
beepr::beep(0)
summary(m1.stan_glm, digits = 3)

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
      beta_Breach*Breach + 
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC,
      beta_SAV, 
      beta_SB, 
      beta_DO,
      beta_Temp, 
      beta_Breach,
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
plot(Goby.m1.no.network)
#plot(Goby.m1, depth = 2)
summary(Goby.m1.no.network)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m1.no.network)

# "no network" ulam coefficients are very close to stan_glmer coefficients - Yay!
# the Zone random effect intercepts are different, but have similar intervals 

save(Goby.m1.no.network, file = "Output/Goby.m1.no.network.RData")

par(mfrow = c(1,1))



#remove Goby ~ Rain + wind + breach since not directly tied to Goby in DAG

### refining to get closer to Dag
### No year to just get causes
Goby.m2.with.year.DAG <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Breach*Breach + 
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_Breach*Breach + 
      beta_Wind*Wind,
    
    #SB model
    SB ~ normal( SB_nu, tau ),
    SB_nu <- a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #SC model
    SC ~ normal( SC_nu, tau),
    SC_nu <- 
      a_SC + 
      beta_DO*DO + 
      beta_SAV*SAV,
    
    #Breach model
    Breach ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_Breach +
      beta_Rain*Rain +
      beta_Wind*Wind, #Ask Darren if include wind?
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain,
    
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC,
      beta_SAV, 
      beta_SB, 
      beta_DO,
      beta_Temp, 
      beta_Breach,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE
)

beepr::beep(0)

precis(Goby.m2.with.year.DAG, depth = 2)
plot(Goby.m2.with.year.DAG)
summary(Goby.m2.with.year.DAG)
stancode(Goby.m2.with.year.DAG)

save(Goby.m2.with.year.DAG, file = "Output/Goby.m2.with.year.DAG.RData")


### refining to get closer to Dag
### No year to just get causes
Goby.m2.no.year.DAG <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
      beta_Substrate*Substrate +
      #beta_Wind*Wind +  
      #beta_Zone*Zone + 
      Area, #offset already logged
    
    #Zone RE model
    a_Goby[Zone] ~ normal(mu_Zone, tau_Zone),
    mu_Zone ~ normal(0, 5),
    tau_Zone ~ exponential(1),
    
    #DO model
    DO ~ normal( DO_nu , tau ),
    DO_nu <- 
      a_DO + 
      beta_Temp*Temp +
      beta_Breach*Breach + 
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_Breach*Breach + 
      beta_Wind*Wind,
    
    #SB model
    SB ~ normal( SB_nu, tau ),
    SB_nu <- a_SB + 
      beta_DO*DO + 
      beta_SAV*SAV +
      Area, #offset

    #SC model
    SC ~ normal( SC_nu, tau),
    SC_nu <- 
      a_SC + 
      beta_DO*DO + 
      beta_SAV*SAV +
      Area, #offset
    
    #Breach model
    Breach ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_Breach +
      beta_Rain*Rain +
      beta_Wind*Wind, #Ask Darren if include wind?
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain,
      
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      #beta_Year, 
      beta_Rain,
      beta_SC,
      beta_SAV, 
      beta_SB, 
      beta_DO,
      beta_Temp, 
      beta_Breach,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE
)

beepr::beep(0)

precis(Goby.m2.no.year.DAG, depth = 2)
plot(Goby.m2.no.year.DAG)
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year.DAG)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2.no.year.DAG)

save(Goby.m2.no.year.DAG, file = "Output/Goby.m2.no.year.DAG.RData")


plot(Goby.m1.no.network)
plot(Goby.m2.no.year.DAG)
plot(Goby.m2.with.year.DAG)

### Simulate Goby with intervention of Wind * Temp --> Goby------------------
# set Wind to -0.25

post <- extract.samples(Goby.m2.no.year.DAG)

#Effect of Rain on Breach
quantile( with(post, beta_Rain*beta_Breach)) ##-0.01
#effect of breach on DO
quantile( with(post, beta_DO*beta_Breach)) ## zero
#effect of Wind on DO
quantile( with(post, beta_SAV*beta_SC))  ## zero



#### SB as counts instead of scaled
#### Modify dat b/c Need to use SB counts!
### Trial runs had all divergences...not sure why?




#####------------try for breach = 0,1,2,
##              SB and SC = neg.bin, but made coefficients really tiny

### refining to get closer to Dag
### No year to just get causes


Goby.m2.no.year.DAG.SC.SB_counts <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      #a_Goby[Substrate] + #random effect 
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + 
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
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
      beta_Breach*Breach + 
      beta_Wind*Wind,
    
    #Temp model
    Temp ~ normal( Temp_nu , tau ),
    Temp_nu <- 
      a_Temp + 
      beta_Breach*Breach + 
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
    # Breach ~ poisson(Breach_nu),
    # log(Breach_nu) <-
    #   a_Breach +
    #   beta_Rain*Rain +
    #   beta_Wind*Wind, 
    
    
    #Breach model scaled
    Breach ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_Breach +
      beta_Rain*Rain +
      beta_Wind*Wind, 
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain,
    
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      #beta_Year, 
      beta_Rain,
      #beta_SC,
      beta_SC_count,
      beta_SAV, 
      #beta_SB,
      beta_SB_count,
      beta_DO,
      beta_Temp, 
      beta_Breach,
      #beta_Breach_count,
      beta_Wind,
      #beta_Zone,
      beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE
)

beepr::beep(0)

precis(Goby.m2.no.year.DAG.SC.SB_counts, depth = 2)
plot(Goby.m2.no.year.DAG.SC.SB_counts)
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year.DAG.SC.SB_counts)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2.no.year.DAG.SC.SB_counts)




##### - SB and SC neg-bin, breach = poisson

### refining to get closer to Dag
### No year to just get causes
Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), 
    log(mu) <- 
      a_Goby[Zone] + #random effect 
      #a_Goby[Substrate] + #random effect 
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC_count*SC_count + #updated
      beta_SAV*SAV + 
      beta_SB_count*SB_count + #competition
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
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
      beta_Wind*Wind +
      beta_SAV*SAV,
    
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
      beta_Rain*Rain +
      beta_Wind*Wind + 
      Area,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain +
      beta_Temp*Temp,
    
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      #beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
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
  data=dat , chains=4 , cores=4 , iter=1000 , cmdstan=TRUE
)

beepr::beep(0)


precis(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, depth = 2)
plot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, 
     pars = c("beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"))
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois)
par(mfrow = c(1,2))

### Simulate Goby with intervention of Wind * Temp --> Goby------------------
# set Wind to -0.25

post <- extract.samples(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois)

#Effect of Rain on Breach
quantile( with(post, beta_Rain*beta_Breach_count), probs = c(0.1, 0.5, 0.9)) ## ns
#effect of breach on DO
quantile( with(post, beta_DO*beta_Breach_count), probs = c(0.1, 0.5, 0.9)) ## ns
#effect of SAV on SC
quantile( with(post, beta_SAV*beta_SC_count), probs = c(0.1, 0.5, 0.9))  ## neg

quantile( with(post, beta_Temp*beta_SAV), probs = c(0.1, 0.5, 0.9))  ## positive

quantile( with(post, beta_DO*beta_Temp), probs = c(0.1, 0.5, 0.9)) ## neg

quantile( with(post, beta_Wind, beta_DO), probs = c(0.1, 0.5, 0.9)) ## positive





#### add year for trend.
Goby.m2.year.DAG.SC.SB_counts.Breach.pois <-  ulam(
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
      #beta_Breach*Breach + 
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
      beta_Wind*Wind +
      beta_SAV*SAV,
    
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
      beta_Rain*Rain +
      beta_Wind*Wind,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain +
      beta_Temp*Temp,
    
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
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
  data=dat , chains=4 , cores=4 , iter=1000 , cmdstan=TRUE
)

beepr::beep(0)


precis(Goby.m2.year.DAG.SC.SB_counts.Breach.pois)
plot(Goby.m2.year.DAG.SC.SB_counts.Breach.pois, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient")

#plot(Goby.m1, depth = 2)
summary(Goby.m2.year.DAG.SC.SB_counts.Breach.pois)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.year.DAG.SC.SB_counts.Breach.pois)


save(Goby.m2.year.DAG.SC.SB_counts.Breach.pois, file = "Output/Goby.m2.year.DAG.SC.SB_counts.Breach.pois.RData")

ggplot(dat, aes(Year_int+1995, y = (Goby/Area), color = as.factor(Zone))) + 
  geom_point() +
  geom_smooth() +
  ylab("density/log(m2)") +
  xlab("Year") +
  theme_classic(base_size=22)




####----------
#### try with imputed data
### need priors for data and the n = 349 data list 

### add year for trend.
Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing <-  ulam(
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
      #beta_Breach*Breach + 
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
      beta_Wind*Wind +
      beta_SAV*SAV,
    
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
      beta_Rain*Rain +
      beta_Wind*Wind,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Rain*Rain +
      beta_Temp*Temp,
    
    #missing data priors
    c(Temp, DO) ~ normal (0,1),
    
    #fixed effects priors
    c(a_Goby, a_Breach, a_DO, a_SB, a_SC, a_SAV, a_Temp,
      beta_Year, 
      beta_Rain,
      beta_SC_count,
      beta_SAV, 
      beta_SB_count,
      beta_DO,
      beta_Temp, 
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
  data=dat.missing , chains=4 , cores=4 , iter=1000 , cmdstan=TRUE
)

beepr::beep(0)


precis(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing, depth = 2)
plot(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient")

#plot(Goby.m1, depth = 2)
summary(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing)

save(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing, 
     file = "Output/Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing.RData")


par(mfrow = c(1,2))

plot(Goby.m2.year.DAG.SC.SB_counts.Breach.pois, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient")

plot(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_Breach_count", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient")


