
# data

set.seed(1908)
N <- 200 # number of pairs
U <- rnorm(N) # simulate confounds
# birth order and family sizes
B1 <- rbinom(N,size=1,prob=0.5) # 50% first borns
M <- rnorm( N , 2*B1 + U )
B2 <- rbinom(N,size=1,prob=0.5)
D <- rnorm( N , 2*B2 + U + 0*M ) # change the 0 to turn on causal influence of mom


# full-luxury bayesian inference
library(rethinking)
library(cmdstanr)
library(lme4)
library(tidyverse)

dat <- list(N=N,M=M,D=D,B1=B1,B2=B2)
set.seed(1908)
flbi <- ulam(
  alist(
    # mom model
    M ~ normal( mu , sigma ),
    mu <- a1 + b*B1 + k*U[i],
    # daughter model
    D ~ normal( nu , tau ),
    nu <- a2 + b*B2 + m*M + k*U[i],
    # B1 and B2
    B1 ~ bernoulli(p),
    B2 ~ bernoulli(p),
    # unmeasured confound
    vector[N]:U ~ normal(0,1),
    # priors
    c(a1,a2,b,m) ~ normal( 0 , 0.5 ),
    c(k,sigma,tau) ~ exponential( 1 ),
    p ~ beta(2,2)
  ), data=dat , chains=4 , cores=4 , iter=2000 , cmdstan=TRUE )


precis(flbi)
#summary(flbi)
traceplot(flbi, ask = FALSE)
stancode(flbi)
plot(flbi)

## end examples ------------

## note code offset using simply "+ log(Area)"

#Goby

library(readr)
goby_master <- read_csv("C:/projects/Goby/data/goby_master.csv")



#goby_master$Since_Breach[is.na(goby_master$Since_Breach)] <- 0
goby_master$Since_Breach <- ifelse(goby_master$Since_Breach == 0.5, 0, 
                                          goby_master$Since_Breach)


dat.temp <- goby_master




## Prep Data
#IVs
Goby   <- dat.temp$Sum_TW  # must be non-negative
Year   <- scale(dat.temp$Year)
SAV    <- scale(dat.temp$SAV)
SB     <- scale(dat.temp$Sum_SB)
SC     <- scale(dat.temp$Sum_SC)
Rain   <- scale(dat.temp$Rain_Sum)
Temp   <- scale(dat.temp$temp_mean)
DO     <- scale(dat.temp$min_DO)
Breach <- scale(dat.temp$Since_Breach)  
Wind   <- scale(dat.temp$u_mean)

Zone   <- as.factor(dat.temp$Zone)
Substrate   <- dat.temp$Dom_substrate
#Random Effects
unique(dat.temp$Zone)
# zone to integer
# Zone   <- as.integer(ifelse(dat.temp$Zone == "E", 1,
#                  ifelse(dat.temp$Zone == "NW", 2,
#                         ifelse(dat.temp$Zone == "W", 3, dat.temp$Zone))))
unique(dat.temp$Dom_substrate)

#pool substrates
Substrate <- as.factor(ifelse(dat.temp$Dom_substrate == "corophium_tubes", "Fine",
                ifelse(dat.temp$Dom_substrate == "mud", "Fine",
                       ifelse(dat.temp$Dom_substrate == "muck", "Fine",
                       ifelse(dat.temp$Dom_substrate == "gravel", "Coarse",
                       ifelse(dat.temp$Dom_substrate == "sand", "Coarse",
                       ifelse(dat.temp$Dom_substrate == "cobble", "Coarse",
                       ifelse(dat.temp$Dom_substrate == "riprap", "Riprap",
                            dat.temp$Dom_substrate))))))))

#Offset
Area <- log(dat.temp$Area)


dat <- data.frame(Goby=Goby, Year=Year, SAV=SAV, SB=SB, SC=SC, Rain=Rain, Temp=Temp, 
            DO=DO, Breach=Breach, Wind=Wind, 
            Zone=Zone, Substrate=Substrate, Area=Area)
#View(dat)
nrow(dat) #n = 349
##remove cases with NAs
dat <- na.omit(dat) #removes 50 cases !  new n = 299, most missingness due to DO
nrow(dat)

#need to convert substrate and zone to numbers
#for now leave out


dat <- tibble(Goby=dat$Goby, 
            Year=dat$Year, 
            SAV=dat$SAV, 
            SB=dat$SB, 
            SC=dat$SC, 
            Rain=dat$Rain, 
            Temp=dat$Temp, 
            DO=dat$DO, 
            Breach=dat$Breach, 
            Wind=dat$Wind, 
            Zone=as.factor(dat$Zone), 
            Substrate=as.factor(dat$Substrate), 
            Area=dat$Area)
dat
summary(dat$Area)
sd(dat$Area)

#test lmer
m1.lmer <- glmer(Goby ~ Year + Rain + SC + SAV + SB + DO + Temp + Breach + Substrate +(1|Zone),
                 data = dat, family = negative.binomial(0.6))
summary(m1.lmer)

#simple model first
set.seed(321)

Goby.m1 <-  ulam(
  alist(
#Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # l  
    a_Goby[Zone] + #random effect  causes chain issues
    beta_Year*Year +
    beta_Rain*Rain + 
    beta_SC*SC + 
    beta_SAV*SAV + 
    beta_SB*SB + 
    beta_DO*DO +
    beta_Temp*Temp +
    beta_Breach*Breach + 
    beta_Substrate*Substrate +
    beta_Wind*Wind +  
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
  
  #SB model
  #SB ~ dgampois( SB_lambda, SB_phi ),
  #SB_lambda <- a_SB + beta_DO*DO + beta_SAV*SAV,
  #SAV model
  
  #SC model
  SC ~ normal( SC_nu, tau),
  SC_nu <- 
    a_SC + 
    beta_DO*DO + 
    beta_SAV*SAV,
  # 
  #Breach model
  Breach ~ normal(Breach_nu, tau),
  Breach_nu <-
    a_Breach +
    beta_Rain*Rain +
    beta_Wind*Wind,
  
  #SAV model
  # SAV ~ normal(SAV_nu, tau),
  # SAV_nu <- 
  #   a_SAV + 
  #   beta_Temp*Temp + 
  #   beta_Rain*Rain +
  #   beta_Wind*Wind,
  
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
    phi ~ normal( 0 , 5 ),
    c(SC_phi, SB_phi) ~ dexp(100) 
), 
data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE)

#beepr::beep(0)

precis(Goby.m1, depth = 2)
plot(Goby.m1)
#plot(Goby.m1, depth = 2)
summary(Goby.m1)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m1)

save(Goby.m1, file = "Output/Goby.m1.RData")

set.seed(100)
### model with no network
Goby.m1.no.network <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  
      #beta_Year*Year +
      beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      beta_Temp*Temp +
      beta_Breach*Breach + 
      #beta_Substrate*Substrate +
      beta_Wind*Wind +  
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
      beta_Wind
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ normal( 0 , 5 ),
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=2000 , cmdstan=TRUE)

precis(Goby.m1.no.network, depth = 2)
plot(Goby.m1.no.network)
#plot(Goby.m1, depth = 2)
summary(Goby.m1.no.network)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m1.no.network)


save(Goby.m1.no.network, file = "Output/Goby.m1.no.network.RData")

par(mfrow = c(1,2))



## remove Goby ~ Rain + wind + breach since not directly tied to Goby in DAG

Goby.m2 <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  causes chain issues
      beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      beta_Temp*Temp +
      #beta_Breach*Breach + 
      #beta_Substrate*Substrate +
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
    
    #SB model
    SB ~ normal( SB_nu, tau ),
    SB_nu <- a_SB + beta_DO*DO + beta_SAV*SAV,
    #SAV model
    
    #SC model
    SC ~ normal( SC_nu, tau),
    SC_nu <- 
      a_SC + 
      beta_DO*DO + 
      beta_SAV*SAV,
    # 
    #Breach model
    Breach ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_Breach +
      beta_Rain*Rain +
      beta_Wind*Wind,
    
    #SAV model
    # SAV ~ normal(SAV_nu, tau),
    # SAV_nu <- 
    #   a_SAV + 
    #   beta_Temp*Temp + 
    #   beta_Rain*Rain +
    #   beta_Wind*Wind,
    
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
      beta_Wind#,
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ normal( 0 , 5 ),
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE)

#beepr::beep(0)

precis(Goby.m2, depth = 2)

plot(Goby.m2)
#plot(Goby.m1, depth = 2)
summary(Goby.m2)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2)

save(Goby.m2, file = "Output/Goby.m2.RData")


### No year to just get causes
Goby.m2.no.year <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  causes chain issues
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      beta_Temp*Temp +
      #beta_Breach*Breach + 
      #beta_Substrate*Substrate +
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
    SB_nu <- a_SB + beta_DO*DO + beta_SAV*SAV,
    #SAV model
    
    #SC model
    SC ~ normal( SC_nu, tau),
    SC_nu <- 
      a_SC + 
      beta_DO*DO + 
      beta_SAV*SAV,
    # 
    #Breach model
    Breach ~ normal(Breach_nu, tau),
    Breach_nu <-
      a_Breach +
      beta_Rain*Rain +
      beta_Wind*Wind,
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
        a_SAV + 
       beta_Temp*Temp + 
       beta_Rain*Rain +
       beta_Wind*Wind,
    
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
      beta_Wind#,
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm(1,10), #to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE)
#beepr::beep(0)

precis(Goby.m2.no.year, depth = 2)

plot(Goby.m2.no.year)
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2.no.year)


### refining to get closer to Dag
### No year to just get causes
Goby.m2.no.year.DAG <-  ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp() constrain to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  causes chain issues
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
      #beta_Substrate*Substrate +
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
      beta_Rain*Rain,# +
      #beta_Wind*Wind, Ask Darren if include wind?
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Temp*Temp +
      beta_DO*DO,
      
    
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
      beta_Wind#,
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm(1,10), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE
)

precis(Goby.m2.no.year.DAG, depth = 2)
plot(Goby.m2.no.year.DAG)
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year.DAG)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2.no.year.DAG)

save(Goby.m2.no.year.DAG, file = "Output/Goby.m2.no.year.DAG.RData")


plot(Goby.m1.no.network)
plot(Goby.m2.no.year.DAG)

### Simulate Goby with intervention of Wind * Temp --> Goby------------------
# set Wind to -0.25

post <- extract.samples(Goby.m2.no.year.DAG)

#Effect of Rain on Breach
quantile( with(post, beta_Rain*beta_Breach)) ##-0.01
#effect of breach on DO
quantile( with(post, beta_DO*beta_Breach)) ## zero
#effect of Wind on DO
quantile( with(post, beta_DO*beta_Wind))  ## zero



##### test ways to avoid negative response values--------------
### refining to get closer to Dag
### No year to just get causes
Goby.m2.no.year.DAG.no_negatives <- ulam(
  alist(
    #Goby model
    Goby ~ dgampois( mu, exp(phi)), # exp(log_scale) to constrain scale to positive reals
    log(mu) <-  # log(mu) from dgampois help
      a_Goby[Zone] + #random effect  causes chain issues
      #beta_Year*Year +
      #beta_Rain*Rain + 
      beta_SC*SC + 
      beta_SAV*SAV + 
      beta_SB*SB + 
      beta_DO*DO +
      #beta_Temp*Temp +
      #beta_Breach*Breach + 
      #beta_Substrate*Substrate +
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
      beta_Rain*Rain,# +
    #beta_Wind*Wind, Ask Darren if include wind?
    
    #SAV model
    SAV ~ normal(SAV_nu, tau),
    SAV_nu <- 
      a_SAV + 
      beta_Temp*Temp +
      beta_DO*DO,
    
    
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
      beta_Wind#,
      #beta_Zone,
      #beta_Substrate  
    ) ~ normal( 0 , 0.5 ),
    tau ~ exponential(1),
    phi ~ dnorm(1,10), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dexp(100) 
  ), 
  data=dat , chains=4 , cores=4 , iter=500 , cmdstan=TRUE
)

precis(Goby.m2.no.year.DAG.no_negatives, depth = 2)
plot(Goby.m2.no.year.DAG.no_negatives)
#plot(Goby.m1, depth = 2)
summary(Goby.m2.no.year.DAG.no_negatives)
#traceplot(Goby.m1, ask = FALSE)
stancode(Goby.m2.no.year.DAG.no_negatives)

save(Goby.m2.no.year.DAG.no_negatives, file = "Output/Goby.m2.no.year.DAG.no_negatives.RData")






