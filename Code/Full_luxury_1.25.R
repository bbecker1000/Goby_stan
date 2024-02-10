

## end examples ------------

## note code offset using simply "+ log(Area)"

#Goby

library(readxl)
library(readr)
library(tidyverse)
library(lme4)
library(rethinking)
library(rstan)
library(rstanarm)
library(sjPlot)
library(marginaleffects)
goby_master <- read_csv("C:/projects/Goby/data/goby_master.csv")

breach_days <- read_excel("C:/projects/Goby/data/RodeoLagoon-Status_WY1995-present.xlsx", 
                                 col_types = c("date", "numeric", "numeric", "text", "text"))
#keep first three columns
breach_days <- breach_days[,c(1:3)]
#only keep value == 1
breach_days_sum <- breach_days %>%
  filter(RodeoLagoonMouth == 1) %>%
  group_by(WaterYear) %>%
  summarize(BreachDays = sum(RodeoLagoonMouth))

#find missing years
#2010, 2014 2020 2021  = 0 breaches
ZeroYears <- tibble(WaterYear = c(2010, 2014, 2020, 2021), BreachDays = c(0,0,0,0))

breach_days_sum2 <- rbind(breach_days_sum, ZeroYears)
breach_days_sum3 <- breach_days_sum2 %>% arrange(WaterYear)
breach_days_sum <- breach_days_sum3

#key for join
breach_days_sum$Year <- breach_days_sum$WaterYear

#Join Breach Days to goby master 
# WY is previous 12 months ...WY1995 = 1994-10-01 to 1995-09-30
# so match WY to goby sampling year...
goby_masterB <- goby_master %>% left_join(breach_days_sum, by = "Year")
goby_master <- goby_masterB

#goby_master$Since_Breach[is.na(goby_master$Since_Breach)] <- 0
goby_master$Since_Breach <- ifelse(goby_master$Since_Breach == 0.5, 1, 
                                          goby_master$Since_Breach)

#make since breach a 0/1
#wonky categorical compared to rstanarm, so need to fix later
#goby_master$Since_Breach <- ifelse(goby_master$Since_Breach ==  0, 0, 1)
                       
goby_master$Since_Breach

hist(goby_master$micro_sum)
summary(goby_master$micro_sum)

dat.temp <- goby_master




## Prep Data
#IVs
Goby   <- dat.temp$Sum_TW  # must be non-negative
Year   <- scale(dat.temp$Year)
Year_int <- as.integer(dat.temp$Year-1995)
SAV    <- scale(dat.temp$SAV)
SB     <- scale(dat.temp$Sum_SB)
SB_count     <- dat.temp$Sum_SB # if using counts
SC     <- scale(dat.temp$Sum_SC)
SC_count     <- dat.temp$Sum_SC
Rain   <- scale(dat.temp$Rain_Sum)
Temp   <- scale(dat.temp$temp_mean)
Temp_2   <- scale(dat.temp$temp_mean^2)
DO     <- scale(dat.temp$min_DO)
Micro  <- scale(dat.temp$micro_sum) #added 2024-01-23, poor fit if raw counts
Breach <- scale(dat.temp$Since_Breach)
Breach_count <- dat.temp$Since_Breach
BreachDays <- scale(dat.temp$BreachDays)
BreachDays_Count <- (dat.temp$BreachDays)
BreachDays_2 <- scale(dat.temp$BreachDays^2)
BreachDays_Count_2 <- (dat.temp$BreachDays^2)
#Breach <- as.factor(dat.temp$Since_Breach)  # need to fix this, but categorical is wonky
Wind   <- scale(dat.temp$u_mean)





Zone   <- as.factor(dat.temp$Zone)
Substrate   <- as.factor(dat.temp$Dom_substrate)
#Random Effects
unique(dat.temp$Zone)
# zone to integer, keep as three zones
Zone   <- as.integer(ifelse(dat.temp$Zone == "E", 0,
                  ifelse(dat.temp$Zone == "NW", 1,
                         ifelse(dat.temp$Zone == "W", 2, dat.temp$Zone))))
unique(dat.temp$Dom_substrate)

plot(dat.temp$Sum_TW, Zone)

#pool substrates
Substrate <- as.integer(ifelse(dat.temp$Dom_substrate == "corophium_tubes", "1",
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
                  SC_count=SC_count, Rain=Rain, Temp=Temp, Temp_2=Temp_2, #added 2024-02-09
            DO=DO, Breach=Breach,  Breach_count=Breach_count, #added 2024-01-23
            BreachDays=BreachDays, BreachDays_2=BreachDays_2,
            BreachDays_Count=BreachDays_Count, BreachDays_Count_2=BreachDays_Count_2,
            Wind=Wind, Micro=Micro,
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
              Temp_2=dat$Temp_2, 
              DO=dat$DO, 
              Breach=dat$Breach, 
              Breach_count=dat$Breach_count, 
              BreachDays=dat$BreachDays,
              BreachDays_2=dat$BreachDays_2,
              BreachDays_Count=dat$BreachDays_Count,
              BreachDays_Count_2=dat$BreachDays_Count_2,
              Wind=dat$Wind, 
              Micro=dat$Micro, #added 2024-01-23
              Zone=as.factor(dat$Zone), 
              Substrate=as.factor(dat$Substrate), 
              Area=dat$Area)
dat.missing

#remove variables with NAs except Temp, DO, and SAV ok
dat.missing <- dat.missing %>%
  filter_at(vars(Year, SAV, SB_count, Goby, SC_count, SB, SC, Substrate, Wind, Zone, Micro, Area), 
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
            Temp_2=dat$Temp_2,
            DO=dat$DO, 
            Breach=dat$Breach, 
            Breach_count=dat$Breach_count, 
            BreachDays = dat$BreachDays,
            BreachDays_2 = dat$BreachDays_2,
            BreachDays_Count = dat$BreachDays_Count,
            BreachDays_Count_2 = dat$BreachDays_Count_2,
            Wind=dat$Wind,
            Micro=dat$Micro,
            Zone=as.factor(dat$Zone), 
            Substrate=as.factor(dat$Substrate), 
            Area=dat$Area)
dat
summary(dat$Area)
sd(dat$Area)


ggplot(dat, aes(x = BreachDays_Count, y = Goby/Area, color = Zone)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~Zone)

ggplot(dat, aes(x = Micro, y = Goby/Area)) +
  geom_point() +
  geom_smooth(method = "loess")





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
  data=dat , chains=4 , cores=4 , iter=2000 , cmdstan=TRUE
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
#### add year for trend.
t0 <- Sys.time()

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS <-  ulam(
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
      beta_Temp*Temp +    # added 2024-02-09
      beta_Temp_2*Temp_2 +  # added 2024-02-09
      beta_Zone*Zone + 
      beta_ZW*Wind*Zone + # interaction between wind and Zone
      #slope of wind piling up algae
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
      # beta_BreachDays*BreachDays + # does breach effect DO?  
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
      beta_Temp_2, 
      beta_Micro,
      #beta_Breach,
      beta_BreachDays,
      beta_BreachDays_2,
      beta_Wind,
      beta_Zone,
      beta_ZW,
      beta_Substrate
    ) ~ normal( 0 , 1 ),
    tau ~ exponential(1),
    phi ~ dnorm( 1, 10 ), #from dgampois help to keep from going negative
    c(SC_phi, SB_phi) ~ dnorm(1, 10)  # use dexp(100) if not neg.bin
  ), 
  data=dat , chains=4 , cores=parallel::detectCores() , iter=20000 , 
  cmdstan=FALSE # to get stanfit object
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
              "beta_Wind",
              "beta_Zone",
              "beta_ZW"),
     xlab = "Beta Coefficient", 
     main = "network model-Breach Direct")


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


Year<-as_tibble(quantile( with(post, beta_Year), probs = PROBS)) ## ns 
Substrate_SC<-as_tibble(quantile( with(post, beta_Substrate*beta_SC_count), probs = PROBS)) ## negative
Rain_DO_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_DO*beta_SAV*beta_SB_count), probs = PROBS))

Breach<-as_tibble(quantile( with(post, beta_BreachDays), probs = PROBS)) #positive
Zone<-as_tibble(quantile( with(post, beta_Zone), probs = PROBS))   #positive
Zone_Wind_int<-as_tibble(quantile( with(post, beta_ZW), probs = PROBS))  #positive
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
names<- c("Rain_Breach_DO", "Breach_DO", "Wind_DO", "Wind_Temp_DO", "Wind_DO_SAV", "Wind_DO_SC",
          "DO_SB", "Year", "Substrate_SC", "Rain_DO_SAV_SB", "Breach", "Zone", "Zone_Wind_int", 
          "Micro", "SAV_SB", "SAV_SC", "SC", "SB","Breach_DO_SC", "Rain_Breach_Temp",
          "Breach_Temp",
          "Wind_Temp",
          "Wind_Temp_SAV",
          "Wind_Temp_SC",
          "Temp_SB",
          "Rain_Temp_SAV_SB")
#add probabilities
plot.posteriors<-rbind(Rain_Breach_DO, Breach_DO, Wind_DO, Wind_Temp_DO, Wind_DO_SAV, Wind_DO_SC,
        DO_SB, Year, Substrate_SC, Rain_DO_SAV_SB, Breach, Zone, Zone_Wind_int,
        Micro, SAV_SB, SAV_SC, SC, SB, Breach_DO_SC, Rain_Breach_Temp,
        Breach_Temp,
        Wind_Temp,
        Wind_Temp_SAV,
        Wind_Temp_SC,
        Temp_SB,
        Rain_Temp_SAV_SB)
#add names
plot.posteriors$names <- rep(names, each=3)
#add probabilities names
plot.posteriors$probability <- rep(c("lower", "median", "upper"), times = length(names))

plot.posteriors.wide <- plot.posteriors %>% 
                            pivot_wider(names_from = c(probability), values_from = value)
#add codes for positive or negative coefficients
plot.posteriors.wide$effect <- ifelse(plot.posteriors.wide$lower<0 & plot.posteriors.wide$median <0 & plot.posteriors.wide$upper <0, "negative",
                                                                  ifelse(plot.posteriors.wide$lower>0 & plot.posteriors.wide$median >0 & plot.posteriors.wide$upper >0, "positive",
                                                                  "neutral"))


print(plot.posteriors.wide, n = 27)

ggplot(plot.posteriors.wide, aes(x = names, y = median, color = effect)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() 



str(rethinking::extract.samples(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct))

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct %>%
  spread_draws(a_Goby[Zone]) %>%
  head(10)



#effects plot
#k <- PSIS(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, pointwise = TRUE)$k
plot(dat$BreachDays, dat$Goby/dat$Area)
ns <- 100
P_seq <- seq( from = -1.2, to = 2.3, length.out = ns)
lambda <- link( Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, data=data.frame(P=P_seq, cid=1))
lmu <- apply( lambda, 2, mean)
lci <- apply( lambda, 2, PI)
lines(P_seq, lmu, lty = 2, lwd = 1.5)
shade(lci, P_seq, xpd = TRUE)



#compare plots with and without direct breach effects()

par(mfrow = c(1,1))
plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_BreachDays", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient", 
     main = "network model")

plot(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, 
     pars = c("beta_Year", "beta_Substrate", "beta_Wind", "beta_BreachDays", "beta_Temp",
              "beta_DO", "beta_SB_count", "beta_SC_count", "beta_SAV", "beta_Rain"),
     xlab = "Beta Coefficient", 
     main = "network model - Breach Direct")



#plot(Goby.m1, depth = 2)
summary(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct)
#traceplot(Goby.m2.no.year.DAG.SC.SB_counts.Breach.pois, ask = FALSE)
stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct)


save(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RData")
#load(file = "Output/Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RData")






ggplot(dat, aes(Year_int+1995, y = (Goby/Area), color = as.factor(Zone))) + 
  geom_point() +
  geom_smooth() +
  ylab("density/log(m2)") +
  xlab("Year") +
  theme_classic(base_size=22)


## corrected effects of covariates on Goby

post <- extract.samples(Goby.m2.year.DAG.SC.SB_counts.BreachDays)

#Effect of Rain --> Breach --> DO --> Goby
quantile( with(post, beta_Rain*beta_BreachDays*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Breach --> DO --> Goby
quantile( with(post, beta_BreachDays*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Wind --> Breach --> DO --> Goby
quantile( with(post, beta_Wind*beta_BreachDays*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Wind --> DO --> Goby
quantile( with(post, beta_Wind*beta_DO), probs = c(0.1, 0.5, 0.9)) ## negative   
#Effect of Wind --> Temp --> DO --> Goby
quantile( with(post, beta_Wind*beta_Temp*beta_DO), probs = c(0.1, 0.5, 0.9)) ## negative  
#Effect of Wind --> Temp --> SAV --> Goby
quantile( with(post, beta_Wind*beta_SAV*beta_DO), probs = c(0.1, 0.5, 0.9)) ## negative  
#Effect of Wind --> DO --> SB --> Goby
quantile( with(post, beta_Wind*beta_DO*beta_SB_count), probs = c(0.1, 0.5, 0.9)) ## ns 
#Effect of DO --> SB --> Goby
quantile( with(post, beta_DO*beta_SB_count), probs = c(0.1, 0.5, 0.9)) ## ns 




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
      a_Goby[Substrate] + #random effect 
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



#Causal Effects calcs:
post <- extract.samples(Goby.m2.year.DAG.SC.SB_counts.Breach.pois.missing)

#Effect of Rain --> Breach --> DO --> Goby
quantile( with(post, beta_Rain*beta_Breach_count*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Breach --> DO --> Goby
quantile( with(post, beta_Breach_count*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Wind --> Breach --> DO --> Goby
quantile( with(post, beta_Wind*beta_Breach_count*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns
#Effect of Wind --> DO --> Goby
quantile( with(post, beta_Wind*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns   
#Effect of Wind --> Temp --> DO --> Goby
quantile( with(post, beta_Wind*beta_Temp*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns  
#Effect of Wind --> Temp --> SAV --> Goby
quantile( with(post, beta_Wind*beta_SAV*beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns  
#Effect of Wind --> DO --> SB --> Goby
quantile( with(post, beta_Wind*beta_DO*beta_SB_count), probs = c(0.1, 0.5, 0.9)) ## ns 
#Effect of DO --> SB --> Goby
quantile( with(post, beta_DO*beta_SB_count), probs = c(0.1, 0.5, 0.9)) ## ns 

quantile( with(post, beta_DO), probs = c(0.1, 0.5, 0.9)) ## ns







