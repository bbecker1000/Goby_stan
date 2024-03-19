# Final Models
source("code/1_DataPrep.R")






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
      beta_SAV_2*SAV_2 + 
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
      beta_Temp_2,
      beta_Micro,
      beta_Wind
      # beta_ZW
      ) ~ 
      normal( 0 , 0.75 ),

      beta_Rain         ~ normal( 0.25 , 0.25 ),  # more rain = more goby
      beta_SC_count     ~ normal(-0.10 , 0.25 ),  # sculpins eat goby larvae
      beta_SAV          ~ normal( 0.00 , 0.25 ),  # goldilocks 
      beta_SAV_2          ~ normal( 0.00 , 0.25 ),  # goldilocks 
      beta_DO           ~ normal( 0.25 , 0.25 ),  # more DO good
      beta_BreachDays   ~ normal( 0.25 , 0.25 ),  # more breach good
      beta_BreachDays_2 ~ normal( 0.25 , 0.25 ),  # goldilocks
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


stancode(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS)

fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS

plot(fit, depth = 2)
trankplot(fit)



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

DO<-as_tibble(quantile( with(post, beta_DO), probs = PROBS)) ## ns 

Rain_Breach_Temp<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays*beta_Temp), probs = PROBS)) ## NS
#Effect of Breach --> Temp --> Goby
Rain_Breach_2<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays_2), probs = PROBS)) ## NS
#Effect of Rain --> Breach --> Goby

Rain_Breach<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays), probs = PROBS)) ## NS
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
Year_2<-as_tibble(quantile( with(post, beta_Year_2), probs = PROBS)) ## ns 


Substrate_SC<-as_tibble(quantile( with(post, beta_Substrate*beta_SC_count), probs = PROBS)) ## negative
Substrate<-as_tibble(quantile( with(post, beta_Substrate), probs = PROBS))
Rain_DO_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_DO*beta_SAV*beta_SB_count), probs = PROBS))

Breach<-as_tibble(quantile( with(post, beta_BreachDays), probs = PROBS)) #positive
Rain_Micro<-as_tibble(quantile( with(post, beta_Rain*beta_Micro), probs = PROBS))
Micro<-as_tibble(quantile( with(post, beta_Micro), probs = PROBS))
SAV_SB<-as_tibble(quantile( with(post, beta_SAV*beta_SB_count), probs = PROBS))
SAV_SC<-as_tibble(quantile( with(post, beta_SAV*beta_SC_count), probs = PROBS))
SAV<-as_tibble(quantile( with(post, beta_SAV), probs = PROBS))
SAV_2<-as_tibble(quantile( with(post, beta_SAV_2), probs = PROBS))
SC<-as_tibble(quantile( with(post, beta_SC_count), probs = PROBS))
SB<-as_tibble(quantile( with(post, beta_SB_count), probs = PROBS))
Breach_DO_SC<-as_tibble(quantile( with(post, beta_BreachDays*beta_DO*beta_SC_count), probs = PROBS)) ##  NS

# 2024-02-01
# ADD RAIN --> BREACH



#names for tibble
names<- c("RAIN -> Breach -> DO", "RAIN -> Breach^2", "BREACH -> DO", "WIND -> DO", "WIND -> Temp -> DO", "WIND -> DO -> SAV", "WIND -> DO -> SC",
          "DO -> SB", "YEAR", "SUBSTRATE -> SC", "SUBSTRATE", "RAIN -> DO -> SAV -> SB", "BREACH", 
          #"DO -> SB", "Year", "Substrate -> SC", "Rain -> DO -> SAV -> SB", "Breach", 
          "MICRO", "SAV -> SB", "SAV -> SC", "SC", "SB","BREACH -> DO -> SC", 
          "DO",
          "RAIN -> Breach -> Temp",
          #"Zone", "Zone -> Wind -> int", 
          "BREACH -> Temp",
          "WIND -> Temp",
          "WIND -> Temp -> SAV",
          "WIND -> Temp -> SC",
          "TEMP -> SB",
          "RAIN -> Temp -> SAV -> SB",
          "TEMP",
          "TEMP^2",
          "BREACH^2",
          "RAIN -> Breach",
          "YEAR^2",
          "SAV",
          "SAV^2")
#add probabilities
plot.posteriors<-rbind(Rain_Breach_DO, Rain_Breach_2, Breach_DO, Wind_DO, Wind_Temp_DO, Wind_DO_SAV, Wind_DO_SC,
        DO_SB, Year, Substrate_SC, Substrate, Rain_DO_SAV_SB, Breach, 
       # DO_SB, Year, Substrate_SC, Rain_DO_SAV_SB, Breach, 
        Micro, SAV_SB, SAV_SC, SC, SB, Breach_DO_SC, 
       DO,
       Rain_Breach_Temp,
        #Zone, Zone_Wind_int,
        Breach_Temp,
        Wind_Temp,
        Wind_Temp_SAV,
        Wind_Temp_SC,
        Temp_SB,
        Rain_Temp_SAV_SB,
        Temp,
        Temp_2,
        Breach_2,
       Rain_Breach,
       Year_2,
       SAV,
       SAV_2)
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

COLORS = c("red", "black", "blue")

plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
ggplot(aes(x = names, y = median, color = effect)) +
  #geom_point(effect = c("red", "black", "blue"))) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Causal Path") + 
  ylab("Causal Effect on Goby Density") +
  scale_color_manual(breaks = c("negative", "neutral", "positive"),
                     values=c("red", "darkgray", "green3")) + 
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme_gray(base_size = 16)


str(rethinking::extract.samples(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS))

Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS %>%
  spread_draws(a_Goby[Zone]) %>%
  head(10)




#effects plot
#k <- PSIS(Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct, pointwise = TRUE)$k
plot(dat$BreachDays, dat$Goby/dat$Area)
ns <- 100
P_seq <- seq( from = -1.2, to = 2.3, length.out = ns)
lambda <- link( Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS, data=data.frame(P=P_seq, cid=1))
lmu <- apply( lambda, 2, mean)
lci <- apply( lambda, 2, PI)
lines(P_seq, lmu, lty = 2, lwd = 1.5)
shade(lci, P_seq, xpd = TRUE)




ggplot(dat, aes(Year, y = (Goby/Area), color = as.factor(Zone))) + 
  geom_point() +
  geom_smooth() +
  ylab("density/log(m2)") +
  xlab("Year") +
  theme_classic(base_size=22)




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









