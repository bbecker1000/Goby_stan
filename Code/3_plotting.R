library(magrittr)
library(dplyr)
library(purrr)
#library(forcats)
library(tidyr)
#library(modelr)
#library(tidybayes)
#(ggplot2)
library(cowplot)
library(rethinking)
library(ggrepel)
library(RColorBrewer)
#library(gganimate)
library(brms)
#library(sjPlot)
#library(marginaleffects)
#library(rstan)
#library(tidybayes.rethinking)
#library(tidyverse)
library(ggthemes)

theme_set(theme_few())


#rename model for plotting
#fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS #m1 is stanfit

names(as_tibble(link(fit)))

post <- as.data.frame(link(fit))
#just get post-warmup values
post <- post[c(3001:6000),]   # if iter = 3000 post[c(2501:5000),]

# Need to plot: 
# Goby(mu)/Area for:
#   Year
#   Goby_lag
#   Temp_2
#   Temp
#   Substrate
#   SC
#   SAV_SC
#   Rain_Breach_Temp
#   Rain_Breach_2
#   Micro
#   Breach_Temp
#   Breach_2
#   Breach

# Goby predictions
d.mu <- post %>% select(starts_with("mu")) %>%
  pivot_longer(
    cols = starts_with("mu"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "mu")

d.DO <- post %>% select(starts_with("DO")) %>%
  pivot_longer(
    cols = starts_with("DO"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "DO") %>%
    select(-'case')

d.Temp <- post %>% 
  #select(-contains("_2")) %>%     #remove_2
  select(starts_with("Temp")) %>%
  pivot_longer(
    cols = starts_with("Temp"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "Temp") %>%
    select(-'case')

d.SB <- post %>% select(starts_with("SB")) %>%
  pivot_longer(
    cols = starts_with("SB"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "SB") %>%
    select(-'case')

d.SC <- post %>% select(starts_with("SC")) %>%
  pivot_longer(
    cols = starts_with("SC"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "SC") %>%
  select(-'case')

hist(d.SC$SC)

d.Breach <- post %>% 
  select(-contains("_2")) %>%     #remove _2
  select(starts_with("Breach")) %>%
  pivot_longer(
    cols = starts_with("Breach"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "Breach") %>%
  select(-'case')

d.SAV <- post %>% 
  select(-contains("_2")) %>%     #remove_2
  select(starts_with("SAV")) %>%
  pivot_longer(
    cols = starts_with("SAV"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "SAV") %>%
  select(-'case')


d.all.post <- bind_cols(d.mu, d.DO, d.Temp, d.SB, d.SC, d.Breach, d.SAV)
## 2024-03-12
## need to add raw data from dat file for unmodeled data
# Year
# Rain
# Temp_2
# BreachDays_2
# Wind
# Micro
# Zone
# Substrate
# Area

d.all.post$Year <- rep(dat$Year, nrow(d.all.post)/nrow(dat))
d.all.post$Year_2 <- rep(dat$Year_2, nrow(d.all.post)/nrow(dat))
d.all.post$Goby_lag <- rep(dat$Goby_lag, nrow(d.all.post)/nrow(dat))
d.all.post$Year_int <- rep(dat$Year_int, nrow(d.all.post)/nrow(dat))
d.all.post$Rain <- rep(dat$Rain, nrow(d.all.post)/nrow(dat))
#d.all.post$Temp <- rep(dat$Temp, nrow(dat))
d.all.post$Temp_2 <- rep(dat$Temp_2, nrow(d.all.post)/nrow(dat))
d.all.post$BreachDays_2 <- rep(dat$BreachDays_2, nrow(d.all.post)/nrow(dat))
d.all.post$Wind <- rep(dat$Wind, nrow(d.all.post)/nrow(dat))
d.all.post$Zone <- rep(dat$Zone, nrow(d.all.post)/nrow(dat))
d.all.post$Substrate <- rep(dat$Substrate, nrow(d.all.post)/nrow(dat))
d.all.post$Micro <- rep(dat$Micro, nrow(d.all.post)/nrow(dat))
d.all.post$SAV_2 <- rep(dat$SAV_2, nrow(d.all.post)/nrow(dat)) 
d.all.post$Area <- rep(dat$Area, nrow(d.all.post)/nrow(dat)) #already logged

## now add some original data for plotting on original scale
d.all.post$Year.unscaled <- 
  rep(dat.for.effects.plot$Year.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$SAV.unscaled <- 
  rep(dat.for.effects.plot$SAV.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$Temp.unscaled <- 
  rep(dat.for.effects.plot$Temp.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$DO.unscaled <- 
  rep(dat.for.effects.plot$DO.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$Rain.unscaled <- 
  rep(dat.for.effects.plot$Rain.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 

d.all.post$Breach.unscaled <- 
  rep(dat.for.effects.plot$Breach.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 


#and also to dat.plot
dat.plot <- dat

dat.plot$Year.unscaled <- 
  rep(dat.for.effects.plot$Year.unscaled)
dat.plot$SAV.unscaled <- 
  rep(dat.for.effects.plot$SAV.unscaled)
dat.plot$Temp.unscaled <- 
  rep(dat.for.effects.plot$Temp.unscaled) 
dat.plot$DO.unscaled <- 
  rep(dat.for.effects.plot$DO.unscaled) 
dat.plot$Rain.unscaled <- 
  rep(dat.for.effects.plot$Rain.unscaled) 
dat.plot$Breach.unscaled <- 
  rep(dat.for.effects.plot$Breach.unscaled) 






# Zone plotting for facets
# 1 = E, 2 = NW, 3 = W
# New facet label names for dose variable
Zone.labs <- c("East", "Northwest", "West")
names(Zone.labs) <- c("1", "2", "3")

#index for which data sample (cases = 314)
d.all.post$SAMPLE <- rep(1:3000, each = nrow(dat))


#BREACH effects plot
#fix names
dat$Breach <- dat$BreachDays

#attach raw data to df for the SC logit plot.
d.all.post$SC_count <- rep(dat$SC_count, times = 3000)
d.all.post$Breachdays <- rep(dat$BreachDays, times = 3000)
d.all.post$SB_count <- rep(dat$SB_count, times = 3000)
d.all.post$DO <- rep(dat$DO, times = 3000)
d.all.post$SAV <- rep(dat$SAV, times = 3000)
d.all.post$Temp_2 <- rep(dat$Temp_2, times = 3000)







#Subset data for the fit lines
# grab random 100 samples near middle of the chain
d <- d.all.post %>% filter(between(SAMPLE, 1950 , 2001) )
d$SB_count


#Breach using the raw breach data as x variable.
p.breach <-  ggplot(data = d, aes(x = Breach.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "grey") + #posterior data
    geom_point(data = dat.plot, aes(x = Breach.unscaled, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data

  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, linewidth=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Annual Breach Days") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.breach

#Year effects plot
p.year <-  ggplot(data = d, aes(x = Year_int+1994, y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.1, color = "gray") + #posterior data

    stat_smooth (data = d, method = "lm", 
                 geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
    geom_point(data = dat.plot, aes(x = Year_int+1994, y = Goby/exp(Area)), alpha = 0.25, 
               color = "blue") + #raw data
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab(" ") +
    xlab("Year") +
    scale_x_continuous(breaks = c(2000, 2010, 2020 )) +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.year
  



#Year^2 effects plot
p.year2 <-  ggplot(data = d, aes(x = Year_int+1995, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "gray") + #posterior data
    geom_point(data = dat.plot, aes(x = Year_int+1995, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data

    stat_smooth (data = d, method = "lm", 
                 formula = y~poly(x,2),
                 geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab(" ") +
    xlab("Year") +
    guides(x = guide_axis(minor.ticks = TRUE)) +
    scale_x_continuous(breaks = c(2000, 2020)) +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs)) 
p.year2

#Goby_lag effects plot
p.Goby_lag <-  ggplot(data = d, aes(x = Goby_lag, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "gray") + #posterior data
    geom_point(data = dat.plot, aes(x = Goby_lag, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", 
               geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  #xlab("Scaled Goby Density(t-1)") +
  #ylab(expression(paste("Scaled Goby Density, t-1, ")                   "))) 
  xlab(bquote("Scaled Goby Density"[t-1])) +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.Goby_lag

#BreachDays_2 effects plot 
#plotting vs breach with a poly to show the breach^2 model
# squared the raw value before centering in model 
p.breach2 <- ggplot(data = d, aes(x = BreachDays_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(data = dat.plot, aes(x = BreachDays_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  geom_point(alpha = 0.1, color = "gray") + #posterior data
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE), alpha=0.05, size=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Breach Days") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.breach2




#DO effects plot

p.DO <-  ggplot(data = d, aes(x = DO.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "gray") + #posterior data
  geom_point(data = dat.plot, aes(x = DO.unscaled, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", 
               geom="line", aes(group = SAMPLE), alpha=0.05, 
               linewidth=0.75, color = "red") +
  ylim(0,200) + 
  ylab(" ") +
  xlab("DO") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.DO




#Wind effects plot 
# no effect on goby, only has a sig coefficient in model
p.wind <- ggplot(data = d, aes(x = Wind, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(data = dat.plot, aes(x = Wind, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("East/West Wind Mean") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.wind


#Rain effects plot 
p.rain <- ggplot(data = d, aes(x = Rain.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "gray") + #posterior data
  #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
  geom_point(data = dat.plot, aes(x = Rain.unscaled, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", 
               #formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), alpha=0.05, size=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Annual Rainfall (cm)") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.rain


#substrate Effects plot
p.substrate <- ggplot(data = d, aes(x = as.factor(Substrate), y = mu/exp(Area))) + #, group = SAMPLE
    geom_jitter(alpha = 0.2, color = "grey", width=0.1) + #posterior data
    #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
    geom_jitter(data = dat.plot, aes(x = as.factor(Substrate), y = Goby/exp(Area), group = Zone), 
                alpha = 0.25, color = "blue", width = 0.1) + #raw data
    stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
    #geom_boxplot(data = dat.plot, aes(x = as.factor(Substrate), y = Goby/exp(Area), group = Substrate), alpha = 0.25) + #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
    #           alpha=0.05, size=0.75, 
    #           color = "red") +stat_summary(geom = "point", fun.y = "mean", col = "red", size = 2, shape = "square", alpha = 0.8) + 
  ylim(0,200) + 
    ylab(" ") +
    xlab("Substrate") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.substrate

#SAV Effects plot
p.sav <- ggplot(data = d, aes(x = SAV.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_jitter(alpha = 0.1, color = "grey", width = 0.05) + #posterior data
  geom_jitter(data = dat.plot, aes(x = SAV.unscaled, y = Goby/exp(Area)), alpha = 0.25, 
              color = "blue", width = 0.05) +
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=1, 
               color = "red") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.75, color = "red") +
 
   #           color = "blue", width = 0.1) + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("SAV Cover Quintile") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.sav

#SAV Effects plot
p.sav2 <- ggplot(data = d, aes(x = SAV_2, y = mu/exp(Area))) + 
  geom_point(data = dat.plot, aes(x = SAV_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "grey") + #posterior data
  # stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
  #              alpha=0.2, size=1, 
  #              color = "red") +
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, color = "red") +#raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("SAV^2") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.sav2



#SC Effects plot
#make sure run binary code before logistic model !
p.SC <- ggplot(data = d, aes(x = as.factor(SC_count), y = mu/exp(Area))) + #, group = SAMPLE
  geom_jitter(alpha = 0.1, color = "gray", width = 0.1) + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
  geom_jitter(data = dat.plot, aes(x = SC_count+1, y = Goby/exp(Area)), 
             alpha = 0.25, color = "blue", width = 0.1) + #raw data
  #geom_boxplot() + 
  #stat_summary(geom = "point", fun.y = "mean", col = "red", size = 2, shape = "square", alpha = 0.8) + 
  ylim(0,200) + 
  ylab(" ") +
  xlab("Sculpin Presence") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.SC

#make sure run binary code before logistic model !
# file 2_Full Luxury, line ~177
p.SB <- ggplot(data = d, aes(x = as.factor(SB_count), y = mu/exp(Area))) + #, group = SAMPLE
  geom_jitter(alpha = 0.1, color = "gray", width = 0.1) + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
  geom_jitter(data = dat.plot, aes(x = SB_count+1, y = Goby/exp(Area)), 
              alpha = 0.25, color = "blue", width = 0.1) + #raw data
  #stat_summary(geom = "point", fun.y = "mean", col = "red", size = 2, shape = "square", alpha = 0.8) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Stickleback Presence") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.SB




#Micro Effects plot
p.micro <- ggplot(data = d, aes(x = Micro, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "gray") + #posterior data
  geom_point(data = dat.plot, aes(x = Micro, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
 # ylab("Goby Density (n/m^2") +
  ylab(expression(paste("Goby Density (", m^2, ")                   "))) + 
  xlab("Microsporidia Count (scaled)") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.micro

#Temp Effects plot
#needs fixing
ggplot(data = d, aes(x = Temp.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "grey") + #posterior data
  # stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=1, 
  #              color = "blue") +
  geom_point(data = dat.plot, aes(x = Temp.unscaled, y = Goby/exp(Area)), alpha = 0.25, color = "blue") +
  stat_smooth (data = d, method = "lm",
               #formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.75, color = "red") +
  #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  #xlab("Water Temperature (C)") +
  xlab(paste0("Water Temperature" , "\u00B0", "C")) + 
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#Temp_2 Effects plot

p.temp2 <- ggplot(data = d, aes(x = Temp.unscaled, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.1, color = "grey") + #posterior data
  geom_point(data = dat.plot, aes(x = Temp.unscaled, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + 
  stat_smooth (data = d, method = "lm",
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.75, color = "red") +
  #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  #xlab("Water Temperature") +
  xlab(paste0("Water Temperature" , " (\u00B0", "C)")) + 
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.temp2 

## panel plot
## want: Breach, Year_2, SB, Micro, Substrate, SAV, Goby_lag, Temp_2, DO) 

p.all.effects <- cowplot::plot_grid(p.year2,
                                    p.SC, 
                                    p.SB, 
                                    
                                    p.micro, 
                                    p.substrate,
                                    p.sav,
                                  
                                    p.Goby_lag,
                                    p.temp2, 
                                    p.breach,
                                    
                                    p.DO, 
                                    p.rain, 
                           ncol=3, labels="AUTO", scale = 0.9, 
                           vjust = 3, hjust = -2.2
                           )
p.all.effects
ggsave("Output/p.all.effects.lag.jpg", width = 35, height = 30, units = "cm")



#random effects groups
fit %>%
  spread_draws(a_Goby[Zone]) %>%
  median_qi()

# plot RE.
fit %>%
  spread_draws(a_Goby[Zone]) %>%
  median_qi() %>%
  ggplot(aes(y = Zone, x = a_Goby, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()





