library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rethinking)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(brms)
library(sjPlot)
library(marginaleffects)
library(cowplot)
library(rstan)
library(tidybayes.rethinking)

theme_set(theme_tidybayes() + panel_border())


#rename model for plotting
#fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS #m1 is stanfit

names(as_tibble(link(fit)))

post <- as.data.frame(link(fit))
#just get post-warmup values
post <- post[c(3001:6000),]   # if itel = 3000 post[c(2501:5000),]

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


# Zone plotting for facets
# 1 = E, 2 = NW, 3 = W
# New facet label names for dose variable
Zone.labs <- c("East", "Northwest", "West")
names(Zone.labs) <- c("1", "2", "3")

#index for which data sample (cases = 301)
d.all.post$SAMPLE <- rep(1:3000, each = nrow(dat))


#BREACH effects plot
#fix names
dat$Breach <- dat$BreachDays

#Subset data for the fit lines
# grab random 100 samples near middle of the chain
d <- d.all.post %>% filter(between(SAMPLE, 1975 , 2001) )


#%>%  #just take a few (50) samples from the data
p.breach <-  ggplot(data = d, aes(x = Breach, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "grey") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, linewidth=0.75, color = "red") +
  geom_point(data = dat, aes(x = Breach, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Breach Days") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.breach

#Year effects plot
p.year <-  ggplot(data = d, aes(x = Year_int+1994, y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.05, color = "gray") + #posterior data

    stat_smooth (data = d, method = "lm", 
                 geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
    geom_point(data = dat, aes(x = Year_int+1994, y = Goby/exp(Area)), alpha = 0.25, 
               color = "blue") + #raw data
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Gobys/m2") +
    xlab("Year") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.year
  
  
#Year^2 effects plot
p.year2 <-  ggplot(data = d, aes(x = Year_2, y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.05, color = "gray") + #posterior data
    stat_smooth (data = d, method = "lm", 
                 formula = y~poly(x,2),
                 geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
    geom_point(data = dat, aes(x = Year_2, y = Goby/exp(Area)), alpha = 0.25, 
               color = "blue") + #raw data
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Goby Density (m2)") +
    xlab("Year") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs)) 
p.year2

#Goby_lag effects plot
p.Goby_lag <-  ggplot(data = d, aes(x = Goby_lag, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  
  stat_smooth (data = d, method = "lm", 
               geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  geom_point(data = dat, aes(x = Goby_lag, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("Goby Density Year-1") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.Goby_lag

#BreachDays_2 effects plot 
#plotting vs breach with a poly to show the breach^2 model
# squared the raw value before centering in model 
p.breach2 <- ggplot(data = d, aes(x = BreachDays_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE), alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = BreachDays_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Breach Days") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.breach2

#Wind effects plot 
# no effect on goby, only has a sig coefficient in model
p.wind <- ggplot(data = d, aes(x = Wind, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = Wind, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("East/West Wind Mean") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.wind


#Rain effects plot 
p.rain <- ggplot(data = d, aes(x = Rain, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = Rain, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", 
               #formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), alpha=0.05, size=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Rainfall") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.rain


#substrate Effects plot
p.substrate <- ggplot(data = d, aes(x = as.factor(Substrate), y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.2, color = "grey") + #posterior data
    #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
    geom_jitter(data = dat, aes(x = as.factor(Substrate), y = Goby/exp(Area), group = Zone), 
                alpha = 0.25, color = "blue", width = 0.2) + #raw data
    #geom_boxplot(data = dat, aes(x = as.factor(Substrate), y = Goby/exp(Area), group = Substrate), alpha = 0.25) + #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Goby Density") +
    xlab("Substrate") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
p.substrate

#SAV Effects plot
p.sav <- ggplot(data = d, aes(x = SAV, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.2, size=1, 
               color = "red") +
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = SAV, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("SAV") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.sav

#SAV Effects plot
p.sav2 <- ggplot(data = d, aes(x = SAV_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "grey") + #posterior data
  # stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
  #              alpha=0.2, size=1, 
  #              color = "red") +
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = SAV_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("SAV^2") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.sav2



#SC Effects plot
#needs debugging
p.SC <- ggplot(data = d, aes(x = SC, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = SC, y = Goby/exp(Area)), 
             alpha = 0.25, color = "blue") + #raw data
  geom_smooth(data = dat, aes(x = SC, y = Goby/exp(Area)),
              method = "lm", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Sculpin Density") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.SC



#Micro Effects plot
p.micro <- ggplot(data = d, aes(x = Micro, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.75, 
               color = "red") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = Micro, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Microsporidia Count") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.micro

#Temp Effects plot
#needs fixing
ggplot(data = d, aes(x = Temp, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "grey") + #posterior data
  # stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=1, 
  #              color = "blue") +
  stat_smooth (data = d, method = "lm",
               #formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = Temp, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Water Temperature") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#Temp_2 Effects plot

p.temp2 <- ggplot(data = d, aes(x = Temp_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "grey") + #posterior data
  stat_smooth (data = d, method = "lm",
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.75, color = "red") +
  geom_point(data = dat, aes(x = Temp_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Water Temperature") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))
p.temp2

## panel plot
p.all.effects <- plot_grid(p.breach, p.temp2, p.sav2, p.rain, p.micro, p.Goby_lag, p.year2)
p.all.effects
ggsave("Output/p.all.effects.lag.png", width = 30, height = 20, units = "cm")


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

# plot RE with stats, takes a long time.
fit %>%
  spread_draws(a_Goby[Zone]) %>%
  ggplot(aes(y = Zone, x = a_Goby)) +
  stat_eye()
beepr::beep()


#effects plot example
mtcars_clean %>%
  group_by(cyl) %>%
  data_grid(hp = seq_range(hp, n = 364)) %>%
  add_predicted_draws(m_mpg) %>%
  ggplot(aes(x = hp, y = mpg, color = cyl, fill = cyl)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = mtcars_clean) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")

#diagnostic plots
plot(fit)
plot(fit, plotfun = "trace", pars = c("beta_BreachDays"), inc_warmup = TRUE)
plot(fit, plotfun = "rhat") + ggtitle("Example of adding title to plot")



