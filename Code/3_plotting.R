library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(tidybayes.rethinking)
library(ggplot2)
library(cowplot)
library(rstan)
library(rethinking)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(brms)
library(sjPlot)
library(marginaleffects)

theme_set(theme_tidybayes() + panel_border())


#rename model for plotting
fit <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct.RS #m1 is stanfit

names(as_tibble(link(fit)))

post <- as.data.frame(link(fit))
#just get post-warmup values
post <- post[c(2501:5000),]

# Need to plot: 
# Goby(mu)/Area for:
#   Year
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

d.Breach <- post %>% 
  select(-contains("_2")) %>%     #remove _2
  select(starts_with("Breach")) %>%
  pivot_longer(
    cols = starts_with("Breach"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "Breach") %>%
  select(-'case')

d.SAV <- post %>% select(starts_with("SAV")) %>%
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
d.all.post$Year_int <- rep(dat$Year_int, 2500)
d.all.post$Rain <- rep(dat$Rain, 2500)
#d.all.post$Temp <- rep(dat$Temp, 2500)
d.all.post$Temp_2 <- rep(dat$Temp_2, 2500)
d.all.post$BreachDays_2 <- rep(dat$BreachDays_2, 2500)
d.all.post$Wind <- rep(dat$Wind, 2500)
d.all.post$Zone <- rep(dat$Zone, 2500)
d.all.post$Substrate <- rep(dat$Substrate, 2500)
d.all.post$Micro <- rep(dat$Micro, 2500)
d.all.post$Area <- rep(dat$Area, 2500) #already logged


# Zone plotting for facets
# 1 = E, 2 = NW, 3 = W
# New facet label names for dose variable
Zone.labs <- c("East", "Northwest", "West")
names(Zone.labs) <- c("1", "2", "3")

#index for which data sample (cases = 314)
d.all.post$SAMPLE <- rep(1:2500, each = 314)

#BREACH effects plot
#fix names
dat$Breach <- dat$BreachDays

#Subset data for the fit lines
d <- d.all.post %>% filter(SAMPLE <= 200) 

#%>%  #just take a few (50) samples from the data
  ggplot(data = d, aes(x = Breach, y = mu/exp(Area))) + #, group = SAMPLE
  #geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.03, size=0.5, color = "red") +
  geom_point(data = dat, aes(x = Breach, y = Goby/exp(Area)), alpha = 0.25, 
             color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("Breach Days (centered)") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#Year effects plot
  ggplot(data = d, aes(x = Year_int+1994, y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.05, color = "gray") + #posterior data

    stat_smooth (data = d, method = "lm", 
                 geom="line", aes(group = SAMPLE), alpha=0.03, size=0.5, color = "red") +
    geom_point(data = dat, aes(x = Year_int+1994, y = Goby/exp(Area)), alpha = 0.25, 
               color = "blue") + #raw data
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Gobys/m2") +
    xlab("Year") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  
  
  
#Year^2 effects plot
  ggplot(data = d, aes(x = Year_2, y = mu/exp(Area))) + #, group = SAMPLE
    geom_point(alpha = 0.05, color = "gray") + #posterior data
    stat_smooth (data = d, method = "lm", 
                 formula = y~poly(x,2),
                 geom="line", aes(group = SAMPLE), alpha=0.03, size=0.5, color = "red") +
    geom_point(data = dat, aes(x = Year_2, y = Goby/exp(Area)), alpha = 0.25, 
               color = "blue") + #raw data
    #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Gobys/m2") +
    xlab("Year") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs)) 
  

#BreachDays_2 effects plot 
# squared the raw value before centering in model 
ggplot(data = d, aes(x = BreachDays_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  stat_smooth (data = d, method = "lm", 
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE), alpha=0.03, size=0.5, color = "red") +
  geom_point(data = dat, aes(x = BreachDays_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("Breach Days^2") +
  facet_wrap(.~Zone)


#Wind effects plot 
# no effect on goby, only has a sig coefficient in model
ggplot(data = d, aes(x = Wind, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.5, color = "red") +
  geom_point(data = dat, aes(x = Wind, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("East/West Wind Mean") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#Rain effects plot 
ggplot(data = d, aes(x = Rain, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = Rain, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  stat_smooth (data = d, method = "lm", 
               #formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("Rainfall") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#substrate Effects plot
ggplot(data = d, aes(x = Substrate, y = mu/exp(Area))) + #, group = SAMPLE
    #geom_point(alpha = 0.05, color = "blue") + #posterior data
    #stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), alpha=0.05, size=0.5) +
    geom_jitter(data = dat, aes(x = Substrate, y = Goby/exp(Area), group = Zone), alpha = 0.25, color = "blue") + #raw data
    geom_boxplot(data = dat, aes(x = Substrate, y = Goby/exp(Area), group = Substrate), alpha = 0.25) + #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
    ylim(0,200) + 
    ylab("Gobys/m2") +
    xlab("Substrate") +
    facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))  


#SAV Effects plot
ggplot(data = d, aes(x = SAV, y = mu/exp(Area))) + #, group = SAMPLE
  #geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.2, size=1, 
               color = "red") +
  stat_smooth (data = d, method = "lm", 
               #formula = y~poly(x,2), 
               geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = SAV, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Gobys/m2") +
  xlab("SAV") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))


#SC Effects plot
ggplot(data = d, aes(x = SC_count/exp(Area), y = mu/exp(Area))) + #, group = SAMPLE
  #geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=1, 
               color = "blue") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = SC_count/exp(Area), y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,300) + 
  xlim(0,30) +
  ylab("Goby Density") +
  xlab("Sculpin Density") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))


#Micro Effects plot
ggplot(data = d, aes(x = Micro, y = mu/exp(Area))) + #, group = SAMPLE
  #geom_point(alpha = 0.05, color = "blue") + #posterior data
  stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
               alpha=0.05, size=1, 
               color = "blue") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  geom_point(data = dat, aes(x = Micro, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,300) + 
  ylab("Goby Density") +
  xlab("Microsporidia Count") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))


#Temp Effects plot
#needs fixing
ggplot(data = d, aes(x = Temp, y = mu/exp(Area))) + #, group = SAMPLE
  #geom_point(alpha = 0.05, color = "blue") + #posterior data
  # stat_smooth (data = d, method = "lm", geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=1, 
  #              color = "blue") +
  stat_smooth (data = d, method = "lm",
               #formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.5, color = "red") +
  geom_point(data = dat, aes(x = Temp, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Water Temperature") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))

#Temp_2 Effects plot

ggplot(data = d, aes(x = Temp_2, y = mu/exp(Area))) + #, group = SAMPLE
  geom_point(alpha = 0.05, color = "grey") + #posterior data
  stat_smooth (data = d, method = "lm",
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE),
               alpha=0.05, size=0.5, color = "red") +
  geom_point(data = dat, aes(x = Temp_2, y = Goby/exp(Area)), alpha = 0.25, color = "blue") + #raw data
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab("Goby Density") +
  xlab("Water Temperature") +
  facet_wrap(.~Zone, labeller = labeller(Zone = Zone.labs))



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



