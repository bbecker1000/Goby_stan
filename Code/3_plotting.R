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




par(mfrow=c(1,3))
for ( s in -2:2 ) {
  idx <- which( dat$Zone==1 )
  plot( dat$BreachDays[dat$Zone==1] , dat$Goby[dat$Zone==1], xlim = c(-2,2), ylim = c(0,2500) ,
    xlab="BreachDays", ylab="Goby", pch = 16, col = rangi2 )
    mu <- link( fit, data = data.frame(Zone=1 , BreachDays=-2:2 ) ) 
  for (i in 1:10) lines( -2:2, mu[i,], col = col.alpha("black, 0.3") )
}

names(as_tibble(link(fit)))


names(post)

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

d.Temp <- post %>% select(starts_with("Temp")) %>%
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

d.Breach <- post %>% select(starts_with("Breach")) %>%
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
d.all.post$Year <- rep(dat$Year, 4500)
d.all.post$Year_int <- rep(dat$Year_int, 4500)
d.all.post$Rain <- rep(dat$Rain, 4500)
d.all.post$Temp_2 <- rep(dat$Temp_2, 4500)
d.all.post$BreachDays_2 <- rep(dat$BreachDays_2, 4500)
d.all.post$Wind <- rep(dat$Wind, 4500)
d.all.post$Zone <- rep(dat$Zone, 4500)
d.all.post$Substrate <- rep(dat$Substrate, 4500)
d.all.post$Area <- rep(dat$Area, 4500) #already logged


#index for which data sample (cases = 314)
d.all.post$SAMPLE <- rep(1:4500, each = 314)

#Breach effects plot
d.all.post %>% filter(SAMPLE <= 50) %>%
ggplot(aes(x = Breach, mu, group = SAMPLE)) +
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,1500)

#Year effects plot
d.all.post %>% filter(SAMPLE <= 50) %>%
  ggplot(aes(x = Year_int+1994, mu/Area, group = SAMPLE, color = as_factor(Zone))) +
  geom_point(alpha = 0.2) + 
  stat_smooth (geom="line", color = "black", alpha=0.1, size=1) +
  #geom_smooth(aes(alpha = 0.1), method = "loess", se = FALSE) +
  ylim(0,1000)

#BreachDays_2 effects plot 
# check if squared the raw value before centering in model !!
d.all.post %>% filter(SAMPLE <= 50) %>%
  ggplot(aes(x = BreachDays_2, mu/Area, group = SAMPLE, color = as_factor(Zone))) +
  geom_point(alpha = 0.2) + 
  stat_smooth (geom="line", color = "gray", alpha=0.3, size=1) +
  #geom_smooth(aes(alpha = 0.1), method = "loess", se = FALSE) +
  ylim(0,1000)

#Wind effects plot 
d.all.post %>% filter(SAMPLE <= 50) %>%
  ggplot(aes(x = Wind, mu/Area, group = SAMPLE, color = as_factor(Zone))) +
  geom_point(alpha = 0.2) + 
  stat_smooth (geom="line", color = "black", alpha=0.1, size=1) +
  #geom_smooth(aes(alpha = 0.1), method = "loess", se = FALSE) +
  ylim(0,1000)

#Rain effects plot 
# check if squared the raw value before centering in model !!
d.all.post %>% filter(SAMPLE <= 50) %>%
  ggplot(aes(x = Rain, mu/Area, group = SAMPLE, color = as_factor(Zone))) +
  geom_point(alpha = 0.2) + 
  stat_smooth (geom="line", color = "black", alpha=0.1, size=1) +
  #geom_smooth(aes(alpha = 0.1), method = "loess", se = FALSE) +
  ylim(0,1000)

#substrate Effects plot
d.all.post %>% filter(SAMPLE <= 50) %>%
  ggplot(aes(x = as.factor(Substrate), mu/Area)) +
  #geom_boxplot() + 
  geom_jitter(alpha = 0.1, size = 0.1) 





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


#effects plot
mtcars_clean %>%
  group_by(cyl) %>%
  data_grid(hp = seq_range(hp, n = 364)) %>%
  add_predicted_draws(m_mpg) %>%
  ggplot(aes(x = hp, y = mpg, color = cyl, fill = cyl)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = mtcars_clean) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")


dat %>%
  data_grid(BreachDays = seq_range(BreachDays, n = 101)) %>%
  add_linpred_draws(fit) %>%
  ggplot(aes(x = BreachDays, y = .linpred)) +
  stat_lineribbon(color = "red") +
  scale_fill_brewer(palette = "Greys")


dat %>%
  group_by(Zone) %>%
  data_grid(BreachDays = seq_range(BreachDays, n = 364)) %>%
  add_predicted_draws(fit) %>%
  ggplot(aes(x = BreachDays, y = Goby, color = Zone, fill = Zone)) +
  #stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = dat) #+
  #scale_fill_brewer(palette = "Set2") +
  #scale_color_brewer(palette = "Dark2")



plot(fit)
plot(fit, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit, plotfun = "hist", pars = "theta", include = FALSE)
plot(fit, plotfun = "trace", pars = c("beta_BreachDays"), inc_warmup = TRUE)
plot(fit, plotfun = "rhat") + ggtitle("Example of adding title to plot")



m %>%
  spread_draws(intercept[condition]) %>%
  compare_levels(intercept, by = condition) %>%
  ggplot(aes(y = condition, x = intercept)) +
  stat_halfeye()

