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

theme_set(theme_tidybayes() + panel_border())

m <- Goby.m2.year.DAG.SC.SB_counts.BreachDays.direct

summary(m)

ggplot(dat, aes(Breach_count, y=Goby/Area, color = Zone)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(.~Zone)
  

m@stanfit

str(rethinking::extract.samples(m))

samples <- (rethinking::extract.samples(m))


m %>%
  spread_draws(a_Goby[Zone]) %>%
  head(10)

m %>%
  recover_types(dat) %>%
  spread_draws(a_Goby[Zone]) %>%
  head(10)

m %>%
  spread_draws(beta_BreachDays) %>%
  median_qi()


m %>%
  spread_draws(a_Goby[Zone]) %>%
  median_qi()

dat %>%
  data_grid(Goby, BreachDays) %>%
  add_linpred_draws(m) %>%
  head(10)





dat %>%
  group_by(Zone) %>%
  #data_grid(Year = seq_range(Year, n = 101)) %>%
  #add_linpred_draws(m, ndraws = 100) %>%
  ggplot(aes(x = BreachDays, y = Goby/Area, color = Zone)) +
  geom_smooth(method = "loess") +
  geom_point() +
  scale_color_brewer(palette = "Dark2")



dat %>%
  group_by(as.factor(Zone)) %>%
  data_grid(Year = seq_range(Year, n = 51)) %>%
  add_predicted_draws(m) %>%
  ggplot(aes(x = Year, y = Goby/Area, color = as.factor(Zone))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = dat) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")


#N.B. the syntax for compare_levels is experimental and may change
m %>%
  spread_draws(a_Goby[Zone]) %>%
  compare_levels(a_Goby, by = Zone) %>%
  ggplot(aes(y = Zone, x = a_Goby)) +
  stat_halfeye()


dat %>%
  group_by(Zone) %>%
  data_grid(Year = seq_range(Year, n = 101)) %>%
  add_linpred_draws(m) %>%
  ggplot(aes(x = Year, y = Goby/Area, color = Zone)) +
  stat_lineribbon(aes(y = .linpred)) +
  geom_point(data = dat) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")


set.seed(123456)
# to keep the example small we use 20 frames, 
# but something like 100 would be better
ndraws = 20


