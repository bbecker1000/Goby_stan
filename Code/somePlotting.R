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

library(rethinking)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(rstan)

theme_set(theme_tidybayes() + panel_border())

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m <- Goby.m2.year.DAG.SC.SB_counts.Breach.pois
str(m)
summary(m)
m <- rethinking::extract.samples(m) # (these are the posterior betas)
m

sim.m <- rethinking::sim(m, n = 100)
sim.m


# convert to a long-tibble 

colnames(sim.m) <- paste0("C",1:ncol(sim.m))
sim.m[1:5,1:5]

as_tibble(sim.m)

as_tibble(sim.m) %>%
  mutate(row_id=factor(row_number()))

mat_df <- as_tibble(sim.m) %>%
  mutate(row_id=factor(row_number())) %>%
  relocate(row_id)

mat_df <- mat_df %>% 
  pivot_longer(-row_id,
               names_to="sample_id",
               values_to="vals") 

#add rownumbers to dat
dat_df <- as_tibble(dat) %>%
  mutate(row_id=factor(row_number())) %>%
  relocate(row_id)

# now join with our raw data
dat_plus_sim <- left_join(mat_df, dat_df, by = "row_id")

ggplot(dat_plus_sim, aes(Year_int+1995, vals, color = Zone)) +
  geom_smooth(method = "gam")


ggplot(dat_plus_sim, aes(SAV, vals, color = Zone)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(SC_count/Area, vals, color = Zone)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(SB_count/Area, vals, color = Zone)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(Breach_count, vals)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(Rain, vals)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(Temp, vals)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(DO, vals)) +
  geom_smooth(method = "lm")

ggplot(dat_plus_sim, aes(Wind, vals)) +
  geom_smooth(method = "lm")
  
    