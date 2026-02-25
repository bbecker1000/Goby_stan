# stan

library(cmdstanr)
library(bayesplot)

# Goby.m2.year.DAG.SC.SB_logistic.BreachDays.direct.RS.lag.2025_06_07
# string file
# fixed + to * for beta in Year_2 coefficient on 2025-09-13
stan_program <- "
data{
    array[314] int BreachDays_Count_2;
    array[314] int BreachDays_Count;
    array[314] int Breach_count;
     vector[314] Breach;
     vector[314] SC;
     vector[314] SB;
     vector[314] SAV_2;
    array[314] int Year_int;
    array[314] int Goby;
     vector[314] Area;
     vector[314] Goby_lag;
     vector[314] Temp_2;
     vector[314] BreachDays_2;
     vector[314] Micro;
     vector[314] Year_2;
     vector[314] Year;
    array[314] int Zone;
     vector[314] Wind;
    array[314] int SB_count;
    array[314] int SC_count;
    array[314] int Substrate;
     vector[314] BreachDays;
     vector[314] Rain;
     vector[314] SAV;
     vector[314] Temp;
     vector[314] DO;
}
parameters{
     real mu_Zone;
     real<lower=0> tau_Zone;
     real beta_Wind;
     real beta_Temp;
     real beta_Year_2;
     real beta_Year;
     real beta_SB_count;
     real a_Temp;
     real a_SAV;
     real a_SC;
     real a_SB;
     real a_DO;
     real a_BreachDays;
     vector[3] a_Goby;
     real beta_Temp_2;
     real beta_Micro;
     real beta_Rain;
     real beta_SC_count;
     real beta_SAV;
     real beta_SAV_2;
     real beta_DO;
     real beta_BreachDays;
     real beta_BreachDays_2;
     real beta_Substrate;
     real beta_Goby_lag;
     real<lower=0> tau;
     real phi;
}
model{
     vector[314] mu;
     vector[314] DO_nu;
     vector[314] Temp_nu;
     vector[314] SB_mu;
     vector[314] SC_mu;
     vector[314] Breach_nu;
     vector[314] SAV_nu;
    phi ~ normal( 1 , 5 );
    tau ~ exponential( 1 );
    beta_Goby_lag ~ normal( 0.25 , 0.25 );
    beta_Substrate ~ normal( 0.25 , 0.25 );
    beta_BreachDays_2 ~ normal( -0.1 , 0.25 );
    beta_BreachDays ~ normal( 0.25 , 0.25 );
    beta_DO ~ normal( 0.25 , 0.25 );
    beta_SAV_2 ~ normal( -0.1 , 0.25 );
    beta_SAV ~ normal( 0 , 0.25 );
    beta_SC_count ~ normal( -0.1 , 0.25 );
    beta_Rain ~ normal( 0.25 , 0.25 );
    beta_Micro ~ normal( 0.25 , 0.25 );
    beta_Temp_2 ~ normal( -0.1 , 0.25 );
    a_Goby ~ normal( 0 , 0.5 );
    a_BreachDays ~ normal( 0 , 0.5 );
    a_DO ~ normal( 0 , 0.5 );
    a_SB ~ normal( 0 , 0.5 );
    a_SC ~ normal( 0 , 0.5 );
    a_SAV ~ normal( 0 , 0.5 );
    a_Temp ~ normal( 0 , 0.5 );
    beta_SB_count ~ normal( 0 , 0.5 );
    beta_Year ~ normal( 0 , 0.5 );
    beta_Year_2 ~ normal( 0 , 0.5 );
    beta_Temp ~ normal( 0 , 0.5 );
    beta_Wind ~ normal( 0 , 0.5 );
    for ( i in 1:314 ) {
        SAV_nu[i] = a_SAV + beta_DO * DO[i] + beta_Temp * Temp[i];
    }
    SAV ~ normal( SAV_nu , tau );
    for ( i in 1:314 ) {
        Breach_nu[i] = a_BreachDays + beta_Rain * Rain[i];
    }
    BreachDays ~ normal( Breach_nu , tau );
    for ( i in 1:314 ) {
        SC_mu[i] = a_SC + beta_Substrate * Substrate[i] + beta_DO * DO[i] + beta_SAV * SAV[i];
        SC_mu[i] = inv_logit(SC_mu[i]);
    }
    SC_count ~ binomial( 1 , SC_mu );
    for ( i in 1:314 ) {
        SB_mu[i] = a_SB + beta_DO * DO[i] + beta_SAV * SAV[i];
        SB_mu[i] = inv_logit(SB_mu[i]);
    }
    SB_count ~ binomial( 1 , SB_mu );
    for ( i in 1:314 ) {
        Temp_nu[i] = a_Temp + beta_BreachDays * BreachDays[i] + beta_Wind * Wind[i];
    }
    Temp ~ normal( Temp_nu , tau );
    for ( i in 1:314 ) {
        DO_nu[i] = a_DO + beta_Temp * Temp[i] + beta_Wind * Wind[i];
    }
    DO ~ normal( DO_nu , tau );
    tau_Zone ~ exponential( 1 );
    mu_Zone ~ normal( 0 , 0.5 );
    a_Goby ~ normal( mu_Zone , tau_Zone );
    for ( i in 1:314 ) {
        mu[i] = a_Goby[Zone[i]] + beta_Year * Year[i] + beta_Year_2 * Year_2[i] + beta_SC_count * SC_count[i] + beta_SAV * SAV[i] + beta_SB_count * SB_count[i] + beta_DO * DO[i] + beta_Micro * Micro[i] + beta_BreachDays * BreachDays[i] + beta_BreachDays_2 * BreachDays_2[i] + beta_Substrate * Substrate[i] + beta_Wind * Wind[i] + beta_Temp * Temp[i] + beta_Temp_2 * Temp_2[i] + beta_Goby_lag * Goby_lag[i] + Area[i];
        mu[i] = exp(mu[i]);
    }
    Goby ~ neg_binomial_2( mu , exp(phi) );
}
"

# write to stan file
f <- write_stan_file(stan_program)

# compile model 
# about 3 minutes
mod <- cmdstan_model(f)
# look
mod$print()

data_list <- as.list(dat)


# Create a folder named "stan_output" in your current working directory
if (!dir.exists("stan_output")) {
  dir.create("stan_output")
}

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 3,
  iter_warmup = 2000,
  iter_sampling = 2000,
  parallel_chains = 3,
  refresh = 100, # print update every 500 iters
  output_dir = "stan_output", # where to save the file
  save_warmup = TRUE # needed for "read_stan_csv+ to work
)

rstan_fit <- read_stan_csv(fit$output_files())
fit$diagnostic_summary()

# for saving and seeing
print(fit$summary(), n = 30)
fit$metadata()$model_params
# print only
fit$cmdstan_summary()


# extract posterior draws for plotting and fun
# default is a 3-D draws_array object from the posterior package
# iterations x chains x variables
# use draws_df from stan for custom ggplots.
draws_df <- fit$draws(format = "df") # or format="array" (Only post warmup draws saved)
str(draws_df)


mcmc_hist(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_dens(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_dens_overlay(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_violin(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_intervals(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_areas(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))
mcmc_areas_ridges(fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")))

mcmc_combo(
  fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate")),
  combo = c("dens_overlay", "trace"),
  # pars = c("beta_BreachDays", "beta_Temp_2"),
  # transformations = list(sigma = "log"),
  gg_theme = legend_none()
)


mcmc_pairs(
  fit$draws(c("beta_BreachDays", "beta_Temp_2", "beta_Substrate"))
)

mcmc_scatter(
  fit$draws(c("beta_BreachDays", "beta_Temp_2"))
)

fit$diagnostic_summary()

# save object
fit$save_object(file = "Output/fit.RDS")

library(posterior)

summarise_draws(fit)
sjPlot::plot_model(fit)

# some plotting from bayesplot with help from gemini

plot_predictions(fit, newdata = dat, condition = "Substrate") +
  geom_point(data = dat, aes(x = x, y = y), alpha = 0.4) +
  theme_minimal() +
  labs(title = "Proof: `plot_predictions` on a `cmdstanr` Object")

rstantools::posterior_predict(fit)

library(tidyverse)


post <- draws_df
PROBS <- c(0.055, 0.5, 0.945) ## 89 % CIs
# Effect of Rain --> Breach --> DO --> Goby
Rain_Breach_DO <- as_tibble(quantile(with(post, beta_Rain * beta_BreachDays * beta_DO), probs = PROBS)) ## NS
# Effect of Breach --> DO --> Goby
Breach_DO <- as_tibble(quantile(with(post, beta_BreachDays * beta_DO), probs = PROBS)) ##  NS
# Effect of Wind --> Breach --> DO --> Goby
Wind_DO <- as_tibble(quantile(with(post, beta_Wind * beta_DO), probs = PROBS)) ## ns
# Effect of Wind --> Temp --> DO --> Goby
Wind_Temp_DO <- as_tibble(quantile(with(post, beta_Wind * beta_Temp * beta_DO), probs = PROBS)) ## ns
# Effect of Wind --> Temp --> SAV --> Goby
Wind_DO_SAV <- as_tibble(quantile(with(post, beta_Wind * beta_SAV * beta_DO), probs = PROBS)) ## ns
# Effect of Wind --> DO --> SB --> Goby
Wind_DO_SC <- as_tibble(quantile(with(post, beta_Wind * beta_DO * beta_SC_count), probs = PROBS)) ## ns
# Effect of DO --> SB --> Goby
DO_SB <- as_tibble(quantile(with(post, beta_DO * beta_SB_count), probs = PROBS)) ## ns

DO <- as_tibble(quantile(with(post, beta_DO), probs = PROBS)) ## ns

Rain_Breach_Temp <- as_tibble(quantile(with(post, beta_Rain * beta_BreachDays * beta_Temp), probs = PROBS)) ## NS
# Effect of Breach --> Temp --> Goby
Rain_Breach_2 <- as_tibble(quantile(with(post, beta_Rain * beta_BreachDays_2), probs = PROBS)) ## NS
# Effect of Rain --> Breach --> Goby

Rain_Breach <- as_tibble(quantile(with(post, beta_Rain * beta_BreachDays), probs = PROBS)) ## NS
# Effect of Rain --> Breach --> Goby

Breach_Temp <- as_tibble(quantile(with(post, beta_BreachDays * beta_Temp), probs = PROBS)) ##  NS
# Effect of Wind --> Breach --> Temp --> Goby
Wind_Temp <- as_tibble(quantile(with(post, beta_Wind * beta_Temp), probs = PROBS)) ## ns

# Effect of Wind --> Temp --> SAV --> Goby
Wind_Temp_SAV <- as_tibble(quantile(with(post, beta_Wind * beta_SAV * beta_Temp), probs = PROBS)) ## ns
# Effect of Wind --> Temp --> SB --> Goby
Wind_Temp_SC <- as_tibble(quantile(with(post, beta_Wind * beta_Temp * beta_SC_count), probs = PROBS)) ## ns
# Effect of Temp --> SB --> Goby
Temp_SB <- as_tibble(quantile(with(post, beta_Temp * beta_SB_count), probs = PROBS)) ## ns
Rain_Temp_SAV_SB <- as_tibble(quantile(with(post, beta_Rain * beta_Temp * beta_SAV * beta_SB_count), probs = PROBS))

Temp <- as_tibble(quantile(with(post, beta_Temp), probs = PROBS)) ## ns
Temp_2 <- as_tibble(quantile(with(post, beta_Temp_2), probs = PROBS)) ## ns

Breach_2 <- as_tibble(quantile(with(post, beta_BreachDays_2), probs = PROBS)) ## ns

Year <- as_tibble(quantile(with(post, beta_Year), probs = PROBS)) ## ns
Year_2 <- as_tibble(quantile(with(post, beta_Year_2), probs = PROBS)) ## ns


Substrate_SC <- as_tibble(quantile(with(post, beta_Substrate * beta_SC_count), probs = PROBS)) ## negative
Substrate <- as_tibble(quantile(with(post, beta_Substrate), probs = PROBS))
Rain_DO_SAV_SB <- as_tibble(quantile(with(post, beta_Rain * beta_DO * beta_SAV * beta_SB_count), probs = PROBS))

Breach <- as_tibble(quantile(with(post, beta_BreachDays), probs = PROBS)) # positive
Rain_Micro <- as_tibble(quantile(with(post, beta_Rain * beta_Micro), probs = PROBS))
Micro <- as_tibble(quantile(with(post, beta_Micro), probs = PROBS))
SAV_SB <- as_tibble(quantile(with(post, beta_SAV * beta_SB_count), probs = PROBS))
SAV_SC <- as_tibble(quantile(with(post, beta_SAV * beta_SC_count), probs = PROBS))
SAV <- as_tibble(quantile(with(post, beta_SAV), probs = PROBS))
SAV_2 <- as_tibble(quantile(with(post, beta_SAV_2), probs = PROBS))
SC <- as_tibble(quantile(with(post, beta_SC_count), probs = PROBS))
SB <- as_tibble(quantile(with(post, beta_SB_count), probs = PROBS))
Breach_DO_SC <- as_tibble(quantile(with(post, beta_BreachDays * beta_DO * beta_SC_count), probs = PROBS)) ##  NS
Goby_Lag <- as_tibble(quantile(with(post, beta_Goby_lag), probs = PROBS))
# 2024-02-01
# ADD RAIN --> BREACH



# names for tibble
names <- c(
  "RAIN → Breach → DO", "RAIN -> Breach^2", "BREACH -> DO", "WIND -> DO", "WIND -> Temp -> DO", "WIND -> DO -> SAV", "WIND -> DO -> SC",
  "DO -> SB",
  "YEAR",
  "SUBSTRATE -> SC", "SUBSTRATE", "RAIN -> DO -> SAV -> SB", "BREACH",
  # "DO -> SB", "Year", "Substrate -> SC", "Rain -> DO -> SAV -> SB", "Breach",
  "MICRO", "SAV -> SB", "SAV -> SC", "SC", "SB", "BREACH -> DO -> SC",
  "DO",
  "RAIN -> Breach -> Temp",
  # "Zone", "Zone -> Wind -> int",
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
  "SAV^2",
  "Goby_Lag"
)
# add probabilities
plot.posteriors <- rbind(
  Rain_Breach_DO, Rain_Breach_2, Breach_DO, Wind_DO, Wind_Temp_DO, Wind_DO_SAV, Wind_DO_SC,
  DO_SB,
  Year,
  Substrate_SC, Substrate, Rain_DO_SAV_SB, Breach,
  # DO_SB, Year, Substrate_SC, Rain_DO_SAV_SB, Breach,
  Micro, SAV_SB, SAV_SC, SC, SB, Breach_DO_SC,
  DO,
  Rain_Breach_Temp,
  # Zone, Zone_Wind_int,
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
  SAV_2,
  Goby_Lag
)
# add names
plot.posteriors$names <- rep(names, each = 3)
# add probabilities names
plot.posteriors$probability <- rep(c("lower", "median", "upper"), times = length(names))

plot.posteriors.wide <- plot.posteriors %>%
  group_by(probability) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = probability, values_from = value) %>%
  select(-row)

# add codes for positive or negative coefficients
plot.posteriors.wide$effect <- ifelse(
  plot.posteriors.wide$lower < 0 &
    plot.posteriors.wide$median < 0 &
    plot.posteriors.wide$upper < 0,
  "negative", ifelse(plot.posteriors.wide$lower > 0 &
    plot.posteriors.wide$median > 0 &
    plot.posteriors.wide$upper > 0,
  "positive",
  "neutral"
  )
)


print(plot.posteriors.wide, n = 31)

COLORS <- c("red", "black", "blue")

plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
  ggplot(aes(x = names, y = median, color = effect)) +
  # geom_point(effect = c("red", "black", "blue"))) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Causal Path") +
  ylab("Causal Effect on Goby Density") +
  scale_color_manual(
    breaks = c("negative", "neutral", "positive"),
    values = c("red", "darkgray", "green3")
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme_classic(base_size = 16)

ggsave("Output/forest.plot.Yearlag-stan.png", width = 20, height = 30, units = "cm")

# now plot with priors
priors <- read_csv("Data/Priors.csv")
# join
plot.posteriors.wide.priors <- left_join(plot.posteriors.wide, priors, by = "names")

## FIX THE LEGEND FOR THE PRIORS !
plot.posteriors.wide.priors %>% # plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
  ggplot(aes(x = names, y = median, color = effect)) +
  # geom_point(effect = c("red", "black", "blue"))) +

  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8,
    position = position_nudge(x = 0.1),
    show.legend = TRUE
  ) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(
    data = plot.posteriors.wide.priors,
    aes(
      x = names, y = prior,
      ymin = prior_lo, ymax = prior_hi
    ),
    color = "blue3", shape = "triangle", size = 0.5, alpha = 0.7,
    show.legend = FALSE,
    position = position_nudge(x = -0.1)
  ) +

  # geom_pointrange(aes(ymin = lower, ymax = upper)) +
  # geom_hline(yintercept = 0, lty = 2) +
  # geom_pointrange(data = ForestWithPriors,
  #                 aes(x = names, y = prior,
  #                     ymin = prior_lo, ymax = prior_hi), color = "lightblue") +

  xlab("Causal Path") +
  ylab("Causal Effect on Goby Density") +
  labs(color = "Posterior") +
  scale_color_manual(
    breaks = c("negative", "neutral", "positive"),
    values = c("red2", "darkgray", "green4")
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  annotate("pointrange", x = 15, y = 0.7, ymin = 0.5, ymax = 0.9,
            colour = "blue3", shape = "triangle", size = 0.5, linewidth = 0.5, alpha = 0.7) +
  annotate("text", x = 15, y = 1.1, label = "Prior", size = 6) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(.85, .5))
ggsave("Output/forest.plot.Yearlag.priors.png", width = 20, height = 30, units = "cm")

theme_few(base_size = 16)


# On to effects plots

# install.packages("marginaleffects")


library(marginaleffects)
library(posterior)

print(summarise_draws(fit), n = 30)


# Extract posterior draws
draws <- as_draws_df(fit$draws())
str(draws)

library(tidyverse)
x <- seq(from = -2, to = 2, by = 0.1)
y <- exp(-0.913 * x)
df <- tibble(x, y)

ggplot(df, aes(x, y)) +
  geom_line()
