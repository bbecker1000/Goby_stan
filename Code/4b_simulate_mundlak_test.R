# Load necessary libraries
library(MASS)
library(dplyr)
library(broom)
library(posterior)
set.seed(123)


#get simulated_data from 4_simulate_mundlak.R
simulated_data <- read_csv("Output/Data/simulated_data.csv")



# Define beta coefficients
betas <- c(
  Rain = 0.5,
  Wind = -0.5,
  Substrate = -0.3,
  Temp = 0.4,
  Temp_2 = -0.4,
  DO = 0.6,
  SAV = 0.4,
  SAV_2 = -0.2,
  SC_count = -0.5,
  SB_count =  0.5,
  BreachDays = 0.3,
  BreachDays_2 = -0.1,
  Micro = 0.2,
  Year = 0.2,
  Year_2 = -0.2,
  Goby_lag = 0.1
)


library(glmmTMB)

model = glmmTMB(
  Goby ~ 
     Year +  
       Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
    #SAV_2 +
     Goby_lag +  
    DO +
     Temp +  
    Temp_2 +
    SC_count +
    Wind +
     BreachDays +  
    BreachDays_2 +
     Rain +
    (1 | Zone),
  offset = Area,
  family = nbinom2,
  data = simulated_data
)


options(scipen = 999)
summary(model)
options(scipen = 0)


# brms model
library(brms)
model.brms.sim <- brm(  # asks to install 2 items everytime run, say yes and wait ~120 seconds
  Goby ~ 
    Year +  
    Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
    #SAV_2 +
    Goby_lag +  
    DO +
    Temp +  
    Temp_2 +
    SC_count +
    Wind +
    BreachDays +  
    BreachDays_2 +
    Rain +
    offset(Area) + 
    (1 | Zone),
  family = negbinomial(link = "log", link_shape = "log"),
  control = list(adapt_delta = 0.95),
  iter = 10000, 
  warmup = 9000,
  data = simulated_data,
  #360 divergences when only 6000
  cores = 4
)

summary(model.brms.sim)

saveRDS(model.brms.sim, file = "Output/Models/model.brms.sim.rds")

#need to add some dummy data not used but named in the stan model (unused hold over variables that I never cleaned up)
simulated_data$BreachDays_Count <- 0
simulated_data$BreachDays_Count_2 <- 0
simulated_data$Breach_count <- 0

simulated_data$SC <- 0
simulated_data$SB <- 0
simulated_data$Year_int <- 0
simulated_data$Breach <- 0
#df_scaled$SAV_2 <- 0
head (simulated_data)


## now run the stan model on this data set to see if we can recover the coefficients.

#make data into a list
list_df_sim <- as.list(simulated_data)

# Need to Convert dataframe to list and flatten any 1-column matrices
list_df_sim <- lapply(as.list(simulated_data), function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(as.vector(x))
  } else {
    return(x)
  }
})

#model "mod.sim" must be compiled from 4a_stancode_simulate...
#updated beta priors to be 0 or closer to extreme coefficients like wind, SC, and SB

library(cmdstanr)

fit.sim.mundlak_MV_U <- mod.sim$sample(  # model from 2B_stan_code.R. or somewhere similar
  data = list_df_sim,
  seed = 123,
  chains = 1,
  iter_warmup = 5000, #need a lot to get past the divergences
  iter_sampling = 1000,
  parallel_chains = 1,
  refresh = 100, # print update every X iters
  adapt_delta = 0.93, #0.8 = lots divergences
  output_dir = "stan_output", # where to save the file
  save_warmup = FALSE # needed for "read_stan_csv+ to work, but can hangup if you save too much data
)

saveRDS(fit.sim.mundlak_MV_U, file = "Output/Models/fit.sim.mundlak_MV_U.rds")

rstan_fit.sim <- read_stan_csv(fit.sim.mundlak_MV_U$output_files())
fit.sim.mundlak_MV_U$diagnostic_summary()

# for saving and seeing
print(fit.sim.mundlak_MV_U$summary(), n = 30)

## get summary
coefficient_summary.sim <- summarise_draws(
  fit.sim.mundlak_MV_U,
  "mean",
  ~quantile(.x, probs = c(0.045, 0.955)) # For a 95% credible interval
)

#filter betas for joint model
beta.joint.sim <- coefficient_summary.sim %>%
  filter(str_starts(variable, "beta_")) %>%
  mutate(variable = str_remove(variable, "^.*?_"))   #remove up to the first underscore
beta.joint.sim$model <- "Shared (Mundlak)"
print(beta.joint.sim, n = 30) 

#get brm glmm model betas
coefficient_summary.sim.brm <- summarise_draws(
  model.brms.sim,
  "mean",
  ~quantile(.x, probs = c(0.045, 0.955)) # For a 95% credible interval
)
beta.sim.brm <- coefficient_summary.sim.brm %>%
  filter(str_starts(variable, "b_")) %>%
  filter(variable != "b_Intercept") %>%   # remove intercept
  mutate(variable = str_remove(variable, "^.*?_"))  #remove up to the first underscore
beta.sim.brm$model <- "Unshared"
print(beta.sim.brm, n = 30) 

#true from line 64
df.true <- rownames_to_column(as.data.frame(betas), var = "variable") 
df.true <- as_tibble(df.true)

df.true$model <- "True"

df.true$`4.5%` <- NA
df.true$`95.5%` <- NA
df.true <- df.true %>%
  rename(mean = betas)

df.true



df.sim <- rbind(df.true, beta.joint.sim, beta.sim.brm)
df.sim <- df.sim %>%
  rename(lower_4.5 = `4.5%`,
         upper_95.5 = `95.5%`)

#check variable names ok
unique(df.sim$variable)

p.sim.mundlak <- 
  df.sim %>%
  filter(variable != "SAV_2") %>%  #was not in the stan model...not in a submodel, so ok to leave in the unshared model.

ggplot(aes(x = variable, y = mean, color = model)) +
  geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), 
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  coord_flip() +
  #ylim(-2, 2) +
  #ylim(-1.5, 1.5) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic(base_size = 18) +
  theme(legend.position = c(0.80, 0.96),
        legend.title = element_blank(),
        legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = NA)) + 
  guides(color = guide_legend(reverse = TRUE)) + 
  ylab("beta estimate") +
  xlab(" ")

p.sim.mundlak

ggsave("Output/sim_mundlak_forest.png", width = 20, height = 30, units = "cm")


# now compare the stan fit with glmmTMB or create a BRMS model with no network.

