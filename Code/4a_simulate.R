# Load necessary libraries
library(MASS)
library(dplyr)
library(broom)
library(posterior)
set.seed(123)





# Number of observations
n <- 314

# Simulate Year as integers from 1 to 27 (representing 1997 to 2023)
year_levels <- 1:27
Year <- sample(year_levels, size = n, replace = TRUE)

# Simulate locations with random intercepts
locations <- c(1, 2, 3)
Location <- sample(locations, size = n, replace = TRUE)
location_intercepts <- rnorm(length(locations), mean = 0, sd = 0.5)
names(location_intercepts) <- locations

# Simulate independent variables with normal distributions
DO_raw <- rnorm(n, mean = 11, sd = 5)         # Approx. range 0–22
Wind_raw <- rnorm(n, mean = 5, sd = 2)
Temp_raw <- rnorm(n, mean = 19.5, sd = 2.5)   # Approx. range 14–25
SC_raw <- rbinom(n, 1, 0.5)
SB_raw <- rbinom(n, 1, 0.5)
Substrate_raw <- rbinom(n, 1, 0.5)
Micro_raw <- pmin(rpois(n, lambda = 5), 20)
SAV_raw <- rpois(n, lambda = 5)
PriorGoby_raw <- rpois(n, lambda = 10)
Area_raw <- rnorm(n, mean = 2.25, sd = 0.5)   # Approx. range 1–3.5
BreachDays_raw <- sample(0:62, size = n, replace = TRUE)
Rain_raw <- rnorm(n, mean = 26, sd = 10)      # Approx. range 2–50

# Create data frame
df <- data.frame(
  Year = Year,
  DO = DO_raw,
  Wind = Wind_raw,
  Temp = Temp_raw,
  SC = SC_raw,
  SB = SB_raw,
  Substrate = Substrate_raw,
  Micro = Micro_raw,
  SAV = SAV_raw,
  PriorGoby = PriorGoby_raw,
  Area = Area_raw,
  BreachDays = BreachDays_raw,
  Rain = Rain_raw,
  Location = Location
)

# Add quadratic terms
df$Year2 <- df$Year^2
df$Temp2 <- df$Temp^2
df$BreachDays2 <- df$BreachDays^2
df$logArea <- log(df$Area)

# Scale continuous and integer predictors
df_scaled <- df %>%
  mutate(across(c(Year, Year2, DO, Wind, Temp, Temp2, Micro, SAV, PriorGoby, BreachDays, BreachDays2, Rain), scale))

# Define beta coefficients
betas <- c(
  Year = 0.25,
  Year2 = -0.9,
  SB = 0.9,
  Substrate = 0.4,
  Micro = 0.3,
  SAV = 0.12,
  PriorGoby = 0,
  DO = 0.05,
  Temp = 0.2,
  Temp2 = -0.25,
  SC = -0.5,
  Wind = 0,
  BreachDays = 0.2,
  BreachDays2 = 0.0,
  Rain = 0.1
)

# Linear predictor with random intercepts and offset
eta <- with(df_scaled,
            betas["Year"] * Year +
              betas["Year2"] * Year2 +
              betas["SB"] * SB +
              betas["Substrate"] * Substrate +
              betas["Micro"] * Micro +
              betas["SAV"] * SAV +
              betas["PriorGoby"] * PriorGoby +
              betas["DO"] * DO +
              betas["Temp"] * Temp +
              betas["Temp2"] * Temp2 +
              betas["SC"] * SC +
              betas["Wind"] * Wind +
              betas["BreachDays"] * BreachDays +
              betas["BreachDays2"] * BreachDays2 +
              betas["Rain"] * Rain +
              location_intercepts[Location] +
              df$logArea
)

# Simulate Goby Count
mu <- exp(eta)
theta <- 1.5
df_scaled$Goby_Count <- rnegbin(n, mu = mu, theta = theta)


df_scaled <- df_scaled %>%
  rename(
    BreachDays = BreachDays,
    BreachDays_2 = BreachDays2,
    SAV = SAV,
    Year = Year,
    Goby = Goby_Count,
    Goby_lag = PriorGoby,
    Temp_2 = Temp2,
    Year_2 = Year2,
    Zone = Location,
    SB_count = SB,
    SC_count = SC
  )

head (df_scaled)

#need to add some dummy data not used but named in the stan model
df_scaled$BreachDays_Count <- 0
df_scaled$BreachDays_Count_2 <- 0
df_scaled$Breach_count <- 0

df_scaled$SC <- 0
df_scaled$SB <- 0
df_scaled$Year_int <- 0
df_scaled$Breach <- 0
df_scaled$SAV_2 <- 0
head (df_scaled)

library(glmmTMB)

model <- glmmTMB(
  Goby ~ 
     Year +  
       Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
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
  data = df_scaled
)

summary(model)


#try lme4
model.glmer <- glmer(
  Goby ~ 
    Year +  
    Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
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
  family = negative.binomial(1.73),  #theta from glmmTMB
  data = df_scaled
)

summary(model.glmer)

# brms model
library(brms)
model.brms <- brm(
  Goby ~ 
    Year +  
    Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
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
  data = df_scaled,
  iter = 6000,
  cores = 4
)

summary(model.brms)





## now run the stan model on this data set to see if we can recover the coefficients.

#make data into a list
list_df_sim <- as.list(df_scaled)

# Convert dataframe to list and flatten any 1-column matrices
list_df_sim <- lapply(as.list(df_scaled), function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(as.vector(x))
  } else {
    return(x)
  }
})


fit.sim <- mod$sample(  # model from 2B_stan_code.R.
  data = list_df_sim,
  seed = 123,
  chains = 3,
  iter_warmup = 2000,
  iter_sampling = 2000,
  parallel_chains = 3,
  refresh = 100, # print update every 500 iters
  output_dir = "stan_output", # where to save the file
  save_warmup = TRUE # needed for "read_stan_csv+ to work
)

rstan_fit.sim <- read_stan_csv(fit.sim$output_files())
fit.sim$diagnostic_summary()

# for saving and seeing
print(fit.sim$summary(), n = 30)

## get summary
coefficient_summary.sim <- summarise_draws(
  fit.sim,
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
  model.brms,
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
#fix up variable names and column names
df.true$variable <-
  ifelse(df.true$variable == "Year2", "Year_2",
         ifelse(df.true$variable == "SB", "SB_count",
                ifelse(df.true$variable == "PriorGoby", "Goby_lag",
                       ifelse(df.true$variable == "Temp2", "Temp_2",
                              ifelse(df.true$variable == "SC", "SC_count",
                                 ifelse(df.true$variable == "BreachDays2", "BreachDays_2",
                                     df.true$variable))))))
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

p.sim1 <- 
  df.sim %>%
  filter(variable != "SAV_2") %>%

ggplot(aes(x = variable, y = mean, color = model)) +
  geom_pointrange(aes(ymin = lower_4.5, ymax = upper_95.5), 
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  coord_flip() +
  #ylim(-2, 2) +
  #ylim(-1.5, 1.5) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic(base_size = 18) +
  theme(legend.position = c(0.2, 0.3)) + 
  ylab("beta estimate") 

p.sim1

ggsave("Output/sim1_forest.png", width = 20, height = 30, units = "cm")


# now compare the stan fit with glmmTMB or create a BRMS model with no network.
