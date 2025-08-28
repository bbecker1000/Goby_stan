
library(brms)


zinb <- read.csv("http://stats.idre.ucla.edu/stat/data/fish.csv")
zinb$camper <- factor(zinb$camper, labels = c("no", "yes"))
head(zinb)


fit_zinb2 <- brm(bf(count ~ persons + child + camper, zi ~ child), 
                 data = zinb, family = zero_inflated_poisson())



fit1 <- brm(
  count ~ zBase * Trt + (1|patient),
  data = epilepsy, family = poisson(),
  prior = prior(normal(0, 10), class = b) +
    prior(cauchy(0, 2), class = sd)
)



goby.brm <- brm(bf(
               Goby ~ Year + BreachDays + (1 | Zone)),
               BreachDays ~ Rain),
               data=dat, 
               family = negbinomial(link = "log", link_shape = "log"), 
               chains=3, 
               iter=3000
)

summary(goby.brm)


FORMULA <- bf(Goby ~ Year + BreachDays + (1 | Zone),
              BreachDays ~ Rain)
           
FAMILY <- c("negbinomial", "gaussian")

goby.brm.luxury <- brm(FORMULA,
                data=dat, 
                family = negbinomial,
                chains=3, 
                iter=3000
)

               


# Load the brms library
library(brms)

# Create some sample data
set.seed(123) # for reproducibility
dat <- data.frame(
  subject_id = factor(1:20),
  treatment = rep(c("A", "B"), each = 10),
  age = rnorm(20, mean = 30, sd = 5)
)

# Create two different outcomes
dat$reaction_time <- rnorm(20, mean = 500 - 10 * dat$age + 50 * (dat$treatment == "B"), sd = 50)
dat$errors_made <- rpois(20, lambda = exp(1 + 0.05 * dat$age - 0.5 * (dat$treatment == "B")))



# Define the formula for the continuous outcome (Gaussian)
bf_reaction <- bf(reaction_time ~ treatment + age, 
                  family = gaussian())

# Define the formula for the count outcome (Poisson)
bf_errors <- bf(errors_made ~ treatment + age, 
                family = poisson())

# The brm function will automatically combine them if you add them together
multi_formula <- bf_reaction + bf_errors

# Fit the multivariate model
multi_model <- brm(
  formula = multi_formula,
  data = dat,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 456,
  backend = "cmdstanr" # Recommended for speed
)


#try with goby

bf_Goby <- bf(Goby ~ Year + BreachDays + (1 | Zone), 
              family = negbinomial(link = "log", link_shape = "log"))

bf_BreachDays <- bf(BreachDays ~ Rain, family = gaussian)

multi_formula <- bf_Goby + bf_BreachDays


goby.brm <- brm(
  formula = multi_formula,
  data=dat, 
  chains=3, 
  iter=3000,
  backend = "cmdstanr"
)

summary(goby.brm)




