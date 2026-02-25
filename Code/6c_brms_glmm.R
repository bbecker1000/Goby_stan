#compare the simulations

library(brms)
model.brms.sim <- brm(  # asks to install 2 items everytime run, say yes and wait ~90 seconds
  Goby ~ 
    Year +  
    Year_2 +
    SB_count + 
    Substrate +  
    Micro +  
    SAV +
    SAV_2 +
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
  control = list(adapt_delta = 0.9),
  iter = 4000, 
  warmup = 3000,
  data = sim_data,
  chains = 1,
  #360 divergences when only 6000
  cores = 4
)

summary(model.brms.sim)

saveRDS(model.brms.sim, file = "Output/Models/model.brms.sim.rds")
