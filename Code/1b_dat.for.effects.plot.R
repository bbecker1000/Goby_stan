## Data Prep for plotting effects on regular scale

##remove cases with NAs
dat.for.effects.plot <- as_tibble(dat.for.effects.plot)
nrow(dat.for.effects.plot) #364 rows

## add original scale covariates
# Year
# SAV
# Goby_Lag
# Temp
# breach_days
# DO
# Rain


dat.for.effects.plot$Year.unscaled <- rep(dat.temp$Year)
dat.for.effects.plot$SAV.unscaled <- rep(dat.temp$SAV)
dat.for.effects.plot$Temp.unscaled <- rep(dat.temp$temp_mean)
dat.for.effects.plot$DO.unscaled <- rep(dat.temp$min_DO)
dat.for.effects.plot$Year.unscaled <- rep(dat.temp$Rain_Sum)



















dat.for.effects.plot <- na.omit(dat.for.effects.plot) #, most missingness due to DO
nrow(dat.for.effects.plot) #314

# now can join other covariates from 
