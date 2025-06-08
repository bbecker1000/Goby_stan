#data_prep

## end examples ------------

## note code offset using simply "+ log(Area)"

#Goby

library(readxl)
library(readr)
library(tidyverse)
library(lme4)
library(rethinking)
library(rstan)
#library(rstanarm)
library(sjPlot)
library(marginaleffects)
goby_master <- read_csv("Data/goby_master_2023.csv")

hist(goby_master$Area)

breach_days <- read_excel("Data/RodeoLagoon-Status_WY1995_2023.xlsx", 
                          col_types = c("date", "numeric", "numeric", "text", "text"))
#keep first three columns
breach_days <- breach_days[,c(1:3)]
#only keep value == 1
breach_days_sum <- breach_days %>%
  filter(RodeoLagoonMouth == 1) %>%
  group_by(WaterYear) %>%
  summarize(BreachDays = sum(RodeoLagoonMouth))

#find missing years
#2010, 2014 2020 2021  = 0 breaches
ZeroYears <- tibble(WaterYear = c(2010, 2014, 2020, 2021), BreachDays = c(0,0,0,0))

breach_days_sum2 <- rbind(breach_days_sum, ZeroYears)
breach_days_sum3 <- breach_days_sum2 %>% arrange(WaterYear)
breach_days_sum <- breach_days_sum3

#key for join
breach_days_sum$Year <- breach_days_sum$WaterYear

#Join Breach Days to goby master 
# WY is previous 12 months ...WY1995 = 1994-10-01 to 1995-09-30
# so match WY to goby sampling year...
goby_masterB <- goby_master %>% left_join(breach_days_sum, by = "Year")
goby_master <- goby_masterB

#goby_master$Since_Breach[is.na(goby_master$Since_Breach)] <- 0
goby_master$Since_Breach <- ifelse(goby_master$Since_Breach == 0.5, 1, 
                                   goby_master$Since_Breach)

#make since breach a 0/1
#wonky categorical compared to rstanarm, so need to fix later
#goby_master$Since_Breach <- ifelse(goby_master$Since_Breach ==  0, 0, 1)

goby_master$Since_Breach

hist(goby_master$micro_sum)
summary(goby_master$micro_sum)

dat.temp <- goby_master


hist(dat.temp$Year)

#2024-03-19 per DF
#remove 2021 due to poor sampling conditions (flood)

#dat.temp <- dat.temp %>% filter(Year != 2021)
hist(dat.temp$Year)

## Prep Data
#IVs
Goby   <- dat.temp$Sum_TW  # must be non-negative
Year   <- scale(dat.temp$Year-1995)
Year_2   <- scale((dat.temp$Year-1995)^2)
Year_int <- as.integer(dat.temp$Year-1995)
SAV    <- scale(dat.temp$SAV)
SAV_2    <- scale(dat.temp$SAV^2)
SB     <- scale(dat.temp$Sum_SB)
SB_count     <- dat.temp$Sum_SB # if using counts
SC     <- scale(dat.temp$Sum_SC)
SC_count     <- dat.temp$Sum_SC
Rain   <- scale(dat.temp$Rain_Sum)
Temp   <- scale(dat.temp$temp_mean)
Temp_2   <- scale(dat.temp$temp_mean^2)
DO     <- scale(dat.temp$min_DO)
Micro  <- scale(dat.temp$micro_sum) #added 2024-01-23, poor fit if raw counts
Breach <- scale(dat.temp$Since_Breach)
Breach_count <- dat.temp$Since_Breach
BreachDays <- scale(dat.temp$BreachDays)
BreachDays_Count <- (dat.temp$BreachDays)
BreachDays_2 <- scale(dat.temp$BreachDays^2)
BreachDays_Count_2 <- (dat.temp$BreachDays^2)
#Breach <- as.factor(dat.temp$Since_Breach)  # need to fix this, but categorical is wonky
Wind   <- scale(dat.temp$u_mean)

Zone   <- (dat.temp$Zone)
Substrate   <- (dat.temp$Dom_substrate)
#Random Effects
unique(dat.temp$Zone)
# zone to integer, keep as three zones
# stan requires groups not be named "0"
Zone   <- as.integer(ifelse(dat.temp$Zone == "E", 1,
                            ifelse(dat.temp$Zone == "NW", 2,
                                   ifelse(dat.temp$Zone == "W", 3, dat.temp$Zone))))
unique(dat.temp$Dom_substrate)

plot(dat.temp$Sum_TW, Zone)

#pool substrates
# stan requires groups not be named "0"
Substrate <- as.integer(ifelse(dat.temp$Dom_substrate == "corophium_tubes", "2",
                               ifelse(dat.temp$Dom_substrate == "mud", "1",
                                      ifelse(dat.temp$Dom_substrate == "muck", "1",
                                             ifelse(dat.temp$Dom_substrate == "gravel", "2",
                                                    ifelse(dat.temp$Dom_substrate == "sand", "2",
                                                           ifelse(dat.temp$Dom_substrate == "cobble", "2",
                                                                  ifelse(dat.temp$Dom_substrate == "riprap", "2",
                                                                         dat.temp$Dom_substrate))))))))

#Offset
Area <- log(dat.temp$Area)


dat <- data.frame(Goby=Goby, Year=Year, Year_2=Year_2, Year_int = Year_int, SAV=SAV, SAV_2=SAV_2,
                  SB=SB, SB_count=SB_count, SC=SC, 
                  SC_count=SC_count, Rain=Rain, Temp=Temp, Temp_2=Temp_2, #added 2024-02-09
                  DO=DO, Breach=Breach,  Breach_count=Breach_count, #added 2024-01-23
                  BreachDays=BreachDays, BreachDays_2=BreachDays_2,
                  BreachDays_Count=BreachDays_Count, BreachDays_Count_2=BreachDays_Count_2,
                  Wind=Wind, Micro=Micro,
                  Zone=Zone, Substrate=Substrate, Area=Area)



## 2024-03-24
## create lagged t-1 overall TWG/AREA density to account for 
## mechanistic larval supply
## in lieu of ar()

dat$GobyDens <- dat$Goby/dat$Area
hist(dat$GobyDens) 
hist(dat$Year_int+1994, n = 27) 

lag_mean_dens <- dat|>
  group_by(Year_int) |>
    summarize(Goby_lag = mean(GobyDens, na.rm=TRUE))
 
lag_mean_dens %>%
 ggplot(aes(Year_int, Goby_lag)) +
  geom_point()

#link to dat
lag_mean_dens$Year_int <- lag_mean_dens$Year_int+1

dat <- left_join(dat, lag_mean_dens, by = "Year_int")

# give year 1 the mean lag density for all years
year_1_lag <- mean(dat$Goby_lag, na.rm=TRUE)
# put in data file for year 1  
dat$Goby_lag <- ifelse(dat$Year_int == 1 | dat$Year_int == 10, year_1_lag, dat$Goby_lag)
#scale
dat$Goby_lag <- scale(dat$Goby_lag)

#remove the survey level density
dat <- dat |> select(-GobyDens)

dat$Goby_lag <- as.numeric(dat$Goby_lag)


str(dat)
#View(dat)
nrow(dat) #n = 364


# use in dat.for.effects.plot script -----------------------
dat.for.effects.plot <- dat  

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
dat.for.effects.plot$Rain.unscaled <- rep(dat.temp$Rain_Sum)
dat.for.effects.plot$Breach.unscaled <- rep(dat.temp$BreachDays)

dat.for.effects.plot <- na.omit(dat.for.effects.plot) #, most missingness due to DO
nrow(dat.for.effects.plot) #314

####---------




#save a file with all cases including missing cases
dat.missing <- tibble(Goby=dat$Goby, 
                      Year=dat$Year, 
                      Year_2=dat$Year_2, 
                      Year_int=dat$Year_int, 
                      SAV=dat$SAV, 
                      SAV_2=dat$SAV_2,
                      SB=dat$SB,
                      SB_count=dat$SB_count,
                      SC=dat$SC, 
                      SC_count=dat$SC_count, 
                      Rain=dat$Rain, 
                      Temp=dat$Temp, 
                      Temp_2=dat$Temp_2, 
                      DO=dat$DO, 
                      Breach=dat$Breach, 
                      Breach_count=dat$Breach_count, 
                      BreachDays=dat$BreachDays,
                      BreachDays_2=dat$BreachDays_2,
                      BreachDays_Count=dat$BreachDays_Count,
                      BreachDays_Count_2=dat$BreachDays_Count_2,
                      Wind=dat$Wind, 
                      Micro=dat$Micro, #added 2024-01-23
                      Zone=as.factor(dat$Zone), 
                      Substrate=as.factor(dat$Substrate), 
                      Area=dat$Area,
                      Goby_lag=dat$Goby_lag)
dat.missing

#remove variables with NAs except Temp, DO, and SAV ok
dat.missing <- dat.missing %>%
  filter_at(vars(Year, SAV, SB_count, Goby, SC_count, SB, SC, Substrate, Wind, Zone, Micro, Area), 
            all_vars(!is.na(.)))

##remove cases with NAs
dat <- na.omit(dat) #, most missingness due to DO
nrow(dat) #314




dat.list <- as.list(dat)  # for stan


#need to convert substrate and zone to numbers


dat <- tibble(Goby=dat$Goby, 
              Year=dat$Year, 
              Year_2=dat$Year_2, 
              Year_int=dat$Year_int, 
              SAV=dat$SAV, 
              SAV_2=dat$SAV_2,
              SB=dat$SB,
              SB_count=dat$SB_count,
              SC=dat$SC, 
              SC_count=dat$SC_count, 
              Rain=dat$Rain, 
              Temp=dat$Temp, 
              Temp_2=dat$Temp_2,
              DO=dat$DO, 
              Breach=dat$Breach, 
              Breach_count=dat$Breach_count, 
              BreachDays = dat$BreachDays,
              BreachDays_2 = dat$BreachDays_2,
              BreachDays_Count = dat$BreachDays_Count,
              BreachDays_Count_2 = dat$BreachDays_Count_2,
              Wind=dat$Wind,
              Micro=dat$Micro,
              Zone=as.integer(dat$Zone), 
              Substrate=as.integer(dat$Substrate), 
              Area=dat$Area,
              Goby_lag=dat$Goby_lag)
dat
summary(dat$Area)
sd(dat$Area)
dat$Zone


dat |>
  ggplot(aes(x = SC_count/Area, y = Goby/Area)) +
  geom_point() +
  geom_smooth()








