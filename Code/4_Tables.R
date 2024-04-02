#tables


####summary data table

names(dat.temp)

#get same values as in analysis table dat in script 1
d.table1 <- dat.temp %>% 
  filter_at(vars(Year, SAV, Sum_TW, Sum_SB, Sum_SC, u_mean, micro_sum, min_DO,Area, temp_mean, Dom_substrate), 
          all_vars(!is.na(.)))


# median and range
#filter cases that are not complete to get same as working dataset.

MEANS <- d.table1 %>%
  summarize_at(c("Year", "Sum_TW", "BreachDays", "micro_sum", "SAV", "Sum_SB", "Sum_SC", 
                 "Rain_Sum", "min_DO", "u_mean", "temp_mean", "Area"), mean, na.rm = TRUE)

RANGES <- dat.temp %>%
  summarize_at(c("Year", "Sum_TW", "BreachDays", "micro_sum", "SAV", "Sum_SB", "Sum_SC", 
                 "Rain_Sum", "min_DO", "u_mean", "temp_mean", "Area"), range, na.rm = TRUE)

SUMMARY_TABLE <- dplyr::bind_rows(MEANS, RANGES)

SUMMARY_TABLE <-t(SUMMARY_TABLE)
#send to google sheet


#model output table
#  coefficients, rhat, ESS

# done in google sheet from model output in 2_Full_luxury....
#  causal coefficients

#model statistics table

#estimated density at beginning, middle and end of time series.
#data from 3_plotting

mean(dat$Goby)
dat$Area
mean(exp(dat$Area))
mean(dat$Goby) / mean(exp(dat$Area))

unique(dat$Year)

density_mean <- d %>% 
  group_by(Year_int) %>%
  summarize(MEAN = mean(mu/exp(Area)))

density_LOW <- d %>% 
  group_by(Year_int) %>%
  summarize(LOW = mean(mu/exp(Area)))

density_HI <- d %>% 
  group_by(Year_int) %>%
  summarize(LOW = mean(mu/exp(Area)))

library(meantables)

MEAN_TABLE <- d %>%
  mutate(GobyDens = mu/exp(Area)) %>%
  group_by(Year_int) %>%
  mean_table(GobyDens)

MEAN_TABLE$year <- MEAN_TABLE$group_cat+1995

d %>%
  ggplot(aes(Year_int+1995, mu/exp(Area),color = as.factor(Zone))) +
  geom_jitter(alpha = 0.1) +
  #geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd), color = "blue") +
  ylab("Goby Density") +
  geom_smooth (method = "lm", se=TRUE,
               formula = y~poly(x,2),
               alpha=0., linewidth=0.75) +
  theme_few(base_size = 16)




