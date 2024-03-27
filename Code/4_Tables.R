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