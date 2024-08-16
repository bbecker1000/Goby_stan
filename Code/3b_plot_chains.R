#plot chains

chain.data <- as_tibble(extract.samples(fit))

#remove warmup
chain.data <- chain.data[c(3001:6000),]
#add rownames
chain.data$sample <- 1:3000

MEAN <- chain.data %>% summarize(MEAN = mean(beta_BreachDays))
chain.data %>%
  ggplot(aes(x=sample, y=beta_BreachDays)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Wind))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Wind)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_SAV))
chain.data %>%
  ggplot(aes(x=sample, y=beta_SAV)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_DO))
chain.data %>%
  ggplot(aes(x=sample, y=beta_DO)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Micro))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Micro)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Substrate))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Substrate)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Goby_lag))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Goby_lag)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_SB_count))
chain.data %>%
  ggplot(aes(x=sample, y=beta_SB_count)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_SC_count))
chain.data %>%
  ggplot(aes(x=sample, y=beta_SC_count)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)



#acf plots
par(mfrow=c(3,4))

acf(chain.data$beta_Substrate)
acf(chain.data$beta_BreachDays)
acf(chain.data$beta_SB_count)
acf(chain.data$beta_SC_count)
acf(chain.data$beta_Goby_lag)
acf(chain.data$beta_Micro)
acf(chain.data$beta_DO)
acf(chain.data$beta_Temp)
acf(chain.data$beta_SAV)
acf(chain.data$beta_Wind)

par(mfrow=c(1,1))

#Pairs plots
cowplot::plot_grid(
  
  chain.data %>%
  ggplot(aes(x=beta_BreachDays, y=beta_SC_count)) +
    geom_point(color = "blue",alpha = 0.2)
,

chain.data %>%
  ggplot(aes(x=beta_SB_count, y=beta_Goby_lag)) +
    geom_point(color = "blue",alpha = 0.2)
,

chain.data %>%
  ggplot(aes(x=beta_BreachDays, y=beta_SC_count)) +
    geom_point(color = "blue",alpha = 0.2)
,

chain.data %>%
  ggplot(aes(x=phi, y=tau)) +
    geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_Micro, y=beta_BreachDays)) +
    geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_Wind, y=beta_BreachDays)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_Year, y=beta_BreachDays)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_Year, y=beta_Rain)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=mu_Zone, y=tau_Zone)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_SAV, y=beta_Substrate)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_DO, y=beta_Temp)) +
  geom_point(color = "blue",alpha = 0.2)
,
chain.data %>%
  ggplot(aes(x=beta_Goby_lag, y=beta_Year)) +
  geom_point(color = "blue",alpha = 0.2)

)










