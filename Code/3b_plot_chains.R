#plot chains

chain.data <- as_tibble(extract.samples(fit))

#remove warmup
chain.data <- chain.data[c(3001:6000),]
#add rownames
chain.data$sample <- 1:3000

chain.data %>%
  ggplot(aes(x=sample, y=beta_BreachDays)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_Wind)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_SAV)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_DO)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_Micro)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_Substrate)) +
  geom_line(alpha = 0.2)

chain.data %>%
  ggplot(aes(x=sample, y=beta_Goby_lag)) +
  geom_line(alpha = 0.2)
