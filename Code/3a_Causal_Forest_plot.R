
#causal plot
library(ggthemes)

post <- extract.samples(fit)
PROBS = c(0.055, 0.5, 0.945) ##89 % CIs
#Effect of Rain --> Breach --> DO --> Goby
Rain_Breach_DO<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays*beta_DO), probs = PROBS)) ## NS
#Effect of Breach --> DO --> Goby
Breach_DO<-as_tibble(quantile( with(post, beta_BreachDays*beta_DO), probs = PROBS)) ##  NS
#Effect of Wind --> Breach --> DO --> Goby
Wind_DO<-as_tibble(quantile( with(post, beta_Wind*beta_DO), probs = PROBS)) ## ns
#Effect of Wind --> Temp --> DO --> Goby
Wind_Temp_DO<-as_tibble(quantile( with(post, beta_Wind*beta_Temp*beta_DO), probs = PROBS)) ## ns  
#Effect of Wind --> Temp --> SAV --> Goby
Wind_DO_SAV<-as_tibble(quantile( with(post, beta_Wind*beta_SAV*beta_DO), probs = PROBS)) ## ns  
#Effect of Wind --> DO --> SB --> Goby
Wind_DO_SC<-as_tibble(quantile( with(post, beta_Wind*beta_DO*beta_SC_count), probs = PROBS)) ## ns 
#Effect of DO --> SB --> Goby
DO_SB<-as_tibble(quantile( with(post, beta_DO*beta_SB_count), probs = PROBS)) ## ns 

DO<-as_tibble(quantile( with(post, beta_DO), probs = PROBS)) ## ns 

Rain_Breach_Temp<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays*beta_Temp), probs = PROBS)) ## NS
#Effect of Breach --> Temp --> Goby
Rain_Breach_2<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays_2), probs = PROBS)) ## NS
#Effect of Rain --> Breach --> Goby

Rain_Breach<- as_tibble(quantile( with(post, beta_Rain*beta_BreachDays), probs = PROBS)) ## NS
#Effect of Rain --> Breach --> Goby

Breach_Temp<-as_tibble(quantile( with(post, beta_BreachDays*beta_Temp), probs = PROBS)) ##  NS
#Effect of Wind --> Breach --> Temp --> Goby
Wind_Temp<-as_tibble(quantile( with(post, beta_Wind*beta_Temp), probs = PROBS)) ## ns

#Effect of Wind --> Temp --> SAV --> Goby
Wind_Temp_SAV<-as_tibble(quantile( with(post, beta_Wind*beta_SAV*beta_Temp), probs = PROBS)) ## ns  
#Effect of Wind --> Temp --> SB --> Goby
Wind_Temp_SC<-as_tibble(quantile( with(post, beta_Wind*beta_Temp*beta_SC_count), probs = PROBS)) ## ns 
#Effect of Temp --> SB --> Goby
Temp_SB<-as_tibble(quantile( with(post, beta_Temp*beta_SB_count), probs = PROBS)) ## ns 
Rain_Temp_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_Temp*beta_SAV*beta_SB_count), probs = PROBS))

Temp<-as_tibble(quantile( with(post, beta_Temp), probs = PROBS)) ## ns 
Temp_2<-as_tibble(quantile( with(post, beta_Temp_2), probs = PROBS)) ## ns 

Breach_2<-as_tibble(quantile( with(post, beta_BreachDays_2), probs = PROBS)) ## ns 

Year<-as_tibble(quantile( with(post, beta_Year), probs = PROBS)) ## ns 
Year_2<-as_tibble(quantile( with(post, beta_Year_2), probs = PROBS)) ## ns 


Substrate_SC<-as_tibble(quantile( with(post, beta_Substrate*beta_SC_count), probs = PROBS)) ## negative
Substrate<-as_tibble(quantile( with(post, beta_Substrate), probs = PROBS))
Rain_DO_SAV_SB<-as_tibble(quantile( with(post, beta_Rain*beta_DO*beta_SAV*beta_SB_count), probs = PROBS))

Breach<-as_tibble(quantile( with(post, beta_BreachDays), probs = PROBS)) #positive
Rain_Micro<-as_tibble(quantile( with(post, beta_Rain*beta_Micro), probs = PROBS))
Micro<-as_tibble(quantile( with(post, beta_Micro), probs = PROBS))
SAV_SB<-as_tibble(quantile( with(post, beta_SAV*beta_SB_count), probs = PROBS))
SAV_SC<-as_tibble(quantile( with(post, beta_SAV*beta_SC_count), probs = PROBS))
SAV<-as_tibble(quantile( with(post, beta_SAV), probs = PROBS))
SAV_2<-as_tibble(quantile( with(post, beta_SAV_2), probs = PROBS))
SC<-as_tibble(quantile( with(post, beta_SC_count), probs = PROBS))
SB<-as_tibble(quantile( with(post, beta_SB_count), probs = PROBS))
Breach_DO_SC<-as_tibble(quantile( with(post, beta_BreachDays*beta_DO*beta_SC_count), probs = PROBS)) ##  NS
Goby_Lag<-as_tibble(quantile( with(post, beta_Goby_lag), probs = PROBS))
# 2024-02-01
# ADD RAIN --> BREACH



#names for tibble
names<- c("RAIN -> Breach -> DO", "RAIN -> Breach^2", "BREACH -> DO", "WIND -> DO", "WIND -> Temp -> DO", "WIND -> DO -> SAV", "WIND -> DO -> SC",
          "DO -> SB", 
          "YEAR", 
          "SUBSTRATE -> SC", "SUBSTRATE", "RAIN -> DO -> SAV -> SB", "BREACH", 
          #"DO -> SB", "Year", "Substrate -> SC", "Rain -> DO -> SAV -> SB", "Breach", 
          "MICRO", "SAV -> SB", "SAV -> SC", "SC", "SB","BREACH -> DO -> SC", 
          "DO",
          "RAIN -> Breach -> Temp",
          #"Zone", "Zone -> Wind -> int", 
          "BREACH -> Temp",
          "WIND -> Temp",
          "WIND -> Temp -> SAV",
          "WIND -> Temp -> SC",
          "TEMP -> SB",
          "RAIN -> Temp -> SAV -> SB",
          "TEMP",
          "TEMP^2",
          "BREACH^2",
          "RAIN -> Breach",
          "YEAR^2",
          "SAV",
          "SAV^2",
          "Goby_Lag")
#add probabilities
plot.posteriors<-rbind(Rain_Breach_DO, Rain_Breach_2, Breach_DO, Wind_DO, Wind_Temp_DO, Wind_DO_SAV, Wind_DO_SC,
                       DO_SB, 
                       Year, 
                       Substrate_SC, Substrate, Rain_DO_SAV_SB, Breach, 
                       # DO_SB, Year, Substrate_SC, Rain_DO_SAV_SB, Breach, 
                       Micro, SAV_SB, SAV_SC, SC, SB, Breach_DO_SC, 
                       DO,
                       Rain_Breach_Temp,
                       #Zone, Zone_Wind_int,
                       Breach_Temp,
                       Wind_Temp,
                       Wind_Temp_SAV,
                       Wind_Temp_SC,
                       Temp_SB,
                       Rain_Temp_SAV_SB,
                       Temp,
                       Temp_2,
                       Breach_2,
                       Rain_Breach,
                       Year_2,
                       SAV,
                       SAV_2,
                       Goby_Lag)
#add names
plot.posteriors$names <- rep(names, each=3)
#add probabilities names
plot.posteriors$probability <- rep(c("lower", "median", "upper"), times = length(names))

plot.posteriors.wide <- plot.posteriors %>% 
  group_by(probability) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = probability, values_from = value) %>%
  select(-row)

#add codes for positive or negative coefficients
plot.posteriors.wide$effect <- ifelse(
  plot.posteriors.wide$lower<0 & 
    plot.posteriors.wide$median <0 & 
    plot.posteriors.wide$upper <0, 
  "negative", ifelse(plot.posteriors.wide$lower>0 
                     & plot.posteriors.wide$median >0 & 
                       plot.posteriors.wide$upper >0, 
                     "positive",
                     "neutral"))


print(plot.posteriors.wide, n = 31)

COLORS = c("red", "black", "blue")

plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
  ggplot(aes(x = names, y = median, color = effect)) +
  #geom_point(effect = c("red", "black", "blue"))) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Causal Path") + 
  ylab("Causal Effect on Goby Density") +
  scale_color_manual(breaks = c("negative", "neutral", "positive"),
                     values=c("red", "darkgray", "green3")) + 
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme_few(base_size = 16)

ggsave("Output/forest.plot.Yearlag.png", width = 20, height = 30, units = "cm")
str(rethinking::extract.samples(fit))

##2025-06-06
# add the priors to the plot
# new file made from exported posteriors and adding priors
library(readr)
ForestWithPriors <- read_csv("Data/ForestWithPriors.csv")
head(ForestWithPriors)

#plot with priors

#(fix priors in excel file)
ForestWithPriors[29,6] <- -0.1
ForestWithPriors[30,6] <- -0.1
ForestWithPriors[34,6] <- -0.1

ForestWithPriors[29,8] <- -0.35
ForestWithPriors[30,8] <- -0.35
ForestWithPriors[34,8] <- -0.35

ForestWithPriors[29,9] <- 0.15
ForestWithPriors[30,9] <- 0.15
ForestWithPriors[34,9] <- 0.15
#remove two rows of NAs
ForestWithPriors <- ForestWithPriors[-c(36:37),]



ForestWithPriors %>% # plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
  ggplot(aes(x = names, y = median, color = effect)) +
  #geom_point(effect = c("red", "black", "blue"))) +

  geom_pointrange(aes(ymin = lower, ymax = upper),position = position_nudge(x = 0.1)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(data = ForestWithPriors, 
                  aes(x = names, y = prior, 
                      ymin = prior_lo, ymax = prior_hi), 
                  color = "blue", shape = "square", alpha = 0.3, 
                      show.legend = FALSE,
                    position = position_nudge(x = -0.1)) +

  # geom_pointrange(aes(ymin = lower, ymax = upper)) +
  # geom_hline(yintercept = 0, lty = 2) +
  # geom_pointrange(data = ForestWithPriors, 
  #                 aes(x = names, y = prior, 
  #                     ymin = prior_lo, ymax = prior_hi), color = "lightblue") +

  xlab("Causal Path") + 
  ylab("Causal Effect on Goby Density") +
  labs(color='Posterior') +
  scale_color_manual(breaks = c("negative", "neutral", "positive"),
                     values=c("red", "darkgray", "green3")) + 
  coord_flip() +
  scale_x_discrete(limits=rev) +
  geom_point(aes(x = 14, y = 1, prior = 0.1), color = "lightblue", 
             shape = "square", size = 2.5) +
  annotate("text", x = 15, y = 1.15, label = "Prior", size = 6) +
  theme_few(base_size = 16) + 
  theme(legend.position = c(.8, .6))
ggsave("Output/forest.plot.Yearlag.priors.png", width = 20, height = 30, units = "cm")

  theme_few(base_size = 16)






