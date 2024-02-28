#adding 2023 data

library(readxl)
library(tidyverse)

counts_2023 <- read_excel("C:/projects/Goby/data/Goby2023data_Rodeo.xlsx", 
                          sheet = "Fishdat")

d1 <- counts_2023 %>% 
  group_by(Station, SPECIES) %>%
  summarize(Count = sum(Numbers))



counts_2023 <- counts_2023 %>%
  mutate(microsporidium = ifelse(`Microsporidian Cyst/tumor` == "N", 0, 1))


d1.micro_sum <- counts_2023 %>% 
                filter(SPECIES == "TW") %>%
                filter(!is.na(microsporidium)) %>%
                group_by(Station) %>%
                summarise_at(vars(microsporidium),
                             list(micro_sum = sum))
  





            