# script for counting sites and studies

# package to read in 
library(dplyr)

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# calculating numbers of sources, sites, and studies
length(unique(PREDICTS_pollinators$Source_ID))
length(unique(PREDICTS_pollinators$SSBS))
length(unique(PREDICTS_pollinators$SS))

# remove additional columns - number of unique species = 5100
pollinator_subset <- PREDICTS_pollinators %>% 
  select(Kingdom, Phylum, Class, Order, Family, Best_guess_binomial) %>%
  filter(Best_guess_binomial != "") %>%
  filter(!grepl("\\d", Best_guess_binomial)) %>%
  unique()

# number of species total
pollinator_subset %>%
  group_by(Order, Class) %>%
  tally() %>% 
  ungroup() %>%
  summarise(total = sum(n))

# number pf species for each class
pollinator_subset %>%
  group_by(Class) %>%
  tally()

# number of vertebrate species
pollinator_subset %>%
  filter(Class %in% c("Insecta")) %>%
  tally()

# number of invertebrate speces
pollinator_subset %>%
  filter(Class %in% c("Aves", "Mammalia", "Reptilia")) %>%
  tally()
