## read in the taxonomic references and convert to just the list of references

# load in packages
library(dplyr)

# clean the taxonomic reference data
taxa_references <- read.csv("Data/global clade extrapolation/clade_extrapolation.csv") %>%
  select(additional_citations) %>%
  rename("Taxonomic_references" = "additional_citations") %>%
  filter(Taxonomic_references != "") %>%
  mutate(Taxonomic_references = gsub("Diptera; Culicidae", "Diptera, Culicidae", Taxonomic_references)) %>%
  tidyr::separate_rows(Taxonomic_references, sep = ";") %>%
  unique() %>%
  dplyr::filter(!grepl("Michael", Taxonomic_references)) %>%
  arrange(Taxonomic_references)

# write the taxonomic reference data to csv
write.csv(taxa_references, "pollinat_taxonomic_references.csv")
