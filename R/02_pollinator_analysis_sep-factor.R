### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "lme4",  "MASS")
lapply(packages, require, character.only = TRUE)

# extra packages called in place -- plyr and StatisticalModels

source("R/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# filter cannot decide factors
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use != "Cannot decide") %>%
  dplyr::filter(Use_intensity != "Cannot decide") %>%
  dplyr::filter(Predominant_land_use != "Secondary vegetation (indeterminate age)") %>%
  droplevels()

# bind together the intensity and type into a single variable
PREDICTS_pollinators$LUI <- paste(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity, sep = "-")

# filter out the factors with low representation, and filter for just minimal use to test effect of land-use type
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(LUI != "Mature secondary vegetation-Intense use") %>%
  dplyr::filter(LUI != "Intermediate secondary vegetation-Intense use") %>%
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI", "Use_intensity", "Predominant_land_use"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# convert to factor and reorder factor levels
pollinator_metrics <- droplevels(pollinator_metrics)

# multi panel plot of abundance, richness, and diversity
# add 1 for abundance and simpson diversity
pollinator_metrics$Total_abundance <- pollinator_metrics$Total_abundance + 1
pollinator_metrics$Simpson_diversity <- pollinator_metrics$Simpson_diversity + 1

### richness - because only 1 fixed effect, select model random effect structure on basis of AIC
model_1b <- glmer(Species_richness ~ Predominant_land_use * Use_intensity + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612
summary(model_1b)
anova(model_1b, test="Chisq")

### abundance
# lmer with log transformation
model_2a <- lmerTest::lmer(log(Total_abundance) ~ Predominant_land_use * Use_intensity + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
summary(model_2a)
anova(model_2a)

### diversity
model_3a <- lmerTest::lmer(log(Simpson_diversity) ~ Predominant_land_use * Use_intensity + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

summary(model_3a)
anova(model_3a)

