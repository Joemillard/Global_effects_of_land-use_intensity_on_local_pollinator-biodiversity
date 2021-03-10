# script for asymptotic species richness as a response variable

### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "lme4", "MASS", "raster", "spdep")
lapply(packages, require, character.only = TRUE)

# source in scripts for deriving predictions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

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

# filter out the factors with low representation
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(LUI != "Mature secondary vegetation-Intense use") %>%
  dplyr::filter(LUI != "Intermediate secondary vegetation-Intense use") %>%
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI"),
                                  sites.are.unique = TRUE,
                                  srEstimators = "Chao")

# convert to factor and reorder factor levels
pollinator_metrics$LUI <- factor(pollinator_metrics$LUI, levels = c("Primary vegetation-Minimal use", 
                                                                    "Primary vegetation-Light use", 
                                                                    "Primary vegetation-Intense use", 
                                                                    "Mature secondary vegetation-Minimal use", 
                                                                    "Mature secondary vegetation-Light use", 
                                                                    "Intermediate secondary vegetation-Minimal use",
                                                                    "Intermediate secondary vegetation-Light use",
                                                                    "Young secondary vegetation-Minimal use",
                                                                    "Young secondary vegetation-Light use",
                                                                    "Young secondary vegetation-Intense use",
                                                                    "Plantation forest-Minimal use",
                                                                    "Plantation forest-Light use",
                                                                    "Plantation forest-Intense use",
                                                                    "Pasture-Minimal use", 
                                                                    "Pasture-Light use", 
                                                                    "Pasture-Intense use", 
                                                                    "Cropland-Minimal use", 
                                                                    "Cropland-Light use", 
                                                                    "Cropland-Intense use",
                                                                    "Urban-Minimal use",
                                                                    "Urban-Light use",
                                                                    "Urban-Intense use"))
pollinator_metrics <- droplevels(pollinator_metrics)

# convert all factor levels to a short acronynm for spacing
pollinator_metrics <- pollinator_metrics %>%
  dplyr::mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "PVMU",
                                           "Primary vegetation-Light use" = "PVLU",
                                           "Primary vegetation-Intense use" = "PVIU",
                                           "Young secondary vegetation-Minimal use" = "YSVMU", 
                                           "Young secondary vegetation-Light use" = "YSVLU",
                                           "Young secondary vegetation-Intense use" = "YSVIU", 
                                           "Intermediate secondary vegetation-Minimal use" = "ISVMU", 
                                           "Intermediate secondary vegetation-Light use" = "ISVLU",
                                           "Mature secondary vegetation-Minimal use" = "MSVMU", 
                                           "Mature secondary vegetation-Light use" = "MSVLU", 
                                           "Plantation forest-Minimal use" = "PFMU",
                                           "Plantation forest-Light use" = "PFLU",
                                           "Plantation forest-Intense use" = "PFIU",
                                           "Pasture-Minimal use" = "PMU", 
                                           "Pasture-Light use" = "PLU", 
                                           "Pasture-Intense use" = "PIU", 
                                           "Cropland-Minimal use" = "CMU", 
                                           "Cropland-Light use" = "CLU", 
                                           "Cropland-Intense use" = "CIU", 
                                           "Urban-Minimal use" = "UMU",
                                           "Urban-Light use" = "ULU",
                                           "Urban-Intense use" = "UIU")))

# multi panel plot of abundance, richness, and diversity
# add 1 for abundance and simpson diversity
pollinator_metrics$Total_abundance <- pollinator_metrics$Total_abundance + 1
pollinator_metrics$Simpson_diversity <- pollinator_metrics$Simpson_diversity + 1
pollinator_metrics$ChaoR <- round(pollinator_metrics$ChaoR)

model_1b <- glmer(ChaoR ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612

richness_metric <- predict_effects(iterations = 1000, 
                                   model = model_1b, 
                                   model_data = pollinator_metrics, 
                                   response_variable = "ChaoR",
                                   fixed_number = 1,
                                   fixed_column = "LUI",
                                   neg_binom = FALSE)

ggsave("Chao_richness_supp_figure.png", scale = 1, dpi = 350)
