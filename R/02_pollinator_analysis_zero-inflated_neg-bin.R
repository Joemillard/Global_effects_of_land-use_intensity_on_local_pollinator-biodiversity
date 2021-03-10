# pollinator analysis - zero inflated negative binomial

### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "gt", "webshot", "lme4", "broom.mixed", "MASS", "roquefort")
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

# check the distribution of factors for intensity and type
table(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity)

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
                                  extra.cols = c("SSB", "SSBS", "LUI", "UN_region", "Country", "Source_for_predominant_land_use"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

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
  mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "PVMU",
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

# richness glmer model
model_1b <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) 

# derive predictions for species richness
richness_metric <- predict_effects(iterations = 1000, 
                                   model = model_1b, 
                                   model_data = pollinator_metrics, 
                                   response_variable = "Species_richness",
                                   fixed_number = 1,
                                   fixed_column = "LUI",
                                   neg_binom = FALSE)

# total abundacne glmmTMB model
model_2a <- glmmTMB::glmmTMB(Total_abundance ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics, ziformula = ~1, family = list(family = "nbinom2", link = "log"))

# derive predictions for abundance
abundance_metric <- predict_effects(iterations = 1000, 
                                    model = model_2a, 
                                    model_data = pollinator_metrics, 
                                    response_variable = "Total_abundance",
                                    fixed_number = 1,
                                    fixed_column = "LUI",
                                    neg_binom = TRUE)

# Simpson diversity glmmTMB model
model_3a <- glmmTMB::glmmTMB(Simpson_diversity ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics, ziformula = ~0, family = list(family = "nbinom2", link = "log")) 

 # derive predictions for simpson diversity
diversity_metric <- predict_effects(iterations = 1000, 
                                    model = model_3a, 
                                    model_data = pollinator_metrics, 
                                    response_variable = "Simpson_diversity",
                                    fixed_number = 1,
                                    fixed_column = "LUI", 
                                    neg_binom = TRUE)

# combined the plots into a single figure
plot_grid((richness_metric + ggtitle("A")), 
          (abundance_metric + guides(shape = FALSE) + ggtitle("B")), 
          (diversity_metric + guides(shape = FALSE) + ggtitle("C")),
          NULL,
          ncol = 2)

ggsave("all_metrics_resampling_approach_all-glmer_2.png", dpi = 400, scale = 1.6)
