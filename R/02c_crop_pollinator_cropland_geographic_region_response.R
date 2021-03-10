### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "lme4", "patchwork", "gt", "broom.mixed", "MASS")
lapply(packages, require, character.only = TRUE)

# extra packages called in place -- plyr and StatisticalModels

source("R/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# bind together the intensity and type into a single variable
PREDICTS_pollinators$LUI <- paste(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity, sep = "-")

PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(LUI %in% c("Primary vegetation-Minimal use","Cropland-Minimal use", "Cropland-Light use", "Cropland-Intense use")) %>%
  dplyr::filter(Class %in% c("Mammalia", "Insecta", "Aves")) %>%
  dplyr::filter(Order %in% c( "Hymenoptera", "Diptera", "Coleoptera", 
                              "Thysanoptera", "Chiroptera", "Passeriformes", 
                              "Psittaciformes", "Apodiformes", "Columbiformes")) %>%
  dplyr::filter(!Family %in% c("Siricidae", "Tenthredinidae")) %>%
  droplevels()
  
# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# convert to factor and reorder factor levels
pollinator_metrics$LUI <- factor(pollinator_metrics$LUI, levels = c("Primary vegetation-Minimal use",
                                                                    "Cropland-Minimal use",
                                                                    "Cropland-Light use",
                                                                    "Cropland-Intense use"))
pollinator_metrics <- droplevels(pollinator_metrics)

# convert all factor levels to a short acronynm for spacing
pollinator_metrics <- pollinator_metrics %>%
  mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "Primary vegetation",
                                    "Cropland-Minimal use" = "Minimal use",
                                    "Cropland-Light use" = "Light use",
                                    "Cropland-Intense use" = "Intense use")))

# assign new variable for tropical/temperate, convert to factor, and filter out NA
pollinator_metrics$zone <- ifelse(pollinator_metrics$Latitude >= -23.5 & pollinator_metrics$Latitude <= 23.5, "Tropics", "Temperate")
pollinator_metrics$zone <- factor(pollinator_metrics$zone, levels = c("Temperate", "Tropics"))
pollinator_metrics <- pollinator_metrics %>%
  filter(!is.na(zone))

# multi panel plot of abundance, richness, and diversity
# add 1 for abundance and simpson diversity
pollinator_metrics$Total_abundance <- pollinator_metrics$Total_abundance + 1
pollinator_metrics$Simpson_diversity <- pollinator_metrics$Simpson_diversity + 1

model_1b <- glmer(Species_richness ~ LUI * zone + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612
summary(model_1b)
StatisticalModels::GLMEROverdispersion(model_1b)

# predict values for species richness
richness_metric <- predict_effects(iterations = 1000,
                                   model = model_1b, 
                                   model_data = pollinator_metrics, 
                                   response_variable = "Species_richness",
                                   fixed_number = 2,
                                   factor_number_1 = 2,
                                   factor_number_2 = 4,
                                   fixed_column = c("zone", "LUI"),
                                   neg_binom = FALSE)

model_2a <- lmer(log(Total_abundance) ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

# predict values from the model for abundance                               
abundance_metric <- predict_effects(iterations = 1000,
                                    model = model_2a, 
                                    model_data = pollinator_metrics, 
                                    response_variable = "Total_abundance",
                                    fixed_number = 2,
                                    factor_number_1 = 2,
                                    factor_number_2 = 4,
                                    fixed_column = c("zone", "LUI"),
                                    neg_binom = FALSE)

# combine the separate plots into a single figure
(richness_metric + xlab("") + ggtitle("A") + scale_x_discrete(labels = c("Non-tropical", "Tropical"))) + 
  (abundance_metric + ggtitle("B") + xlab("") + scale_x_discrete(labels = c("Non-tropical", "Tropical")) + guides(colour = FALSE)) + 
  plot_layout(ncol = 1) & 
  theme(text = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        aspect.ratio = 0.7/2)

ggsave("crop_pollinators_resampling_cropland_zone_8.png", dpi = 400, scale = 1)
