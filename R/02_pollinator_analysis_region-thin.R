# jack-knife the model by UN_region to check that particular areas aren't having a significant effect on the model

# packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "patchwork", "gt", "webshot", "lme4", "broom.mixed", "MASS")
lapply(packages, require, character.only = TRUE)

source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# check use_intensity and use_type factors
table(PREDICTS_pollinators$Use_intensity)
table(PREDICTS_pollinators$Predominant_land_use)
table(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity)

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

# set up the string for all the UN_regions
UN_regions <- levels(PREDICTS_pollinators$UN_region)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI", "UN_region"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# create empty list for each iteration
jacked_plots <- list()

# for each UN_region calculate the model downstream
for(i in 1:length(UN_regions)){

  # create the indexed UN_region vector for that iteration
  UN_regions_filt <- UN_regions[-i]
  
  # filter the pollinator data set for each set of UN_regions
  pollinator_metrics_filt <- pollinator_metrics[pollinator_metrics$UN_region %in% UN_regions_filt,]
  
  # drop any leftover factor levels
  pollinator_metrics_filt <- droplevels(pollinator_metrics_filt)
  
  # print the row number of the new dataframe
  print(nrow(pollinator_metrics_filt))

  # check that the same number of factor combinations are all represented for LUI -
  print(length(unique(pollinator_metrics_filt$LUI))) # all confirmed as 33
  
  # convert to factor and reorder factor levels
  pollinator_metrics_filt$LUI <- factor(pollinator_metrics_filt$LUI, levels = c("Primary vegetation-Minimal use", 
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
  # drop the extra levels
  pollinator_metrics_filt <- droplevels(pollinator_metrics_filt)
  
  # convert all factor levels to a short acronynm for spacing
  pollinator_metrics_filt <- pollinator_metrics_filt %>%
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
  
  # multi panel plot of abundance, richness, and diversity
  # add 1 for abundance and simpson diversity
  pollinator_metrics_filt$Total_abundance <- pollinator_metrics_filt$Total_abundance + 1
  pollinator_metrics_filt$Simpson_diversity <- pollinator_metrics_filt$Simpson_diversity + 1
  
  ### richness - because only 1 fixed effect, select model random effect structure on basis of AIC
  model_1b <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics_filt, family = poisson) # best model - failed to converge with 0.00144612
  richness_metric <- predict_effects(iterations = 1000, 
                                     model = model_1b, 
                                     model_data = pollinator_metrics_filt, 
                                     response_variable = "Species_richness",
                                     fixed_number = 1,
                                     fixed_column = "LUI",
                                     neg_binom = FALSE)
  
  ### abundance
  model_2a <- lmer(log(Total_abundance) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics_filt) # best model
  abundance_metric <- predict_effects(iterations = 1000, 
                                      model = model_2a, 
                                      model_data = pollinator_metrics_filt, 
                                      response_variable = "Total_abundance",
                                      fixed_number = 1,
                                      fixed_column = "LUI",
                                      neg_binom = FALSE)
  
  ### diversity
  model_3a <- lmer(log(Simpson_diversity) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics_filt) # best model
  diversity_metric <- predict_effects(iterations = 1000, 
                                      model = model_3a, 
                                      model_data = pollinator_metrics_filt, 
                                      response_variable = "Simpson_diversity",
                                      fixed_number = 1,
                                      fixed_column = "LUI",
                                      neg_binom = FALSE)
  
  # build the plot for each jackknife
  if(i <= 5){
    jacked_plots[[i]] <- plot_grid((richness_metric + 
                                      coord_cartesian(ylim = c(-42, 110)) +
                                      xlab(NULL) +
                                      ggtitle(paste("Excluding", UN_regions[i])) +
                                      scale_y_continuous("Richness diff. (%)") +
                                      guides(shape = FALSE) + 
                                      theme(axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank())),
              (abundance_metric + 
                 coord_cartesian(ylim = c(-50, 250)) +
                 xlab(NULL) +
                 ggtitle(" ") +
                 scale_y_continuous("Abundance diff. (%)") +
                 guides(shape = FALSE) +
                 theme(axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank())), ncol = 2)

  }
}

{(jacked_plots[[1]] + jacked_plots[[2]] + jacked_plots[[3]] + jacked_plots[[4]] + jacked_plots[[5]])} + plot_layout(ncol = 1)

ggsave("UN_region_jack.png", dpi = 350, scale = 1.5)
