# confidence sensitivity analysis for predicts pollinators
# packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "gt", "webshot", "lme4", "broom.mixed", "MASS")
lapply(packages, require, character.only = TRUE)

source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_5.rds")

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
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# create new column just for 1-4 for confidence
PREDICTS_pollinators$conf_summary <- PREDICTS_pollinators$confidence
PREDICTS_pollinators$conf_summary[PREDICTS_pollinators$confidence == 5.1] <- 1
PREDICTS_pollinators$conf_summary[PREDICTS_pollinators$conf_summary == 5.2] <- 2
PREDICTS_pollinators$conf_summary[PREDICTS_pollinators$conf_summary == 5.3] <- 3
PREDICTS_pollinators$conf_summary[PREDICTS_pollinators$conf_summary == 5.4] <- 4

# set up the vector for all the confidences to remove
confidence_vector <- list(NA, 1, c(1:2))

# create empty list for each iteration
jacked_plots <- list()

# succesivly remove the higher confidence levels
for(i in 1:length(confidence_vector)){
  
  # filter the pollinator data set for each set of UN_regions
  PREDICTS_pollinators_filt <- PREDICTS_pollinators[!PREDICTS_pollinators$conf_summary %in% confidence_vector[[i]],]
  
  # drop any leftover factor levels
  PREDICTS_pollinators_filt <- droplevels(PREDICTS_pollinators_filt)
  
  # print the row numbers after filtering the dataset
  nrow(PREDICTS_pollinators_filt)
  
  # calculate site metrics including all species (confirmed and not confirmed pollinator)
  pollinator_metrics_filt <- SiteMetrics(diversity = PREDICTS_pollinators_filt,
                                    extra.cols = c("SSB", "SSBS", "LUI", "UN_region"),
                                    sites.are.unique = TRUE, srEstimators = TRUE)

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
                                                                                "Intermediate secondary vegetation-Intense use",
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
                                      "Intermediate secondary vegetation-Intense use" = "ISVIU",
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
                                      coord_cartesian(ylim = c(-42, 80)) +
                                      xlab(NULL) +
                                      scale_y_continuous("Richness diff. (%)") +
                                      guides(shape = FALSE) + 
                                      theme(axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank())), 
                                   (abundance_metric + 
                                      coord_cartesian(ylim = c(-50, 180)) +
                                      xlab(NULL) +
                                      scale_y_continuous("Abundance diff. (%)") +
                                      guides(shape = FALSE) +
                                      theme(axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank())),
                                   ncol = 2)}
}

# set theme space to remove space between rows of plots
plot_title <- function(region){ggdraw() +
    draw_label(
      paste("Including", region, sep = " "),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(plot.margin = margin(60, 0, 0, 25))
}

# combine all the sets of plot for all jacknives
plot_grid(plot_grid(plot_title(region = "confidence 1, 2, 3, and 4"), jacked_plots[[1]], ncol = 1), 
          plot_grid(plot_title(region = "confidence 2, 3, and 4"), jacked_plots[[2]], ncol = 1), 
          plot_grid(plot_title(region = "confidence 3 and 4"), jacked_plots[[3]], ncol = 1), nrow = 3) 

ggsave("confidence_sensitivity_validation_exclude-high.png", scale = 1.5, dpi = 350)
