# analysis in cropland between temperate and tropical, split between low adn high forest cover baseline

### packages for analysis,and suppress messages
# extra packages called in place -- plyr and StatisticalModels
packages <- c("dplyr", "yarg", "lme4", "cowplot", "raster", "ggplot2", "MASS")
lapply(packages, require, character.only = TRUE)

source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# bring in the forest data and merge with the PREDICTS sites
## calculate forest extent for each predicts primary vegetation site

# read in the forest data
hansen_tree_cover <- raster(here::here("Data/forest_data/Hansen_full.tif"))

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# bind together the intensity and type into a single variable
PREDICTS_pollinators$LUI <- paste(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity, sep = "-")

PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(LUI %in% c("Primary vegetation-Minimal use","Cropland-Minimal use", "Cropland-Light use", "Cropland-Intense use")) %>%	 
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE) %>%
  filter(!is.na(Longitude))

# convert the coordinates for primary site to sp points
convert_spat <- function(data_file) {
  data_fin <- data_file %>%
    dplyr::select(Longitude, Latitude) %>%
    SpatialPoints()
  
  return(data_fin)
}

# run function for converting coordinates to spatial
prim_spat <- convert_spat(pollinator_metrics)

# extract the value of forest cover for each primary vegetation site
prim_cover <- extract(hansen_tree_cover, prim_spat, na.rm = TRUE)

# bind the coordinates back onto the extracted coordinates
pollinator_metrics_cover <- pollinator_metrics %>%
  cbind(prim_cover) %>%
  rename(forest_cover = prim_cover)

# create factors for high and low forest cover
pollinator_metrics_cover$forest_fact[pollinator_metrics_cover$forest_cover >= 60] <- "high_cover"
pollinator_metrics_cover$forest_fact[pollinator_metrics_cover$forest_cover <= 40] <- "low_cover"

# assign new variable for tropical/temperate, convert to factor, and filter out NA
pollinator_metrics_cover$zone <- ifelse(pollinator_metrics_cover$Latitude >= -23.5 & pollinator_metrics_cover$Latitude <= 23.5, "Tropical", "Non-tropical")
pollinator_metrics_cover$zone <- factor(pollinator_metrics_cover$zone, levels = c("Non-tropical", "Tropical"))
pollinator_metrics_cover <- pollinator_metrics_cover %>%
  filter(!is.na(zone))

# convert to factor and reorder factor levels
pollinator_metrics_cover$LUI <- factor(pollinator_metrics_cover$LUI, levels = c("Primary vegetation-Minimal use",
                                                                    "Cropland-Minimal use",
                                                                    "Cropland-Light use",
                                                                    "Cropland-Intense use"))
pollinator_metrics_cover <- droplevels(pollinator_metrics_cover)

# convert all factor levels to a short acronynm for spacing
pollinator_metrics_cover <- pollinator_metrics_cover %>%
  mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "Primary vegetation",
                                    "Cropland-Minimal use" = "Minimal use",
                                    "Cropland-Light use" = "Light use",
                                    "Cropland-Intense use" = "Intense use")))

# add 1 for abundance and simpson diversity
pollinator_metrics_cover$Total_abundance <- pollinator_metrics_cover$Total_abundance + 1
pollinator_metrics_cover$Simpson_diversity <- pollinator_metrics_cover$Simpson_diversity + 1

# filter for just cropland data set to bind onto primary vegetation
pollinator_metrics_cropland <- pollinator_metrics_cover %>%
  filter(LUI %in% c("Minimal use", "Light use", "Intense use")) %>%
  droplevels()

# create vector for baseline forest cover
forest_cover <- c("low_cover", "high_cover")

# set up plot objects for each metric
richness_object <- list()
abundance_object <- list()
diversity_object <- list()

# set up loop, with each iteration removing one forest cover baseline
for(i in 1:length(forest_cover)){
  
  # filter for primary vegetation high or low cover and drop the leftover levels
  pollinator_metrics_cover_filt <- pollinator_metrics_cover %>%
    filter(LUI == "Primary vegetation") %>%
    filter(forest_fact == !!forest_cover[i]) %>%
    droplevels()

  # bind the filtered primary vegetation data onto the cropland data
  pollinator_metrics_cover_filt <- rbind(pollinator_metrics_cover_filt, pollinator_metrics_cropland)
  
  # print the number of factor combinations
  print(table(pollinator_metrics_cover_filt$zone, pollinator_metrics_cover_filt$LUI))
  
  model_1b <- glmer(Species_richness ~ LUI * zone + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics_cover_filt, family = poisson) # best model - failed to converge with 0.00144612
    
  richness_object[[i]] <- predict_effects(iterations = 1000,
                                      model = model_1b, 
                                      model_data = pollinator_metrics_cover_filt, 
                                      response_variable = "Species_richness",
                                      fixed_number = 2,
                                      factor_number_1 = 2,
                                      factor_number_2 = 4,
                                      fixed_column = c("zone", "LUI"), 
                                      neg_binom = FALSE)
  
  model_2a <- lmer(log(Total_abundance) ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics_cover_filt) # best model
    
  # predict values from the model                                
  abundance_object[[i]] <- predict_effects(iterations = 1000,
                                           model = model_2a,                                              
                                           model_data = pollinator_metrics_cover_filt, 
                                           response_variable = "Total_abundance",
                                           fixed_number = 2,
                                           factor_number_1 = 2,
                                           factor_number_2 = 4,
                                           fixed_column = c("zone", "LUI"), 
                                           neg_binom = FALSE)
   
  model_3a <- lmer(log(Simpson_diversity) ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics_cover_filt) # best model
    
  diversity_object[[i]] <- predict_effects(iterations = 1000,
                                           model = model_3a, 
                                           model_data = pollinator_metrics_cover_filt, 
                                           response_variable = "Simpson_diversity",
                                           fixed_number = 2,
                                           factor_number_1 = 2,
                                           factor_number_2 = 4,
                                           fixed_column = c("zone", "LUI"), 
                                           neg_binom = FALSE)
}


abundance_object[[1]]
richness_object[[1]]

abundance_object[[2]]
richness_object[[2]] + scale_y_continuous("")

plot_grid((richness_object[[1]] + xlab("")  + ggtitle("Low forest cover baseline (cover <= 40%)") + theme(legend.position = "none")),
          (richness_object[[2]] + xlab("") + scale_y_continuous("") + ggtitle("High forest cover baseline (cover >= 60%)") + theme(legend.position = "none")),
          abundance_object[[1]] + xlab("")  + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
          abundance_object[[2]] + xlab("") + ylab("NULL") + scale_y_continuous("") + theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = "black"), plot.title = element_text(hjust=0.5)), ncol = 2)
          

ggsave("variable_baseline_temperate_tropical_3.png", scale = 1.1, dpi = 350)
