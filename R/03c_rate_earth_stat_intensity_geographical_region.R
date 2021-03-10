### packages for analysis,and suppress messages
packages <- c("patchwork", "dplyr", "yarg", "cowplot",
              "lme4", "raster", "rworldmap", "rworldxtra", "viridis", "ggplot2", "broom.mixed", "gt")

lapply(packages, require, character.only = TRUE)

# source in prediction functions
source("R/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS(here::here("outputs/PREDICTS_pollinators_8_exp.rds"))

# read in the fertiliser data
fert_data <- raster(here::here("outputs/fertiliser_application_rate.tif"))

# subset for unqiue sites
sites.sub_xy <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland")) %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  unique() %>%
  SpatialPoints()

# extract the fertislier values and join back onto the coordinates
sites.sub_xy$fert <- extract(fert_data, sites.sub_xy, na.rm = FALSE)

# turn the spatial points into a dataframe with the fertiliser data
fert_dat <- data.frame(coords = sites.sub_xy@coords, fert = sites.sub_xy@data$fert)

# PREDICTS data compilation

# filter cannot decide factors
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland")) %>%
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# create object as list of pollinator subsets, and then use apply to calculate metrics
# split diversity data into list of four for each order
diversityOrder <- split(x = PREDICTS_pollinators, f = PREDICTS_pollinators$Order)

# drop unused levels from each list
diversityOrder <- lapply(diversityOrder, function(x) return(droplevels(x)))

# Remove empty rows
diversityOrder <- Filter(function(x) dim(x)[1] > 0, diversityOrder)

# calculate site metrics for each of the four order level subsets
order.sites.div <- do.call('rbind',lapply(X = diversityOrder, FUN = SiteMetrics,
                                          extra.cols = c("SSB", "SSBS","Biome", "Sampling_method",
                                                         "Study_common_taxon", "Sampling_effort",
                                                         "Sampling_effort_unit", "Realm",
                                                         "Predominant_land_use", "Order", "Order_use", "confidence_fct"),
                                          sites.are.unique = TRUE, srEstimators = TRUE))

# assign new variable for tropical/temperate, convert to factor, and filter out NA - not currently used
order.sites.div$zone <- ifelse(order.sites.div$Latitude >= -23.5 & order.sites.div$Latitude <= 23.5, "Tropics", "Temperate")
order.sites.div$zone <- factor(order.sites.div$zone, levels = c("Temperate", "Tropics"))

# add one too total abundance and diversity for log models
order.sites.div$Total_abundance <- order.sites.div$Total_abundance + 1
order.sites.div$Simpson_diversity <- order.sites.div$Simpson_diversity + 1

# merge the fert data with the PREDICTS sites
PREDICTS_fert <- inner_join(order.sites.div, fert_dat, by = c("Longitude" = "coords.Longitude", "Latitude" = "coords.Latitude"))

# build model predicting species richness against fertiliser and tropical/temperate zone
model_1 <- glmer(Species_richness ~ log10(fert) * zone + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")
model_1a <- glmer(Species_richness ~ log10(fert) * zone + (1|SS) + (1|SSB), data = PREDICTS_fert, family = "poisson")
model_1b <- glmer(Species_richness ~ log10(fert) * zone + (1|SS), data = PREDICTS_fert, family = "poisson")

AIC(model_1, model_1a, model_1b)

model_1_fert <- glmer(Species_richness ~ log10(fert) + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")
model_1_int <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")

AIC(model_1, model_1_fert)

# build table of best fitting model
model_1_table <- tidy(model_1) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  gt() %>%
  tab_header(
    title = "Figure 4 model summary - species richness")

# predict for tropics/temperate for species richness
richness_model <- predict_continuous(model = model_1,
                                     model_data = PREDICTS_fert, 
                                     response_variable = "Species_richness",
                                     categorical_variable = c("zone"),
                                     continuous_variable = c("fert"),
                                     continuous_transformation = log10,
                                     random_variable = c("SS", "SSB", "SSBS"))

### 
model_2 <- lmer(log(Total_abundance) ~ log10(fert) * zone + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2a <- lmer(log(Total_abundance) ~ log10(fert) * zone + (1|SS), data = PREDICTS_fert)

AIC(model_2, model_2a)

model_2_fert <- lmer(log(Total_abundance) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2_int <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_2, model_2_fert)

# build table of best fitting model
model_2_table <- tidy(model_2) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  gt() %>%
  tab_header(
    title = "Figure 4 model summary - total abundance")

# derive predictions for abundance
abundance_model <- predict_continuous(model = model_2,
                                      model_data = PREDICTS_fert, 
                                      response_variable = "Total_abundance",
                                      categorical_variable = c("zone"),
                                      continuous_variable = c("fert"),
                                      continuous_transformation = log10,
                                      random_variable = c("SS", "SSB", "SSBS"))

### simpson diversity
model_3 <- lmer(log(Simpson_diversity) ~ log10(fert) * zone + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3a <- lmer(log(Simpson_diversity) ~ log10(fert) * zone + (1|SS), data = PREDICTS_fert)

AIC(model_3, model_3a)

model_3_fert <- lmer(log(Simpson_diversity) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3_int <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_3, model_3_fert)

# build table of best fitting model
model_3_table <- tidy(model_3) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  gt() %>%
  tab_header(
    title = "Figure 4 model summary - Simpson diversity")

# derive predictions for simpson diversity
diversity_model <- predict_continuous(model = model_3,
                                      model_data = PREDICTS_fert, 
                                      response_variable = "Simpson_diversity",
                                      categorical_variable = c("zone"),
                                      continuous_variable = c("fert"),
                                      continuous_transformation = log10,
                                      random_variable = c("SS", "SSB", "SSBS"))

# table of AICs
selection_table <- data.frame("Response" = c(rep("Species richness", 3),
                                             rep("Total abundance", 3),
                                             rep("Simpson diversity", 3)),
                              "Model" = c("Species_richness ~ fert * zone + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ fert + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Total_abundance ~ fert * zone + (1|SS) + (1|SSB)", 
                                          "Total_abundance ~ fert + (1|SS) + (1|SSB)",
                                          "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                          "Simpson_diversity ~ fert * zone + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ fert + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ 1 + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(model_1), AIC(model_1_fert), AIC(model_1_int), 
                                        AIC(model_2), AIC(model_2_fert), AIC(model_2_int),
                                        AIC(model_3), AIC(model_3_fert), AIC(model_3_int))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# set up combined dataframe, assign new class variable
combined_metrics <- rbind(richness_model, abundance_model, diversity_model)  %>%   
  group_by(zone) %>%
  filter(fert > quantile(fert, probs = c(0.025))) %>% # for each order, filter for 95% of the data
  filter(fert < quantile(fert, probs = c(0.975))) %>%
  ungroup()

#  and then split by class and metric
separate_frames <- split(combined_metrics, with(combined_metrics, interaction(metric)), drop = TRUE)

# run function for plotting each metric/class combination, and reassign in order consistent with other analysis
plot_list <- lapply(separate_frames, plot_fert_response, categorical_variable = "zone")
plot_list <- list(plot_list[[2]], plot_list[[3]], plot_list[[1]])

# add new scales for amended colours (blue and orange)
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + scale_fill_manual(values = c("#5D69B1", "#E58606"), labels = c("Non-tropical", "Tropical")) + scale_colour_manual(values = c("#5D69B1", "#E58606"), labels = c("Non-tropical", "Tropical"))
}

# retrieve the legend to put at bottom right panel
add_legend <- get_legend(
  
  # create some space to the left of the legend
  plot_list[[1]] +
    theme(legend.background = element_rect(colour = "white"),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.position = "bottom") +
    labs(colour = "Geographical zone") +
    labs(fill = "Geographical zone"))

# build the multiplot with legend at bottom right panel
plot_set_up <- 
  plot_grid((plot_list[[1]] + 
               ylab("Species richness") + 
               xlab("Total fertiliser application rate (kg/ha)") + 
               theme(legend.position = "none") + 
               scale_y_continuous(breaks = c(0, 0.6931472, 1.3862944, 2.0794415), labels = c("1", "2", "4", "8"))),
  (plot_list[[2]] + 
     ylab("Total abundance") + 
     xlab("Total fertiliser application rate (kg/ha)") + 
     theme(legend.position = "NULL") +
     scale_y_continuous(breaks = c(0, 0.6931472, 1.3862944, 2.0794415, 2.7725887, 3.465736, 4.158883), labels = c("1", "2", "4", "8", "16", "32", "64"))),
  (plot_list[[3]] + 
     ylab("Simpson diversity") + 
     xlab("Total fertiliser application rate (kg/ha)") + 
     theme(legend.position = "NULL") +
     scale_y_continuous(limits = c(0.65, 1.63), breaks = c(0.6931472, 1.0986123, 1.3862944, 1.609438), labels = c("2", "3", "4", "5"))),
  add_legend,
  labels = c("A", "B", "C"),
  label_size = 13,
  ncol = 2,
  rel_widths = c(1, 1))

ggsave("fertiliser_response_zone_5.png", scale = 1.1, dpi = 400)
          