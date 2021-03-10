### packages for analysis,and suppress messages
packages <- c("patchwork", "dplyr", "yarg", "lme4", "raster", "rworldmap", "rworldxtra", "viridis", "ggplot2", "broom.mixed", "gt")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in prediction functions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS(here::here("outputs/PREDICTS_pollinators_8_exp.rds"))

## read in the pesticide data to merge onto the sites
pesticide_high <- raster("C:/Users/joeym/Documents/PhD/Aims/Aim 4 - climate land use interaction effects on pollinators/outputs/2015_Pesticide_totalAPR_High_cropped.tif")

# subset for unqiue sites
sites.sub_xy <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland")) %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  unique() %>%
  SpatialPoints()

# extract the fertislier values and join back onto the coordinates
sites.sub_xy$pest <- extract(pesticide_high, sites.sub_xy, na.rm = FALSE)

# turn the spatial points into a dataframe with the fertiliser data
pest_dat <- data.frame(coords = sites.sub_xy@coords, pest = sites.sub_xy@data$pest)

# build map
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot the fertiliser data at coordinates, with a base map underneath
pest_dat %>%
  dplyr::filter(!is.na(pest)) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "grey") +
  geom_point(aes(x = coords.Longitude, y = coords.Latitude, colour = log10(pest + 0.01))) +
  coord_map(projection = "mollweide") +
  scale_colour_viridis(name = "Total pesticide application rate (kg/ha)", breaks = c(-3, -2, -1, 0, 1, 2), labels = c(expression("1x10"^-3), expression("1x10"^-2), expression("1x10"^-1), expression("1x10"^0), expression("1x10"^1), expression("1x10"^2))) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.text.align = 0)

# save plot of pesticide application rate
ggsave("total_high_pesticide-application-rate_3.png", scale = 1.2, dpi = 350)

### combine the pesticide data with PREDICTS for modelling
# filter cannot decide factors
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland")) %>%
  dplyr::filter(Order %in% c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Apodiformes", "Passeriformes")) %>%
  mutate(confidence_fct = factor(confidence)) %>%
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

# add one too total abundance and diversity for log models
order.sites.div$Total_abundance <- order.sites.div$Total_abundance + 1
order.sites.div$Simpson_diversity <- order.sites.div$Simpson_diversity + 1

# convert the predicts coordinates to spatial points
predicts_points <- order.sites.div %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints()

# use the extract function to extract the pesticide values from both the low and hihg estimate
pesticide_high_points <- extract(pesticide_high, predicts_points)

# build data frame of predicts climate data with the fertiliser data
PREDICTS_fert <- data.frame(order.sites.div, "fert" = pesticide_high_points) %>%
  mutate(fert = fert + 0.01)

# species richness
model_1 <- glmer(Species_richness ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert, family = "poisson")
model_1a <- glmer(Species_richness ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert, family = "poisson")
model_1b <- glmer(Species_richness ~ log10(fert) * Order + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")
model_1b_fert <- glmer(Species_richness ~ log10(fert) + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")
model_1b_int <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")

AIC(model_1b, model_1b_fert)

# plot for species richness
richness_model <- predict_continuous(model = model_1b,
                                     model_data = PREDICTS_fert, 
                                     response_variable = "Species_richness",
                                     categorical_variable = c("Order"),
                                     continuous_variable = c("fert"),
                                     continuous_transformation = log10,
                                     random_variable = c("SS", "SSB", "SSBS"))
# models for total abundance
model_2 <- lmer(log(Total_abundance) ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2a <- lmer(log(Total_abundance) ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert)
model_2_fert <- lmer(log(Total_abundance) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2_int <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_2, model_2_fert)

# plot for total abundance
abundance_model <- predict_continuous(model = model_2,
                                      model_data = PREDICTS_fert, 
                                      response_variable = "Total_abundance",
                                      categorical_variable = c("Order"),
                                      continuous_variable = c("fert"),
                                      continuous_transformation = log10,
                                      random_variable = c("SS", "SSB", "SSBS"))

# models for simpson diversity
model_3 <- lmer(log(Simpson_diversity) ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3a <- lmer(log(Simpson_diversity) ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert)
model_3_fert <- lmer(log(Simpson_diversity) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3_int <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_3, model_3_fert)

# plot for simpson diversity
diversity_model <- predict_continuous(model = model_3,
                                      model_data = PREDICTS_fert, 
                                      response_variable = "Simpson_diversity",
                                      categorical_variable = c("Order"),
                                      continuous_variable = c("fert"),
                                      continuous_transformation = log10,
                                      random_variable = c("SS", "SSB", "SSBS"))

# set up combined dataframe, assign new class variable
combined_metrics <- rbind(richness_model, abundance_model, diversity_model) %>%
  mutate(class = if_else(Order %in% c("Apodiformes", "Passeriformes"), "Vertebrate", "Invertebrate")) %>%   
  group_by(Order) %>%
  filter(fert > quantile(fert, probs = c(0.025))) %>% # for each order, filter for 95% of the data
  filter(fert < quantile(fert, probs = c(0.975))) %>%
  ungroup()

# then split by class and metric
separate_frames <- split(combined_metrics, with(combined_metrics, interaction(class, metric)), drop = TRUE)

# run function for plotting each metric/class combination, and reassign in order consistent with other analysis
plot_list <- lapply(separate_frames, plot_fert_response, categorical_variable = "Order")
plot_list <- list(plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[1]], plot_list[[2]])

# combine all the predicted value plots into one multiplot
{plot_list[[1]] +
    ylab("Species richness") +
    scale_fill_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_colour_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-2.7725887, -2.0794415, -1.3862944, -0.6931472, 0, 0.6931472, 1.3862944, 2.0794415), labels = c("0.06", "0.13", "0.25", "0.5", "1", "2", "4", "8")) +
    xlab(NULL) +
    ggtitle("Invertebrates") +
    labs(subtitle = "A") +
    theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())} +
  
{plot_list[[2]] +
    ylab(NULL) +
    scale_fill_manual(name = "", values = c("#009E73", "black")) +
    scale_colour_manual(name = "", values = c("#009E73", "black")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-0.967584, -0.2876821, 0.4054651, 1.0986123), labels = c("0.38", "0.75", "1.5", "3"))+
    xlab(NULL) +
    ggtitle("Vertebrates") +
    labs(subtitle = "B") +
    theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())} +
  
{plot_list[[3]] +
    ylab("Total abundance") +
    scale_fill_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_colour_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-0.6931472, 0, 0.6931472, 1.3862944, 2.0794415, 2.7725887, 3.465736, 4.158883, 4.85203, 5.545177), labels = c("0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256")) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab(NULL) +
    labs(subtitle = "C")} +
  
{plot_list[[4]] +
    ylab(NULL) +
    scale_fill_manual(name = "", limits = c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Apodiformes", "Passeriformes"), 
                      labels = c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Apodiformes", "Passeriformes"), 
                      values = c("#0072B2", "#CC79A7", "#E69F00","#56B4E9", "#009E73", "black")) +
    scale_colour_manual(name = "", limits = c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Apodiformes", "Passeriformes"), 
                        labels = c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Apodiformes", "Passeriformes"),                                
                        values = c("#0072B2", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "black")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0.6931472, 1.3862944, 2.0794415, 2.7725887), labels = c("2", "4", "8", "16")) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    labs(subtitle = "D") +
    xlab(NULL)} +
  
{plot_list[[5]] +
    ylab("Simpson diversity") +
    scale_fill_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_colour_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.6931472, 1.0986123, 1.3862944, 1.609438), labels = c("1", "2", "3", "4", "5")) +
    xlab("Total pesticide application rate (kg/ha)") +
    theme(legend.position = "none") +
    labs(subtitle = "E")} +
  
{plot_list[[6]] +
    ylab(NULL) +
    scale_fill_manual(name = "", values = c("#009E73", "black")) +
    scale_colour_manual(name = "", values = c("#009E73", "black")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0.6931472, 1.0986123, 1.3862944, 1.609438), labels = c("2", "3", "4", "5")) +
    xlab("Total pesticide application rate (kg/ha)") +
    theme(legend.position = "none") +
    labs(subtitle = "F")} +
  
# set the number of columns for the pesticide plot
plot_layout(ncol = 2)

# save the pesticide plot
ggsave("pesticide_application_rate_biodiversity_1.png", scale = 1.2, dpi = 350)
