### packages for analysis,and suppress messages
packages <- c("reshape2", "gt", "patchwork", "dplyr", "yarg", "lme4", "raster", "rworldmap", "rworldxtra", "viridis", "ggplot2", "broom.mixed", "gt")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in prediction functions
source("R/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS(here::here("outputs/PREDICTS_pollinators_8_exp.rds"))

# read in the fertiliser data
fert_data <- raster(here::here("ouputs/fertiliser_application_rate_1.tif"))

# subset for unqiue sites
sites.sub_xy <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland")) %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  unique() %>%
  SpatialPoints()

# extract the fertiliser values and join back onto the coordinates
sites.sub_xy$fert <- extract(fert_data, sites.sub_xy, na.rm = FALSE)

# turn the spatial points into a dataframe with the fertiliser data
fert_dat <- data.frame(coords = sites.sub_xy@coords, fert = sites.sub_xy@data$fert)

# build map
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot the fertiliser data at coordinates, with a base map underneath
fert_dat %>%
  dplyr::filter(!is.na(fert)) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "grey") +
  geom_point(aes(x = coords.Longitude, y = coords.Latitude, colour = log10(fert + 1))) +
  scale_colour_viridis(name = "Total fertiliser application rate (kg/ha)", option = "magma", breaks = c(1, 2, 3), labels = c(expression("1x10"^1), expression("1x10"^2), expression("1x10"^3))) +
  coord_map(projection = "mollweide") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

# save the intensity plot
ggsave("site_fertiliser_rate_2.png", scale = 1.2, dpi = 350)

# PREDICTS data compilation

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

# assign new variable for tropical/temperate, convert to factor, and filter out NA - not currently used
order.sites.div$zone <- ifelse(order.sites.div$Latitude >= -23.5 & order.sites.div$Latitude <= 23.5, "Tropics", "Temperate")
order.sites.div$zone <- factor(order.sites.div$zone, levels = c("Temperate", "Tropics"))

# add one too total abundance and diversity for log models
order.sites.div$Total_abundance <- order.sites.div$Total_abundance + 1
order.sites.div$Simpson_diversity <- order.sites.div$Simpson_diversity + 1

# merge the fert data with the PREDICTS sites
PREDICTS_fert <- inner_join(order.sites.div, fert_dat, by = c("Longitude" = "coords.Longitude", "Latitude" = "coords.Latitude")) %>%
  mutate(fert  = fert + 1)

# table for number of sites in each order
figure_5_table <- PREDICTS_fert %>%
  group_by(Order) %>%
  tally() %>%
  arrange(desc(n)) %>%
  gt() %>%
  tab_header(
    title = "Global cropland site representation (Figure 5)") %>%
  fmt_number(
    columns = vars(n),
    drop_trailing_zeros = TRUE,
    sep_mark = ""
  )

# species richness
model_1 <- glmer(Species_richness ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert, family = "poisson")
model_1a <- glmer(Species_richness ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert, family = "poisson")
model_1b <- glmer(Species_richness ~ log10(fert) * Order + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")

AIC(model_1, model_1a, model_1b)

model_1b_fert <- glmer(Species_richness ~ log10(fert) + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")
model_1b_int <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson")

AIC(model_1b, model_1b_fert)

# build table of best fitting model
model_1b_table <- tidy(model_1b) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:12) %>%
  gt() %>%
  tab_header(
    title = "Model summary - species richness (Figure 5)")

# derive predictions for species richness
richness_model <- predict_continuous(model = model_1b,
                              model_data = PREDICTS_fert, 
                              response_variable = "Species_richness",
                              categorical_variable = c("Order"),
                              continuous_variable = c("fert"),
                              continuous_transformation = log10,
                              random_variable = c("SS", "SSB", "SSBS"))

### 
model_2 <- lmerTest::lmer(log(Total_abundance) ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2a <- lmer(log(Total_abundance) ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert)

summary(model_2)

AIC(model_2, model_2a)

model_2_fert <- lmer(log(Total_abundance) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_2_int <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_2, model_2_fert)

# build table of best fitting model
model_2_table <- tidy(model_2) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:12) %>%
  gt() %>%
  tab_header(
    title = "Model summary - total abundance (Figure 5)")

# derive predictions for total abundance
abundance_model <- predict_continuous(model = model_2,
                                     model_data = PREDICTS_fert, 
                                     response_variable = "Total_abundance",
                                     categorical_variable = c("Order"),
                                     continuous_variable = c("fert"),
                                     continuous_transformation = log10,
                                     random_variable = c("SS", "SSB", "SSBS"))

### simpson diversity
model_3 <- lmerTest::lmer(log(Simpson_diversity) ~ log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3a <- lmer(log(Simpson_diversity) ~ log10(fert) * Order + (1|SS), data = PREDICTS_fert)

summary(model_3)

AIC(model_3, model_3a)

model_3_fert <- lmer(log(Simpson_diversity) ~ log10(fert) + (1|SS) + (1|SSB), data = PREDICTS_fert)
model_3_int <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = PREDICTS_fert)

AIC(model_3, model_3_fert)

# build table of best fitting model
model_3_table <- tidy(model_3) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:12) %>%
  gt() %>%
  tab_header(
    title = "Model summary - Simpson diversity (Figure 5)")

# derive predictions for simpson diversity
diversity_model <- predict_continuous(model = model_3,
                                     model_data = PREDICTS_fert, 
                                     response_variable = "Simpson_diversity",
                                     categorical_variable = c("Order"),
                                     continuous_variable = c("fert"),
                                     continuous_transformation = log10,
                                     random_variable = c("SS", "SSB", "SSBS"))

# function for extracting interaction p value - fertiliser in intercept model
extract_p_value <- function(model){
  summary(model)$coef[1,5]
}

extract_p_value_glmer <- function(model){
  summary(model)$coef[1,4]
}

# set up empty vectors for p values
p_values_model_1b <- c()
p_values_model_2 <- c()
p_values_model_3 <- c()
ref_factors <- c()

# extract interaction significance
for(i in 1:length(levels(PREDICTS_fert$Order))){
  
  # set the reference factor
  PREDICTS_fert$Order <- relevel(PREDICTS_fert$Order, ref = levels(PREDICTS_fert$Order)[i])
  
  # retrieve reference factor for that iteration
  ref_factors[i] <- levels(PREDICTS_fert$Order)[1]
  
  # print the list of factors for the reference order
  print(levels(PREDICTS_fert$Order))
  
  # for each model with different reference factor, rerun that model and extract the p values
  p_values_model_1b[i] <- glmer(Species_richness ~ 0 + log10(fert) * Order + (1|SS) + (1|SSB) + (1|SSBS), data = PREDICTS_fert, family = "poisson") %>%
    extract_p_value_glmer()
  p_values_model_2[i] <- lmerTest::lmer(log(Total_abundance) ~ 0 + log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert) %>%
    extract_p_value()
  p_values_model_3[i] <- lmerTest::lmer(log(Simpson_diversity) ~ 0 + log10(fert) * Order + (1|SS) + (1|SSB), data = PREDICTS_fert) %>%
    extract_p_value()
}

# create dataframe of all the p value interactions
all_p_value <- data.frame("reference" = ref_factors, 
                          "Species_richness" = p_values_model_1b, 
                          "Total_abundance" = p_values_model_2, 
                          "Simpson_diversity" = p_values_model_3)

# melt the dataframe of p values
p_value_table <- melt(all_p_value) %>%
  mutate(significant = ifelse(value < 0.05, "significant", "non-significant"))  %>%
  mutate(reference = factor(reference)) %>%
  mutate(value = factor(variable))

# table of AICs
selection_table <- data.frame("Response" = c(rep("Species richness", 3),
                                             rep("Total abundance", 3),
                                             rep("Simpson diversity", 3)),
                              "Model" = c("Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)", 
                                          "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                          "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                          "Simpson_diversity ~ Order * LUI + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ LUI + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ 1 + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(model_1b), AIC(model_1b_fert), AIC(model_1b_int), 
                                        AIC(model_2), AIC(model_2_fert), AIC(model_2_int),
                                        AIC(model_3), AIC(model_3_fert), AIC(model_3_int))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# set up combined dataframe, assign new class variable
combined_metrics <- rbind(richness_model, abundance_model, diversity_model) %>%
  mutate(class = if_else(Order %in% c("Apodiformes", "Passeriformes"), "Vertebrate", "Invertebrate")) %>%   
  group_by(Order) %>%
  filter(fert > quantile(fert, probs = c(0.025))) %>% # for each order, filter for 95% of the data
  filter(fert < quantile(fert, probs = c(0.975))) %>%
  ungroup()

# merge the predicted values with the significance
combined_metrics <- left_join(combined_metrics, p_value_table, by = c("Order" = "reference", "metric" = "variable"))

# then split by class and metric
separate_frames <- split(combined_metrics, with(combined_metrics, interaction(class, metric)), drop = TRUE)

# calculate percentage change per 1000 kg/ha for each taxonomic order
diff_fert <- function(data_file){
  data_fin <- data_file %>%
    group_by(Order) %>%
    mutate(fert_diff = max(fert) - min(fert)) %>%
    mutate(y_value = exp(y_value)) %>%
    mutate(resp_change = ((max(y_value) - min(y_value)) / min(y_value)) * 100) %>%
    mutate(per_1000 = (resp_change/fert_diff) * 1000) %>%
    ungroup() %>%
    dplyr::select(Order, fert_diff, resp_change, per_1000) %>%
    unique()
  
  return(data_fin)
  
}

# run function to calculate change per 1000 kg/ha
change_per_1000 <- lapply(separate_frames, diff_fert)

# run function for plotting each metric/class combination, and reassign in order consistent with other analysis
plot_list <- lapply(separate_frames, plot_fert_response, categorical_variable = "Order")
plot_list <- list(plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[1]], plot_list[[2]])

# extract all the values for the data for a table of all values
model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[1]] %>%
    dplyr::select(colour, x, y) %>% unique() %>%
    cbind(unique(ggplot_build(model_plot)$data[[2]]) %>% dplyr::select(ymin, ymax))
}

# set up colours for invertebrates
invertebrate_colours <- function(data_1){
  data_1$Order[data_1$colour == "#F8766D"] <- "Coleoptera"
  data_1$Order[data_1$colour == "#7CAE00"] <- "Diptera"
  data_1$Order[data_1$colour == "#00BFC4"] <- "Hymenoptera"
  data_1$Order[data_1$colour == "#C77CFF"] <- "Lepidoptera"
  return(data_1)
}

# set up colours for vertebrates
vertebrate_colours <- function(data_1){
  data_1$Order[data_1$colour == "#F8766D"] <- "Apodiformes"
  data_1$Order[data_1$colour == "#00BFC4"] <- "Passeriformes"
  return(data_1)
}

# bind together all the prediciton data
final_prediction_data <- rbind(model_data(plot_list[[1]]) %>% mutate(metric = "Richness") %>% invertebrate_colours(),
  model_data(plot_list[[2]]) %>% mutate(metric = "Richness") %>% vertebrate_colours(),
  model_data(plot_list[[3]]) %>% mutate(metric = "Abundance") %>% invertebrate_colours(),
  model_data(plot_list[[4]]) %>% mutate(metric = "Abundance") %>% vertebrate_colours(),
  model_data(plot_list[[5]]) %>% mutate(metric = "Diversity") %>% invertebrate_colours(),
  model_data(plot_list[[6]]) %>% mutate(metric = "Diversity") %>% vertebrate_colours()) %>%
  dplyr::select(-colour)

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
    #scale_linetype_discrete() +
    xlab(NULL) +
    ggtitle("Vertebrates") +
    labs(subtitle = "B") +
    guides(fill = FALSE, colour = FALSE) +
    theme(plot.title = element_text(hjust=0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())} +
  
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
    guides(linetype = FALSE) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    labs(subtitle = "D") +
    xlab(NULL)} +
  
{plot_list[[5]] +
    ylab("Simpson diversity") +
    scale_fill_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_colour_manual(name = "", values = c("#56B4E9", "#E69F00", "#0072B2", "#CC79A7")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.6931472, 1.0986123, 1.3862944, 1.609438), labels = c("1", "2", "3", "4", "5")) +
    xlab("Total fertiliser application rate (kg/ha)") +
    theme(legend.position = "none") +
    labs(subtitle = "E")} +
  
{plot_list[[6]] +
    ylab(NULL) +
    scale_fill_manual(name = "", values = c("#009E73", "black")) +
    scale_colour_manual(name = "", values = c("#009E73", "black")) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0.6931472, 1.0986123, 1.3862944, 1.609438), labels = c("2", "3", "4", "5")) +
    xlab("Total fertiliser application rate (kg/ha)") +
    theme(legend.position = "none") +
    labs(subtitle = "F")} +
  
  plot_layout(ncol = 2)

ggsave("fertiliser_rate_7.pdf", scale = 1.1, dpi = 350)
