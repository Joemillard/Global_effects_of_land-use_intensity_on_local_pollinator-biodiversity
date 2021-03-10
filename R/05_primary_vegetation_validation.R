## calculate forest extent for each predicts primar vegetation site

## find the coordinates for all the primary vegetation sites for temperate and tropical realm
packages <- c("dplyr", "yarg", "cowplot", "here",
              "raster", "ggplot2", "rworldmap", "rworldxtra", "data.table", "viridis")
lapply(packages, require, character.only = TRUE)

# read in extra functions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in the forest data
hansen_tree_cover <- raster(here::here("Data/forest_data/Hansen_full.tif"))

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Primary vegetation")) %>%
  dplyr::filter(Use_intensity == "Minimal use") %>%
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Use_intensity", "Biome"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# assign new variable for tropical/temperate, convert to factor, and filter out NA
pollinator_metrics$zone <- ifelse(pollinator_metrics$Latitude >= -23.5 & pollinator_metrics$Latitude <= 23.5, "Tropical", "Non-tropical")
pollinator_metrics$zone <- factor(pollinator_metrics$zone, levels = c("Non-tropical", "Tropical"))

# select only columns of interest
pollinator_metrics <- pollinator_metrics %>%
  filter(!is.na(zone)) %>%
  dplyr::select(Latitude, Longitude, Predominant_land_use, Use_intensity, zone, Biome) %>%
  mutate(zone = factor(zone))

# build map of minimal use primary vegetation sites
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# pollinator frame for map and rename levels
pollinator_map_data <- pollinator_metrics %>%
  dplyr::select(Latitude, Longitude, zone, Biome) %>%
  unique() %>%
  mutate(Biome = factor(Biome, levels = c("Boreal Forests/Taiga", 
                                          "Temperate Conifer Forests",                               
                                          "Temperate Broadleaf & Mixed Forests", 
                                          "Montane Grasslands & Shrublands",                         
                                          "Temperate Grasslands, Savannas & Shrublands", 
                                          "Mediterranean Forests, Woodlands & Scrub",                
                                          "Deserts & Xeric Shrublands", 
                                          "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                          "Tropical & Subtropical Coniferous Forests", 
                                          "Tropical & Subtropical Dry Broadleaf Forests",            
                                          "Tropical & Subtropical Moist Broadleaf Forests", 
                                          "Mangroves" ),
                        labels = c("Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Grassland & shrubland",
                                   "Grassland & shrubland", 
                                   "Forest", 
                                   "Grassland & shrubland", 
                                   "Grassland & shrubland",
                                   "Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Mangroves")))
  

# build map for distribution of sites
site_distribution <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, fill = Biome), data = pollinator_map_data, shape = 21, colour = "white", alpha = 0.5) +
  scale_fill_manual("Biome", values = c("#228833", "#CCBB44", "#66CCEE")) +
  geom_line(aes(x = x, y = y), linetype = "dashed", colour = "black", 
            data = data.frame(x = c(-175, 175), y = c(23.5, 23.5))) +
  geom_line(aes(x = x, y = y), linetype = "dashed", colour = "black", 
            data = data.frame(x = c(-175, 175), y = c(-23.5, -23.5))) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "A") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 12),
        panel.background = element_rect(fill = "white"))

# map of biomes 
ggsave("biome_map_2.png", scale = 1.1, dpi = 350)

# make plot of biome by temperate/tropical
site_count <- pollinator_metrics %>%
  filter(!is.na(zone)) %>%
  group_by(zone, Biome) %>%
  count(Biome) %>%
  ungroup() %>%
  mutate(Biome = factor(Biome, levels = c("Boreal Forests/Taiga", 
                                          "Temperate Conifer Forests",                               
                                          "Temperate Broadleaf & Mixed Forests", 
                                          "Montane Grasslands & Shrublands",                         
                                          "Temperate Grasslands, Savannas & Shrublands", 
                                          "Mediterranean Forests, Woodlands & Scrub",                
                                          "Deserts & Xeric Shrublands", 
                                          "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                          "Tropical & Subtropical Coniferous Forests", 
                                          "Tropical & Subtropical Dry Broadleaf Forests",            
                                          "Tropical & Subtropical Moist Broadleaf Forests", 
                                          "Mangroves" ),
                        labels = c("Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Grassland & shrubland",
                                   "Grassland & shrubland", 
                                   "Forest", 
                                   "Grassland & shrubland", 
                                   "Grassland & shrubland",
                                   "Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Mangroves"))) %>%
  mutate(Biome = forcats::fct_reorder(Biome, -n)) %>%
  ggplot() +
    geom_bar(aes(x = zone, fill = Biome, y = n), stat = "identity") +
    xlab("") +
    ylab("Site count") +
    labs(subtitle = "C") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
    scale_fill_manual("Biome", values = c("#228833", "#CCBB44", "#66CCEE")) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank())

ggsave("biome_temp_tropica_distribution_2.png", scale = 0.9)

# split metrics up by zone
pollinator_metrics_split <- split(pollinator_metrics, f = pollinator_metrics$zone)

# convert the coordinates for primary site to sp points
convert_spat <- function(data_file) {
  data_fin <- data_file %>%
    dplyr::select(Longitude, Latitude) %>%
    unique() %>%
    SpatialPoints()
  
  return(data_fin)
}

# run function for converting coordinates to spatial
prim_spat <- lapply(pollinator_metrics_split, convert_spat)

# extract the value of forest cover for each primary vegetation site
prim_cover <- list()
for(i in 1:length(prim_spat)){
  prim_cover[[i]] <- extract(hansen_tree_cover, prim_spat[[i]], na.rm = FALSE)
}

# bind the coordinates back onto the extracted coordinates
for(i in 1:length(pollinator_metrics_split)){
  pollinator_metrics_split[[i]] <- pollinator_metrics_split[[i]] %>%
    dplyr::select(Latitude, Longitude, Biome) %>%
    unique() %>%
    cbind(prim_cover[[i]]) %>%
    rename(forest_cover = `prim_cover[[i]]`)
}

# build a map of forest cover values
ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, colour = forest_cover), data = rbindlist(pollinator_metrics_split), alpha = 0.4) +
  scale_colour_viridis("Forest cover (%)") +
  geom_line(aes(x = x, y = y), linetype = "dashed", colour = "black", 
            data = data.frame(x = c(-175, 175), y = c(23.5, 23.5))) +
  geom_line(aes(x = x, y = y), linetype = "dashed", colour = "black", 
            data = data.frame(x = c(-175, 175), y = c(-23.5, -23.5))) +
  coord_map(projection = "mollweide") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 12),
        panel.background = element_rect(fill = "white"))

# save the forest cover plot
ggsave("forest_cover_plot_3.png", scale = 1.1, dpi = 350)

# filter out mangrove as only in tropical
bound_forest_cover <- rbindlist(pollinator_metrics_split) %>%
  filter(Biome != "Mangroves")

# assign new variable for tropical/temperate, convert to factor, and filter out NA
bound_forest_cover$zone <- ifelse(bound_forest_cover$Latitude >= -23.5 & bound_forest_cover$Latitude <= 23.5, "Tropical", "Non-tropical")
bound_forest_cover$zone <- factor(bound_forest_cover$zone, levels = c("Non-tropical", "Tropical"))
   
# reassign the biome factor for forest and grassland
bound_forest_cover %>%
  mutate(Biome = factor(Biome, levels = c("Boreal Forests/Taiga", 
                                          "Temperate Conifer Forests",                               
                                          "Temperate Broadleaf & Mixed Forests", 
                                          "Montane Grasslands & Shrublands",                         
                                          "Temperate Grasslands, Savannas & Shrublands", 
                                          "Mediterranean Forests, Woodlands & Scrub",                
                                          "Deserts & Xeric Shrublands", 
                                          "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                          "Tropical & Subtropical Coniferous Forests", 
                                          "Tropical & Subtropical Dry Broadleaf Forests",            
                                          "Tropical & Subtropical Moist Broadleaf Forests"),
                        labels = c("Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Grassland & shrubland",
                                   "Grassland & shrubland", 
                                   "Forest", 
                                   "Grassland & shrubland", 
                                   "Grassland & shrubland",
                                   "Forest", 
                                   "Forest", 
                                   "Forest"))) %>%
  group_by(zone, Biome) %>%
  summarise(mean_val = mean(forest_cover, na.rm = TRUE)) %>%
  ungroup()

# make a boxplot of forest cover
cover_boxplot <- bound_forest_cover %>%
  mutate(Biome = factor(Biome, levels = c("Boreal Forests/Taiga", 
                                          "Temperate Conifer Forests",                               
                                          "Temperate Broadleaf & Mixed Forests", 
                                          "Montane Grasslands & Shrublands",                         
                                          "Temperate Grasslands, Savannas & Shrublands", 
                                          "Mediterranean Forests, Woodlands & Scrub",                
                                          "Deserts & Xeric Shrublands", 
                                          "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                          "Tropical & Subtropical Coniferous Forests", 
                                          "Tropical & Subtropical Dry Broadleaf Forests",            
                                          "Tropical & Subtropical Moist Broadleaf Forests"),
                        labels = c("Forest", 
                                   "Forest", 
                                   "Forest", 
                                   "Grassland & shrubland",
                                   "Grassland & shrubland", 
                                   "Forest", 
                                   "Grassland & shrubland", 
                                   "Grassland & shrubland",
                                   "Forest", 
                                   "Forest", 
                                   "Forest"))) %>% group_by(zone, Biome) %>% tally()
  ggplot() +
    geom_boxplot(aes(x = zone, y = forest_cover, fill = Biome)) +
    geom_segment(y = -4, yend = -4, x = 0.4, xend = 3) +
    geom_segment(y = -4, yend = 110, x = 0.4, xend = 0.4) +
    coord_cartesian(ylim = c(-4, 105)) +
    geom_text(x = 0.8, y = 106, label = "n = 532", size = 3) +
    geom_text(x = 1.2, y = 106, label = "n = 162", size = 3) +
    geom_text(x = 1.8, y = 106, label = "n = 339", size = 3) +
    geom_text(x = 2.2, y = 106, label = "n = 88", size = 3) +
    xlab("") +
    ylab("Forest cover (%)") +
    labs(subtitle = "B") +
    scale_fill_manual("Biome", values = c("#228833", "#CCBB44")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none", panel.border= element_blank())

# save the boxplot
ggsave("forest_cover_boxplot_2.png", scale = 0.9, dpi = 350)

# build combined plot for primary vegetation minimal use
plot_grid(site_distribution, plot_grid(cover_boxplot, NULL, site_count, rel_widths = c(1, 0.05, 1), ncol = 3), ncol = 1)

# save the combined plot for primary vegetation - minimal use
ggsave("primary_veg_minimal_use_2.png", dpi = 350, scale = 1)
