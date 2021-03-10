# script for correlating pesticide against fertiliser application
# packages for analysis,and suppress messages
packages <- c( "dplyr", "ggplot2", "raster")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in prediction functions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS(here::here("outputs/PREDICTS_pollinators_5.rds"))

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

## read in the fertiliser data
fert_data <- raster(here::here("fertiliser_application_rate.tif"))

# extract the fertislier values and join back onto the coordinates
sites.sub_xy$fert <- extract(fert_data, sites.sub_xy, na.rm = FALSE)

# turn the spatial points into a dataframe with the fertiliser data
fert_dat <- data.frame(coords = sites.sub_xy@coords, fert = sites.sub_xy@data$fert) %>%
  dplyr::select(fert)

# bind the pesticide and fertiliser data
bound_pest_fert <- cbind(fert_dat, pest_dat) %>%
  dplyr::select(coords.Longitude, coords.Latitude, fert, pest)

# build a plot for correlation and fertiliser and pesticide application
ggplot(bound_pest_fert) +
  geom_point(aes(x = fert, y = pest), shape = 21) +
  xlab("Total fertiliser application rate (kg/ha)") +
  ylab("Total pesticide application rate (kg/ha)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4500)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("fert-pest_correlation.png", scale = 1, dpi = 350)
