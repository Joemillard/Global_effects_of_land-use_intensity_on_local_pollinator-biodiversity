# script for including climate in the models as a covariate

### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "lme4", "MASS", "raster", "spdep")
lapply(packages, require, character.only = TRUE)

# source in scripts for deriving predictions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

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

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI", "Sample_end_latest"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE) %>%
  dplyr::filter(!is.na(Latitude))

# create new column for id number
pollinator_metrics <- pollinator_metrics %>%
  mutate("id_col" = 1:nrow(pollinator_metrics))

# PREDICTS sites with the month of the recording
PRED_sites <- pollinator_metrics %>% 
  dplyr::select(Latitude, Longitude, Sample_end_latest, id_col) %>%
  dplyr::mutate(Sample_end_latest = paste(substr(Sample_end_latest, start = 1, stop = 7), sep = "", ".tif")) %>%
  dplyr::filter(!is.na(Latitude))

# read in the climate data for each set of predicts dates
# list the climate specific folders in this directory
temp_dirs <- list.dirs("G:/Extra_data_files/world_clim/temp", recursive = FALSE)
precip_dirs <- list.dirs("G:/Extra_data_files/world_clim/precip", recursive = FALSE)

# set up a list for the directories for temperature and precipitation
climate_dirs <- list(temp_dirs, precip_dirs)

# set up list for values
all_max <- list()

# iterate throught the whole process for temperature and precipitation
for(k in 1:length(climate_dirs)){

  # loop through each directory and create a list of all files
  all.ras <- NULL
  crop.files <- list()
  
  for(i in 1:length(climate_dirs[[k]])){
    crop.files[[i]] <- list.files(climate_dirs[[k]][i])
    crop.files[[i]] <- paste0(climate_dirs[[k]][i], "/", crop.files[[i]])
  }
  
  # file paths for each of per hectare application
  unlisted_crops <- unlist(crop.files)
  
  # for each date of the predicts sites, read in the rasters for those dates and extract the values for those coordinates
  # set up empty list for each dataframe
  ind_raster_frame <- list()
  
  # create unique list of the end dates for the predicts sites
  pred_dates <- unique(PRED_sites$Sample_end_latest)
  
  # set up overall list
  raster_max <- list()
  
  # iterate through each unique date, reading in the raster for that date
  system.time(
    for(i in 1:length(pred_dates)){
      
      # filter the raster for only coordinates we have PREDICTS sites for that date
      PRED_sites_filt <- PRED_sites %>%
        dplyr::filter(Sample_end_latest == pred_dates[i]) %>%
        dplyr::select(Longitude, Latitude) %>%
        SpatialPoints()
      
      # identify site ids for merging
      site_ids <- PRED_sites %>%
        dplyr::filter(Sample_end_latest == pred_dates[i]) %>%
        dplyr::select(Longitude, Latitude, id_col)
      
      # find the index of the end sample date to go 11 months previous
      site_index <- which(grepl(pred_dates[i], unlisted_crops) == TRUE)
      
      # select the previous 11 indices - i.e. 11 months previous worth of pixels
      ind_raster <- unlisted_crops[(site_index - 11): site_index]
  
      # for each of the months previous to the end point, read extract the max temp for the coordinates for that date
      for(j in 1:length(ind_raster)){
        
        # read in the raster from that iteration and extract all the values for that date
        rate_rasters <- raster(ind_raster[[j]])
        
        # filter the raster for that date for the locations we have predicts sites
        PRED_coords <- cbind(site_ids, "max_monthly_temp" = extract(rate_rasters, PRED_sites_filt, na.rm = FALSE))
        
        # convert that set of dates to a dataframe
        ind_raster_frame[[j]] <- as.data.frame(PRED_coords)
        
      }  
        
      raster_max[[i]] <- data.table::rbindlist(ind_raster_frame) %>%
        group_by(id_col) %>%
        mutate("max_temp" = max(max_monthly_temp)) %>%
        ungroup()
      
      # print the iteration number
      print(i) 
        
  })
  
  # remove the max monthly temperature and unique for the annual hottest temperature
  all_max[[k]] <- data.table::rbindlist(raster_max) %>%
    dplyr::select(-max_monthly_temp) %>%
    unique()

}

# checking that the temperatures and preciptation look sensible
all_max[[1]] %>%
  ggplot() +
  geom_point(aes(x = Longitude, y = Latitude, colour = max_temp)) + 
  viridis::scale_colour_viridis() +
  coord_map(projection = "mollweide")

all_max[[2]] %>%
  ggplot() +
  geom_point(aes(x = Longitude, y = Latitude, colour = log10(max_temp))) + 
  viridis::scale_colour_viridis() +
  coord_map(projection = "mollweide")

# join the max temps and precipitation back onto the predicts data
temp_joined_metrics <- inner_join(pollinator_metrics, all_max[[1]], by = "id_col")
temp_joined_metrics <- inner_join(temp_joined_metrics, all_max[[2]], by = "id_col")

# convert to factor and reorder factor levels
temp_joined_metrics$LUI <- factor(temp_joined_metrics$LUI, levels = c("Primary vegetation-Minimal use", 
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
temp_joined_metrics <- droplevels(temp_joined_metrics)

# convert all factor levels to a short acronynm for spacing
temp_joined_metrics <- temp_joined_metrics %>%
  dplyr::mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "PVMU",
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
temp_joined_metrics$Total_abundance <- temp_joined_metrics$Total_abundance + 1
temp_joined_metrics$Simpson_diversity <- temp_joined_metrics$Simpson_diversity + 1

### richness - because only 1 fixed effect, select model random effect structure on basis of AIC
model_1b <- glmer(Species_richness ~ LUI + max_temp.x + log(max_temp.y) + (1|SS) + (1|SSB) + (1|SSBS), data = temp_joined_metrics, family = poisson)
model_1b_2 <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = temp_joined_metrics, family = poisson) 

### abundance
model_2a <- lmer(log(Total_abundance) ~ LUI + max_temp.x + log(max_temp.y) + (1|SS) + (1|SSB), data = temp_joined_metrics)
model_2a_2<- lmer(log(Total_abundance) ~ LUI + (1|SS) + (1|SSB), data = temp_joined_metrics)

# set up each vector
covariate_model <- c(fixef(model_1b)[2:22], fixef(model_2a)[2:22])
LUI_model <- c(fixef(model_1b_2)[2:22], fixef(model_2a_2)[2:22])
metric <- c(rep("Species richness", 21), rep("Total abundance", 21))

# extract the standard errors for all models
covariate_st_err <- c(summary(model_1b)$coefficients[, 2][2:22], summary(model_2a)$coefficients[, 2][2:22])
LUI_st_err <- c(summary(model_1b_2)$coefficients[, 2][2:22], summary(model_2a_2)$coefficients[, 2][2:22])

# extract strings from right side
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# set up the covariate dataframe
covar_data <- data.frame("Covariate_model" = covariate_model,
           "LUI_only_model" = LUI_model,
           "metric" = metric,
           "cov_st_err" = covariate_st_err,
           "LUI_st_err" = LUI_st_err,
           "LUI" = gsub("LUI", "", names(rep(fixef(model_1b)[2:22], 2)))) %>%
  mutate(Intensity = substrRight(LUI, 2)) %>%
  mutate(LUI = gsub("LU", "", LUI)) %>%
  mutate(LUI = gsub("MU", "", LUI)) %>%
  mutate(LUI = gsub("IU", "", LUI)) %>%
  mutate(LUI = factor(LUI, levels = c("PV", "MSV", "ISV", "YSV", "PF", "P", "C", "U"), 
                      labels = c("Primary", "MSV", "ISV", "YSV", "Plantation", "Pasture", "Cropland", "Urban")))

# build the covariate scatter plot
ggplot(covar_data) +
      geom_point(aes(x = Covariate_model, y = LUI_only_model, colour = LUI)) +
      geom_errorbar(aes(x = Covariate_model, ymin = LUI_only_model - LUI_st_err, ymax = LUI_only_model + LUI_st_err, colour = LUI)) +
      geom_errorbar(aes(y = LUI_only_model, xmin = Covariate_model - cov_st_err, xmax = Covariate_model + cov_st_err, colour = LUI)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") + 
      scale_colour_manual("Land-use type",  values = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")) +
      xlab("Covariate LUI fixed effects (biodiversity ~ LUI + max temperature + log(max precipitation))") +
      ylab("LUI only fixed effects (biodiversity ~ LUI)") +
      facet_wrap(~metric) +
      theme_bw() +
      theme(strip.text = element_text(size = 14),
            text = element_text(size = 13),
            panel.grid = element_blank())

ggsave("climate_covariate_model_3.png", scale = 1.2, dpi = 350)

