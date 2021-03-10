library(raster)
library(here)
library(ggplot2)
library(rasterVis)

### new script for per hectare (rate application) aggregation to check gives same result as
# Charlie's calculation from total application

# list the crop specific folders in this directory
cropdirs <- list.dirs(here("Data/Charlie_intensity_analogues/FertilizerCropSpecific_Geotiff"), recursive = FALSE)

# loop through each directory and create a list of all files
all.ras <- NULL
crop.files <- list()

for(i in 1:length(cropdirs)){
  crop.files[[i]] <- list.files(cropdirs[i], pattern = "Rate.tif$")
  crop.files[[i]] <- paste0(cropdirs[i], "/", crop.files[[i]])
}

# file paths for each of per hectare application
unlisted_crops <- unlist(crop.files)
rate_rasters <- list()

# read in each of the rasters
for(i in 1:length(unlisted_crops)){
  rate_rasters[[i]] <- raster(unlisted_crops[i])
  print(i)
}

# organise all of the rasters into a stack and sum - change so don't remove 0 values
stacked_rasters <- stack(rate_rasters)
fert.total <- sum(stacked_rasters, na.rm = T)
#values(fert.total)[values(fert.total) == 0] <- NA

writeRaster(fert.total, NAflag = -3.4e+38, "fertiliser_application_rate_1.tif", overwrite = TRUE)

# read in new ate raster
fertiliser_raster <- raster("fertiliser_application_rate_1.tif")

# plot the rate raster
gplot(fertiliser_raster) +
  geom_raster(aes(fill = value), alpha = 0.7, na.rm = T) +
  scale_fill_distiller(direction = 1, name = "Tons") +
  ggtitle("Application rate, 17 crops") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(text = element_text(size = 10)) +
  theme_bw()





