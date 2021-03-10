### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "lme4", "MASS", "DHARMa", "spdep", "StatisticalModels")
lapply(packages, require, character.only = TRUE)

# extra packages called in place -- plyr and StatisticalModels

source("R/00_functions.R")

# Tim's roquefort function for spatial autocorrelation
SpatialAutocorrelationTest<-function(model,all.data,siteModel=TRUE){
  
  if ("data" %in% names(model)){
    model.data<-model$data
    model.data$res<-residuals(model$model)
  } else {
    model.data<-model@frame
    model.data$res<-residuals(model)
  }
  
  model.data$Longitude<-all.data$Longitude[match(model.data$SSBS,all.data$SSBS)]
  model.data$Latitude<-all.data$Latitude[match(model.data$SSBS,all.data$SSBS)]
  
  studies<-character()
  failed<-character()
  moran.i<-numeric()
  moran.p<-numeric()
  
  i=1
  for (ss in unique(model.data$SS)){
    cat(paste("\rProcessing study ",i," of ",length(unique(model.data$SS)),sep=""))
    data.sub<-droplevels(model.data[model.data$SS==ss,])
    
    if(!siteModel){
      resids<-tapply(data.sub$res,data.sub$SSBS,mean)
      long<-tapply(data.sub$Longitude,data.sub$SSBS,mean)
      lat<-tapply(data.sub$Latitude,data.sub$SSBS,mean)
      data.sub<-data.frame(res=resids,Longitude=long,Latitude=lat)
    }
    
    ds.nb<-try(dnearneigh(cbind(data.sub$Longitude,data.sub$Latitude),
                          d1=0.00000001,d2=10),silent=TRUE)
    ds.listw<-try(nb2listw(ds.nb),silent=TRUE)
    mt<-tryCatch(moran.test(data.sub$res,ds.listw),silent=TRUE,error=function(e) e, 
                 warning=function(w) w)
    
    if(class(mt)[1]=="htest"){
      if ((!is.na(mt$statistic))){
        studies<-c(studies,ss)
        moran.i<-c(moran.i,mt$statistic)
        moran.p<-c(moran.p,mt$p.value)
      } else {
        failed<-c(failed,ss)
      }
      
    } else {
      failed<-c(failed,ss)
    }
    
    i<-i+1
  }
  
  return(list(studies=studies,I=moran.i,P=moran.p,failed=failed))
  
}

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
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI", "UN_region", "Country", "Source_for_predominant_land_use"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# convert to factor and reorder factor levels
pollinator_metrics$LUI <- factor(pollinator_metrics$LUI, levels = c("Primary vegetation-Minimal use", 
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
pollinator_metrics <- droplevels(pollinator_metrics)

# convert all factor levels to a short acronynm for spacing
pollinator_metrics <- pollinator_metrics %>%
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
pollinator_metrics$Total_abundance <- pollinator_metrics$Total_abundance + 1
pollinator_metrics$log_abundance <- log(pollinator_metrics$Total_abundance)

# species richness
model_1b <- GLMER(responseVar = "Species_richness", fixedStruct = "LUI", randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", modelData = pollinator_metrics, fitFamily = "poisson") # best model - failed to converge with 0.00144612

# run the roquefort spatial autocorrelation function to test for autocorrelation at the study level
spatial_output_richness <- SpatialAutocorrelationTest(model_1b, pollinator_metrics)

# total abundance
model_2a <- GLMER(responseVar = "log_abundance", fixedStruct = "LUI", randomStruct = "(1|SS) + (1|SSB)", modelData = pollinator_metrics, fitFamily = "gaussian", saveVars = "SSBS") # best model - failed to converge with 0.00144612

# test for spatial autocorrelation
spatial_output_abundance <- SpatialAutocorrelationTest(model_2a, pollinator_metrics)

# plot histogram of p values for spatial autocorrelation
p_values_richness <- data.frame("p_value" = spatial_output_richness$P, "metric" = "Species richness")
p_values_abundance <- data.frame("p_value" = spatial_output_abundance$P, "metric" = "Total abundance")

# count the number of studies below 0.05
significance_percentage <- rbind(p_values_richness, p_values_abundance) %>%
  mutate(significance = ifelse(p_value < 0.05, "significant", "non-significant")) %>%
  group_by(metric, significance) %>%
  tally() %>%
  ungroup() %>%
  mutate(total = sum(n)) %>%
  mutate(percentage = (n/total) * 100) %>%
  filter(significance == "significant") %>%
  mutate(percentage = paste(as.character(round(percentage, 2)), "% of studies", sep = "")) %>%
  mutate(percentage = paste("p < 0.05 in \n", percentage, sep = " "))
                                         
# bind together the p values and plot
rbind(p_values_richness, p_values_abundance) %>% 
  ggplot() +
    geom_histogram(aes(p_value)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    geom_text(x = 0.17, y = 60, label = "p = 0.05") +
    geom_text(aes(x = 0.25, y = 50, label = percentage), data = significance_percentage) +
    facet_wrap(~metric) +
    ylab("Frequency") +
    xlab("p value") +
    geom_vline(xintercept = 0.05, colour = "red", linetype = "dashed") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 13),
          strip.text = element_text(size = 14))

ggsave("spatial_autocorrelation_figure_2.png", scale = 1, dpi = 350)
