### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "lme4", "patchwork", "gt", "broom.mixed", "MASS")
lapply(packages, require, character.only = TRUE)

# extra packages called in place -- plyr and StatisticalModels

source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# check use_intensity and use_type factors
table(PREDICTS_pollinators$Use_intensity)
table(PREDICTS_pollinators$Predominant_land_use)
table(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity)

# check the distribution of factors for intensity and type
table(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity)

# bind together the intensity and type into a single variable
PREDICTS_pollinators$LUI <- paste(PREDICTS_pollinators$Predominant_land_use, PREDICTS_pollinators$Use_intensity, sep = "-")

PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(LUI %in% c("Primary vegetation-Minimal use","Cropland-Minimal use", "Cropland-Light use", "Cropland-Intense use")) %>%	 
  droplevels()

# check the taxa in PREDICTS_pollinators
PREDICTS_pollinators %>% group_by(Order) %>% tally()

# count the number of sites - Mature secondary intense (5), Intermediate secondary intense (56), Urban intense (44)
site_count <- PREDICTS_pollinators %>% dplyr::select(LUI, SSBS) %>% unique()
table(site_count$LUI)

# re-check the number for each factor for LUI
table(PREDICTS_pollinators$LUI)

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
pollinator_metrics <- SiteMetrics(diversity = PREDICTS_pollinators,
                                  extra.cols = c("SSB", "SSBS", "LUI"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# convert to factor and reorder factor levels
pollinator_metrics$LUI <- factor(pollinator_metrics$LUI, levels = c("Primary vegetation-Minimal use",
                                                                    "Cropland-Minimal use",
                                                                    "Cropland-Light use",
                                                                    "Cropland-Intense use"))
pollinator_metrics <- droplevels(pollinator_metrics)

# convert all factor levels to a short acronynm for spacing
pollinator_metrics <- pollinator_metrics %>%
  mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "Primary vegetation",
                              "Cropland-Minimal use" = "Cropland (Minimal use)",
                              "Cropland-Light use" = "Cropland (Light use)",
                              "Cropland-Intense use" = "Cropland (Intense use)")))

# assign new variable for tropical/temperate, convert to factor, and filter out NA
pollinator_metrics$zone <- ifelse(pollinator_metrics$Latitude >= -23.5 & pollinator_metrics$Latitude <= 23.5, "Tropics", "Temperate")
pollinator_metrics$zone <- factor(pollinator_metrics$zone, levels = c("Temperate", "Tropics"))
pollinator_metrics <- pollinator_metrics %>%
  filter(!is.na(zone))

# multi panel plot of abundance, richness, and diversity
# add 1 for abundance and simpson diversity
pollinator_metrics$Total_abundance <- pollinator_metrics$Total_abundance + 1
pollinator_metrics$Simpson_diversity <- pollinator_metrics$Simpson_diversity + 1

## function for randomly sampling sets of 1000 sites from tropical and non-tropical
# first count the number of sites in tropical and temperate 
zone_count <- pollinator_metrics %>%
  mutate(zone = as.character(zone)) %>%
  group_by(zone) %>%
  tally() %>%
  ungroup()

# function for filtering the sampled data
sample_rows <- function(data_file, zone){
  
  # filter for number of rows for that zone
  zone_rows <- zone_count$n[zone_count$zone == zone]
  
  # sample from the rows 1000 times
  row_samp <- sample(1:zone_rows, size = 1000)
  
  # filter the datafile for that number of rows from the zone data
  data_fin <- data_file[data_file$zone == zone,]
  data_fin <- data_fin[row_samp,]
  
  return(data_fin)
}

# take 10 samples from both the temperate and tropical subsets
fin_samp <- list()
for(i in 1:100){
  tropic_samp <- sample_rows(pollinator_metrics, zone = "Tropics")
  temp_samp <- sample_rows(pollinator_metrics, zone = "Temperate")
  
  fin_samp[[i]] <- rbind(tropic_samp, temp_samp)
  
}

# for each set of samples run it over the total abundance model
abundance_metric_samp <- list()
abundance_metric_samp_dat <- list()
model_2a_samp <- list()
for(i in 1:length(fin_samp)){
  model_2a_samp[[i]] <- lmer(log(Total_abundance) ~ LUI * zone + (1|SS) + (1|SSB), data = fin_samp[[i]]) # best model
  
  # predict values from the model                                
  abundance_metric_samp[[i]] <- predict_effects(iterations = 1000,
                                      model = model_2a_samp[[i]], 
                                      model_data = fin_samp[[i]], 
                                      response_variable = "Total_abundance",
                                      fixed_number = 2,
                                      factor_number_1 = 2,
                                      factor_number_2 = 4,
                                      fixed_column = c("zone", "LUI"),
                                      neg_binom = FALSE)
  
  # extract the data from each abundance model
  abundance_metric_samp_dat[[i]] <- ggplot_build(abundance_metric_samp[[i]])$data[[3]] %>%
    dplyr::select(ymin, ymax) %>%
    mutate(difference = abs(ymin - ymax)) %>%
    mutate(zone = c(rep("Temperate", 4), rep("Tropical", 4))) %>%
    filter(!is.na(ymin)) %>%
    mutate(intensity = c(rep(c("Minimal", "Light", "Intense"), 2))) %>%
    mutate(sample_no = i)
}

# plot the distribution of the confidence interval sizes as a boxplot
data.table::rbindlist(abundance_metric_samp_dat) %>%
  mutate(intensity = factor(intensity, levels = c("Minimal", "Light", "Intense"))) %>%
  mutate(zone = factor(zone, levels = c("Temperate", "Tropical"), labels = c("Non-tropical (sites = 1000)", "Tropical (sites = 1000)"))) %>%
  ggplot() +
  geom_violin(aes(x = "", y = difference, fill = intensity), draw_quantiles  = c(0.5)) +
  ylab("Resampled 95% CI range") + 
  xlab("") +
  scale_fill_manual("Cropland use-intensity", values = c("yellow2", "orange", "red")) +
  facet_wrap(~zone) +
  scale_y_continuous(limits = c(0, 400), breaks = c(0, 100, 200, 300, 400), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        text = element_text(size = 13),
        strip.text = element_text(size = 13))

# save the resampled plot
ggsave("sample_variation_2.png", scale = 1.1, dpi = 350)

# build table for site representation
figure_3_table <- pollinator_metrics %>%
  mutate(zone = factor(zone, levels = c("Temperate", "Tropics"), labels = c("Non-tropical", "Tropical"))) %>%
  group_by(zone, LUI) %>%
  tally() %>%
  gt() %>%
  tab_header(
    title = "Non-tropical/tropical cropland site representation (Figure 3)") %>%
  fmt_number(
    columns = vars(n),
    drop_trailing_zeros = TRUE,
    sep_mark = ""
  )

# save the table as png
gtsave(figure_3_table, "figure_3_table_exp.png")

### richness - because only 1 fixed effect, select model random effect structure on basis of AIC
model_1 <- glmer(Species_richness ~ LUI * zone + (1|SS), data = pollinator_metrics, family = poisson)
summary(model_1)
blmeco::dispersion_glmer(model_1) 
StatisticalModels::GLMEROverdispersion(model_1) # model_1 is dlightly overdispersed

model_1a <- glmer(Species_richness ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics, family = poisson)
summary(model_1a)
blmeco::dispersion_glmer(model_1a)
StatisticalModels::GLMEROverdispersion(model_1a)

model_1b <- glmer(Species_richness ~ LUI * zone + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612
summary(model_1b)
blmeco::dispersion_glmer(model_1b)
StatisticalModels::GLMEROverdispersion(model_1b)
StatisticalModels::R2GLMER(model_1b) # psuedo R squared

# model model_1b has lowest AIC - model_1b has the lowest AIC
AIC(model_1, model_1a, model_1b)

# checking null and simpler models
model_1b_intercept <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612
model_1b_LUI <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612

AIC(model_1b, model_1b_intercept, model_1b_LUI)

# build table of best fitting model
model_1b_table <- tidy(model_1b) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  mutate(term = gsub("zone", "zone-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:8) %>%
  gt() %>%
  tab_header(
    title = "Model summary - species richness (Figure 3)")

# save the model table for figure 2
write.table(model_1b_table, "figure_3_model_table_species-richness_exp.txt", sep = ",", row.names = FALSE, quote = FALSE)

## test for singularity - i.e. is the random effects structure too complicated - singularity od models is fine
tt <- getME(model_1b,"theta")
ll <- getME(model_1b,"lower")
min(tt[ll==0])

richness_metric <- predict_effects(iterations = 1000,
                                            model = model_1b, 
                                            model_data = pollinator_metrics, 
                                            response_variable = "Species_richness",
                                            fixed_number = 2,
                                            factor_number_1 = 2,
                                            factor_number_2 = 4,
                                            fixed_column = c("zone", "LUI"),
                                            neg_binom = FALSE)

# extract all the values for the data for a table of all values
model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(pollinator_metrics$LUI), 2)) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# values for richness for table
model_data(richness_metric)

### abundance
model_2 <- lmer(log(Total_abundance) ~ LUI * zone + (1|SS), data = pollinator_metrics) 
summary(model_2)
blmeco::dispersion_glmer(model_2)
StatisticalModels::GLMEROverdispersion(model_2)

model_2a <- lmerTest::lmer(log(Total_abundance) ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
summary(model_2a)
blmeco::dispersion_glmer(model_2a)
StatisticalModels::GLMEROverdispersion(model_2a)
StatisticalModels::R2GLMER(model_2a) # psuedo R squared

# check the AIC values for random effect structures - model_2a is lower
AIC(model_2, model_2a)

# checking fixed effects for that set of random effects
model_2a_intercept <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
model_2a_LUI <- lmer(log(Total_abundance) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

AIC(model_2a, model_2a_intercept, model_2a_LUI)

# build table of best fitting model
model_2a_table <- tidy(model_2a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  mutate(term = gsub("zone", "zone-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:8) %>%
  gt() %>%
  tab_header(
    title = "Model summary - total abundance (Figure 3)")

# save the model table for figure 2
write.table(model_2a_table, "figure_3_model_table_total_abundance_exp.txt", sep = ",", row.names = FALSE, quote = FALSE)

# predict values from the model                                
abundance_metric <- predict_effects(iterations = 1000,
                                             model = model_2a, 
                                             model_data = pollinator_metrics, 
                                             response_variable = "Total_abundance",
                                             fixed_number = 2,
                                             factor_number_1 = 2,
                                             factor_number_2 = 4,
                                             fixed_column = c("zone", "LUI"),
                                             neg_binom = FALSE)

# values for richness for table
model_data(abundance_metric)

### diversity
model_3 <- lmer(log(Simpson_diversity) ~ LUI * zone + (1|SS), data = pollinator_metrics)
summary(model_3)
blmeco::dispersion_glmer(model_3)
StatisticalModels::GLMEROverdispersion(model_3)

model_3a <- lmer(log(Simpson_diversity) ~ LUI * zone + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
blmeco::dispersion_glmer(model_3a)
StatisticalModels::GLMEROverdispersion(model_3a)

AIC(model_3, model_3a)

# checking fixed effects for that set of random effects
model_3a_intercept <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
model_3a_LUI <- lmer(log(Simpson_diversity) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

AIC(model_3a, model_3a_intercept, model_3a_LUI)

# build table of best fitting model
model_3a_table <- tidy(model_3a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  mutate(term = gsub("zone", "zone-", term)) %>%
  gt() %>%
  tab_header(
    title = "Figure 3 model summary - Simpson diversity")

# save the model table for figure 2
gtsave(model_3a_table, "figure_3_model_table_Simpson-diversity.png")

diversity_metric <- predict_effects(iterations = 1000,
                                             model = model_3a, 
                                             model_data = pollinator_metrics, 
                                             response_variable = "Simpson_diversity",
                                             fixed_number = 2,
                                             factor_number_1 = 2,
                                             factor_number_2 = 4,
                                             fixed_column = c("zone", "LUI"),
                                             neg_binom = FALSE)

# table of AICs
selection_table <- data.frame("Response" = c(rep("Species richness", 3),
                                             rep("Total abundance", 3),
                                             rep("Simpson diversity", 3)),
                              "Model" = c("Species_richness ~ LUI * zone (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Total_abundance ~ LUI * zone + (1|SS) + (1|SSB)", 
                                          "Total_abundance ~ LUI + (1|SS) + (1|SSB)", 
                                          "Total_abundance ~ 1 + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ LUI * zone + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ LUI + (1|SS) + (1|SSB)",
                                          "Simpson_diversity ~ 1 + (1|SS) + (1|SSB)"), 
                              "AIC" = c(AIC(model_1b), AIC(model_1b_LUI), AIC(model_1b_intercept), 
                                        AIC(model_2a), AIC(model_2a_LUI), AIC(model_2a_intercept), 
                                        AIC(model_3a), AIC(model_3a_LUI), AIC(model_3a_intercept))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# table for AIC values
gtsave(selection_table, "figure_2_model_table.png")

(richness_metric + xlab("") + ggtitle("A") + scale_x_discrete(labels = c("Non-tropical", "Tropical")) + guides(colour = guide_legend("Land-use intensity"))) +
(abundance_metric + ggtitle("B") + xlab("") + scale_x_discrete(labels = c("Non-tropical", "Tropical")) + guides(colour = FALSE)) + 
  plot_layout(ncol = 1) & 
  theme(text = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        aspect.ratio = 0.7/2)

ggsave("resampling_cropland_zone_8.pdf", dpi = 400, scale = 0.9)
