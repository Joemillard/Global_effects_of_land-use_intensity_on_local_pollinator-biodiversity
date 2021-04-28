### packages for analysis,and suppress messages
packages <- c("dplyr", "yarg", "cowplot", "gt", "webshot", "lme4", "broom.mixed", "MASS", "roquefort")
lapply(packages, require, character.only = TRUE)

# extra packages called in place -- plyr and StatisticalModels

source("R/00_functions.R")

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
                                  extra.cols = c("SSB", "SSBS", "LUI", "UN_region", "Country", "Source_for_predominant_land_use"),
                                  sites.are.unique = TRUE,
                                  srEstimators = TRUE)

# check the distribution of factors for intensity and type
table(pollinator_metrics$LUI)

# convert to factor and reorder factor levels
pollinator_metrics$LUI <- factor(pollinator_metrics$LUI, levels = c("Primary vegetation-Minimal use", 
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
# drop left-over factor levels
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
pollinator_metrics$Simpson_diversity <- pollinator_metrics$Simpson_diversity + 1

# build table for site representation
figure_2_table <- pollinator_metrics %>%
  group_by(LUI) %>%
  tally() %>%
  gt() %>%
  tab_header(
    title = "Global site representation (Figure 2)") %>%
  fmt_number(
    columns = vars(n),
    drop_trailing_zeros = TRUE,
    sep_mark = ""
  )

### richness - because only 1 fixed effect, select model random effect structure on basis of AIC
model_1 <- glmer(Species_richness ~ LUI + (1|SS), data = pollinator_metrics, family = poisson)
summary(model_1)
StatisticalModels::GLMEROverdispersion(model_1) # model_1 is slightly overdispersed

model_1a <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics, family = poisson)
summary(model_1a)
StatisticalModels::GLMEROverdispersion(model_1a)

model_1b <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612
summary(model_1b)
anova(model_1b)
StatisticalModels::GLMEROverdispersion(model_1b)

# check the AIC values - random structure from model_1b is lowest
AIC(model_1, model_1a, model_1b)

# with the given random effect structure, check null and simpler models
model_1b_intercept <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = pollinator_metrics, family = poisson) # best model - failed to converge with 0.00144612

AIC(model_1b, model_1b_intercept)

# build table of best fitting model
model_1b_table <- tidy(model_1b) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:22) %>%
  gt() %>%
  tab_header(
    title = "Model summary - species richness (Figure 2)")

# derive predictions for the richness model
richness_metric <- predict_effects(iterations = 1000, 
                                   model = model_1b, 
                                   model_data = pollinator_metrics, 
                                   response_variable = "Species_richness",
                                   fixed_number = 1,
                                   fixed_column = "LUI",
                                   neg_binom = FALSE)

# new approach for calculating percentage change following comments from reviewers
((10*1.26152911 - 10/1.19782141) / 10) * 100 # urban -- 42.66801
((10*1.27581931 - 10/1.11060222) / 10) * 100 # plantation -- 37.54069
((10*1.02093625 - 10/1.15724824) / 10) * 100 # YSV -- 15.68174

# extract all the values for the data for a table of all values for source data
model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = levels(pollinator_metrics$LUI)) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# extract the values for species richness for supplementary data
model_data(richness_metric)
  
### abundance
model_2 <- lmer(log(Total_abundance) ~ LUI + (1|SS), data = pollinator_metrics) 
summary(model_2)

# lmer with log transformation
model_2a <- lmerTest::lmer(log(Total_abundance) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics) # best model
summary(model_2a)

# check the AIC values - random structure from model_2a is lowest
AIC(model_2, model_2a)

# check null model isn't better with those random effects
model_2a_intercept <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

# build table of best fitting model
model_2a_table <- broom::tidy(model_2a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  mutate(term = gsub("sd_Observation\\.", "", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:22) %>%
  gt() %>%
  tab_header(
    title = "Model summary - total abundance (Figure 2)")

AIC(model_2a, model_2a_intercept)

# derive predictions for the abundance model
abundance_metric <- predict_effects(iterations = 1000, 
                                    model = model_2a, 
                                    model_data = pollinator_metrics, 
                                    response_variable = "Total_abundance",
                                    fixed_number = 1,
                                    fixed_column = "LUI",
                                    neg_binom = FALSE)

# new approach for calculating percentage change following comments from reviewers
((10*1.5355206012 - 10/1.0874067897) / 10) * 100 # urban -- 61.59016
((10*1.7948945791 - 10*1.0481701986) / 10) * 100 # pasture -- 74.67244

# extract the values for species richness for supplementary data
model_data(abundance_metric)

### diversity
model_3 <- lmer(log(Simpson_diversity) ~ LUI + (1|SS), data = pollinator_metrics)
summary(model_3)

model_3a <- lmer(log(Simpson_diversity) ~ LUI + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

# check the AIC values - random structure from model_3a is lowest
AIC(model_3, model_3a)

# check null isn't better for those random effects
model_3a_intercept <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = pollinator_metrics) # best model

# build table of best fitting model
model_3a_table <- tidy(model_3a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  mutate(term = gsub("sd_Observation\\.", "", term)) %>%
  gt() %>%
  tab_header(
    title = "Figure 2 model summary - Simpson diversity")

AIC(model_3a, model_3a_intercept)

# table of AICs
selection_table <- data.frame("Response" = c(rep("Species richness", 2),
                                             rep("Total abundance", 2),
                                             rep("Simpson diversity", 2)),
                              "Model" = c("Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)", 
                                          "Total_abundance ~ LUI + (1|SS) + (1|SSB)", 
                                          "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                          "Simpson_diversity ~ LUI + (1|SS) + (1|SSB)", 
                                          "Simpson_diversity ~ 1 + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(model_1b), AIC(model_1b_intercept), 
                                      AIC(model_2a), AIC(model_2a_intercept),
                                      AIC(model_3a), AIC(model_3a_intercept))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# combine the plots for Simpson diversity and total abundance into single figure
plot_grid((richness_metric + ggtitle("A") + theme(axis.title.y = element_text(size = 14), 
                                                  legend.background = element_blank(), 
                                                  legend.text = element_text(size = 13))), 
          (abundance_metric + guides(shape = FALSE) + ggtitle("B") + theme(axis.title.y = element_text(size = 14))) + 
          NULL,
          NULL,
          ncol = 2)

# save the final plot
ggsave("all_metrics_resampling_approach_7.pdf", dpi = 400, scale = 1.5)
