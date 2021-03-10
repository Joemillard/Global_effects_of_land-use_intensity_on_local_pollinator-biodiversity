### packages for analysis,and suppress messages
# have checked the removal of the braconidae, which has little effect - keep them in
# combine vertebrates in with the invertebrates
packages <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# extra packages called in place -- plyr and StatisticalModels

# source in additional functions
source("Scripts/global_analysis/Land-use_intensity_predicts_differential_effects_on_global_pollinator_biodiversity/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS(here::here("outputs/PREDICTS_pollinators_8_exp.rds"))

PREDICTS_pollinators_intense <- PREDICTS_pollinators

### taxa - use intensity
PREDICTS_pollinators_intense$Order_use <- paste(PREDICTS_pollinators_intense$Order, PREDICTS_pollinators_intense$Use_intensity, sep = "-")
PREDICTS_pollinators_intense$LUI <- paste(PREDICTS_pollinators_intense$Predominant_land_use, PREDICTS_pollinators_intense$Use_intensity, sep = "-")

# filter cannot decide factors
PREDICTS_pollinators_intense <- PREDICTS_pollinators_intense %>%
  dplyr::filter(Use_intensity != "Cannot decide") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  dplyr::filter(LUI != "Primary vegetation-Intense use") %>%
  dplyr::filter(LUI != "Primary vegetation-Light use") %>%
  droplevels()

# reassign for cropland and primary vegetation factors
PREDICTS_pollinators_intense <- PREDICTS_pollinators_intense %>%
  mutate(LUI = plyr::revalue(LUI, c("Primary vegetation-Minimal use" = "Primary vegetation",
                                                  "Cropland-Minimal use" = "Cropland (Minimal use)",
                                                  "Cropland-Light use" = "Cropland (Light use)",
                                                  "Cropland-Intense use" = "Cropland (Intense use)")))

# assign new variable for tropical/temperate, convert to factor, and filter out NA
PREDICTS_pollinators_intense$zone <- ifelse(PREDICTS_pollinators_intense$Latitude >= -23.5 & PREDICTS_pollinators_intense$Latitude <= 23.5, "Tropics", "Temperate")
PREDICTS_pollinators_intense$zone <- factor(PREDICTS_pollinators_intense$zone, levels = c("Temperate", "Tropics"))

# filter the metrics for those with greater than 50 sites
PREDICTS_pollinators_intense <- PREDICTS_pollinators_intense[PREDICTS_pollinators_intense$Order %in% c("Hymenoptera", "Lepidoptera", "Diptera", "Coleoptera", "Passeriformes", "Apodiformes"),]

# drop levels after filtering for taxa of interest
PREDICTS_pollinators_intense <- droplevels(PREDICTS_pollinators_intense)

# correct for sampling effort
PREDICTS_pollinators_intense <- CorrectSamplingEffort(PREDICTS_pollinators_intense)

# create object as list of pollinator subsets, and then use apply to calculate metrics
# split diversity data into list of four for each order
diversityOrder <- split(x = PREDICTS_pollinators_intense, f = PREDICTS_pollinators_intense$Order)

# drop unused levels from each list
diversityOrder <- lapply(diversityOrder, function(x) return(droplevels(x)))

# Remove empty rows
diversityOrder <- Filter(function(x) dim(x)[1] > 0, diversityOrder)

# calculate site metrics for each of the four order level subsets
order.sites.div <- do.call('rbind',lapply(X = diversityOrder, FUN = SiteMetrics,
                                          extra.cols = c("SSB", "SSBS","Biome", "Sampling_method",
                                                         "Study_common_taxon", "Sampling_effort",
                                                         "Sampling_effort_unit", "Realm",
                                                         "Use_intensity", "Order", "LUI", "zone"),
                                          sites.are.unique = TRUE, srEstimators = TRUE))

# reorder the factors
order.sites.div$LUI <- factor(order.sites.div$LUI, levels = c("Primary vegetation", "Cropland (Minimal use)", "Cropland (Light use)", "Cropland (Intense use)"))
order.sites.div$Order <- factor(order.sites.div$Order, levels = c("Hymenoptera", "Diptera", "Lepidoptera", "Coleoptera", "Passeriformes", "Apodiformes"))

# add 1 to abundance and diversity 
order.sites.div$Total_abundance <- order.sites.div$Total_abundance + 1
order.sites.div$Simpson_diversity <- order.sites.div$Simpson_diversity + 1

# build table for site representation
figure_4_table <- order.sites.div %>%
  group_by(Order, LUI) %>%
  tally() %>%
  gt() %>%
  tab_header(
    title = "Global cropland site representation (Figure 4)") %>%
  fmt_number(
    columns = vars(n),
    drop_trailing_zeros = TRUE,
    sep_mark = ""
  )

### models - richness (with interaction)
model_1 <- glmer(Species_richness ~ Order * LUI + (1|SS), data = order.sites.div, family = poisson)
summary(model_1)
StatisticalModels::GLMEROverdispersion(model_1)
AIC(model_1) #

model_1a <- glmer(Species_richness ~ Order * LUI + (1|SS) + (1|SSB), data = order.sites.div, family = poisson)
summary(model_1a)
StatisticalModels::GLMEROverdispersion(model_1a)
AIC(model_1a) #

model_1b <- glmer(Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS), data = order.sites.div, family = poisson)
summary(model_1b)
StatisticalModels::GLMEROverdispersion(model_1b)
AIC(model_1b) # best model

AIC(model_1, model_1a, model_1b)

model_1b_LUI <- glmer(Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS), data = order.sites.div, family = poisson)
model_1b_int <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = order.sites.div, family = poisson)

AIC(model_1b, model_1b_LUI, model_1b_int)

# build table of best fitting model
model_1b_table <- tidy(model_1b) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:23) %>%
  gt() %>%
  tab_header(
    title = "Model summary - species richness (Figure 4)")

# predict corrected effects for species richness
richness_metric <- predict_effects(iterations = 1000, 
                                   model = model_1b, 
                                   model_data = order.sites.div, 
                                   response_variable = "Species_richness",
                                   fixed_number = 2,
                                   fixed_column = c("Order", "LUI"),
                                   factor_number_1 = 6,
                                   factor_number_2 = 4,
                                   neg_binom = FALSE)

# extract all the values for the data for a table of all values
model_data <- function(model_plot){
  ggplot_build(richness_metric)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(order.sites.div$LUI), 6)[1:23]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric)

## abundance models
model_2 <- lmer(log(Total_abundance) ~ Order * LUI + (1|SS), data = order.sites.div)
summary(model_2)

model_2a <- lmerTest::lmer(log(Total_abundance) ~ Order * LUI + (1|SS) + (1|SSB), data = order.sites.div)
summary(model_2a)

AIC(model_2, model_2a) # best model

model_2a_LUI <- lmer(log(Total_abundance) ~ LUI + (1|SS) + (1|SSB), data = order.sites.div)
model_2a_int <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = order.sites.div)

AIC(model_2a, model_2a_LUI, model_2a_int)

# build table of best fitting model
model_2a_table <- tidy(model_2a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:23) %>%
  gt() %>%
  tab_header(
    title = "Model summary - total abundance (Figure 4)")

# predict corrected effects for total abundance
abundance_metric <- predict_effects(iterations = 1000, 
                                    model = model_2a, 
                                    model_data = order.sites.div, 
                                    response_variable = "Total_abundance",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 6,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# extract all the values for the data for a table of all values
model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(order.sites.div$LUI), 6)[1:23]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric)

## diversity models
model_3 <- lmer(log(Simpson_diversity) ~ Order * LUI + (1|SS), data = order.sites.div)
summary(model_3)

model_3a <- lmerTest::lmer(log(Simpson_diversity) ~ Order * LUI + (1|SS) + (1|SSB), data = order.sites.div)
summary(model_3a)

AIC(model_3, model_3a) # best model

model_3a_LUI <- lmer(log(Simpson_diversity) ~ LUI + (1|SS) + (1|SSB), data = order.sites.div)
model_3a_int <- lmer(log(Simpson_diversity) ~ 1 + (1|SS) + (1|SSB), data = order.sites.div)

AIC(model_3a, model_3a_LUI, model_3a_int)

# build table of best fitting model
model_3a_int <- tidy(model_3a) %>% 
  mutate(term = gsub("sd_\\(Intercept\\)\\.", "", term)) %>%
  mutate(term = gsub("Order", "Order-", term)) %>%
  mutate(term = gsub("LUI", "LUI-", term)) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
  slice(1:23) %>%
  gt() %>%
  tab_header(
    title = "Model summary - Simpson diversity (Figure 4)")

# predict corrected effects for diversity models
diversity_metric <- predict_effects(iterations = 1000, 
                                    model = model_3a, 
                                    model_data = order.sites.div, 
                                    response_variable = "Simpson_diversity",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 6,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# abundance data
model_data(diversity_metric)

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
                              "AIC" = c(AIC(model_1b), AIC(model_1b_LUI), AIC(model_1b_int), 
                                        AIC(model_2a), AIC(model_2a_LUI), AIC(model_2a_int),
                                        AIC(model_3a), AIC(model_3a_LUI), AIC(model_3a_int))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# scripts for multi plot
(richness_metric + xlab(NULL) + guides(colour = FALSE) + ggtitle("A") + scale_y_continuous("Species richness diff. (%)") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())) + 
(abundance_metric + ggtitle("B") + xlab(NULL) + scale_y_continuous("Total abundance diff. (%)") + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 10), legend.title = element_text(size = 11)) + guides(colour = guide_legend("Land-use intensity"))) + 
(diversity_metric + ggtitle("C") + xlab(NULL) + guides(colour = FALSE) + scale_y_continuous("Simpson diversity diff. (%)") + theme(axis.text.x = element_text(size = 10))) + plot_layout(ncol = 1)

ggsave("land-use_intensity_response_all-taxa_10.pdf", scale = 1.1, dpi = 350)

