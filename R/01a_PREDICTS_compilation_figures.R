### packages
library(dplyr)
library(ggplot2)
library(forcats)
library(rworldmap)
library(rworldxtra)
library(viridis)
library(patchwork)
library(gt)

# source the functions R scripts
source("R/00_functions.R")

# read in rds for PREDICTS pollinators
PREDICTS_pollinators <- readRDS("outputs/PREDICTS_pollinators_8_exp.rds")

# remove additional columns - number of unique species = 5100
pollinator_subset <- PREDICTS_pollinators %>% 
  select(Kingdom, Phylum, Class, Order, Family, Best_guess_binomial, confidence, clade_rank) %>%
  filter(Best_guess_binomial != "") %>%
  filter(!grepl("\\d", Best_guess_binomial)) %>%
  unique()

# reorder the order factors for the plot
pollinator_subset$Order <- reorder(pollinator_subset$Order, pollinator_subset$Order, FUN = function(x) -length(x))
pollinator_subset$Class <- reorder(pollinator_subset$Class, pollinator_subset$Class, FUN = function(x) -length(x))
pollinator_subset$confidence <- factor(pollinator_subset$confidence, levels = c(1, 2, 3, 4, 5.1, 5.2, 5.3, 5.4))

# table for supplementary information for counts of species
# add `summarise(total = sum(n))` for total number of species
species_count_table <- pollinator_subset %>%
  group_by(Order, Class) %>%
  tally() %>% 
  ungroup() %>%
  gt() %>%
  tab_header(
    title = "Pollinating animal species (Figure 1)") %>%
  fmt_number(
    columns = vars(n),
    drop_trailing_zeros = TRUE,
    sep_mark = ""
  )

# build table of counts
pollinator_subset_tab <- pollinator_subset %>% 
  group_by(Order, Class) %>%
  tally() %>%
  ungroup()

# filter for upper
pollinat_group_upper <- pollinator_subset_tab %>%
  filter(n >= 14)

# filter for lower and assign "other order
pollinat_group_lower <- pollinator_subset_tab %>%
  filter(n < 14) %>%
  mutate(Order = "Other")

# bind lower and upper pollinat together
bound_group_pollinat <- rbind(pollinat_group_upper, pollinat_group_lower)

# build table for number of orders
# plot of species count for PREDICTS pollinator subset by class and order
order_class <- bound_group_pollinat %>%
  ggplot()+
    geom_bar(aes(x = Order, y = n, fill = Class), stat = "identity") +
    theme_bw() +
    ylab("Species count") +
    labs(subtitle = "B") +
    scale_y_continuous(limits = c(0, 2500), expand = c(0, 0)) +
    scale_fill_viridis(discrete = TRUE, option = "magma") +
    guides(fill = FALSE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(), 
          text = element_text(size = 13), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

# plot by class and confidence
class_conf <- pollinator_subset %>% 
  mutate(confidence = as.numeric(as.character(confidence))) %>%
  mutate(group = confidence > 5) %>% 
  group_by(group, Class) %>%
  tally() %>%
  ungroup() %>% 
  mutate(confidence = as.character(group)) %>%
  mutate(confidence = factor(group, levels = c("FALSE", "TRUE"), labels = c("Direct", "Extrapolated"))) %>%
  mutate(n = as.numeric(n)) %>%
  ggplot() +
    geom_bar(aes(x = confidence, fill = Class, y = n), stat = "identity", width = 0.75) +
    theme_bw() +
    ylab("") +
    labs(subtitle = "C") +
    xlab("Pollination confidence") +
    scale_fill_viridis(discrete = TRUE, option = "magma", na.translate = F) +
    scale_y_continuous(limits = c(0, 3500), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(), 
          axis.title.x = element_text(vjust = 14),
          text = element_text(size = 13))

### plot for geographic distribution of sites
# remove additional columns - number of unique species = 5100
pollinator_subset_geog <- PREDICTS_pollinators %>% 
  select(Latitude, Longitude, Use_intensity, Predominant_land_use) %>%
  unique() %>%
  dplyr::filter(Predominant_land_use != "Cannot decide") %>%
  dplyr::filter(Use_intensity != "Cannot decide") %>%
  dplyr::filter(Predominant_land_use != "Secondary vegetation (indeterminate age)") %>%
  mutate(Use_intensity = factor(Use_intensity, levels = c("Minimal use", "Light use", "Intense use"))) %>%
  mutate(Predominant_land_use = factor(Predominant_land_use, levels = c("Primary vegetation", "Mature secondary vegetation", 
                                                                        "Intermediate secondary vegetation", "Young secondary vegetation",
                                                                        "Plantation forest", "Pasture", "Cropland", "Urban")))

# calculate percentage in each UN region
PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use != "Cannot decide") %>%
  dplyr::filter(Use_intensity != "Cannot decide") %>%
  dplyr::filter(Predominant_land_use != "Secondary vegetation (indeterminate age)") %>%
  dplyr::select(SSBS, UN_region) %>%
  unique() %>%
  group_by(UN_region) %>%
  tally() %>%
  mutate(total = sum(n)) %>%
  filter(UN_region != "Americas") %>%
  bind_rows(data.frame("UN_region"= c("North America", "South America and the Caribbean"), "total" = 12170, "n" = c((2471 + 494), (112 + 1390)))) %>%
  mutate(percentage = (n/total * 100))

# build map
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# build map for distribution of sites
site_distribution <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, fill = Predominant_land_use), data = pollinator_subset_geog, shape = 21, colour = "white", alpha = 0.5) +
  scale_fill_manual("Predominant land-use", values = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "           A") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"))

# multiplot for map and histogram
site_distribution + (order_class + class_conf) + plot_layout(ncol = 1, widths = c(2), height = c(2, 1))

ggsave("conf_class_dist_PRED_14.pdf", dpi = 350, scale = 1.15)