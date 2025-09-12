#################################################################################################
##                                                                                             ##
##                 Cell size / Density / Infection Models + Figure 2                           ##
##                                                                                             ##
#################################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(glmmTMB)
library(emmeans)
library(readr)

#################################################################################################
##                                                                                             ##
##                                  Setting up the data                                        ##
##                    This assumes the ASV modeling script has been ran also                   ##
##                                                                                             ##
#################################################################################################

# Define path to feature data
csv_dir <- "/Users/syrenawhitner/Desktop/KML_Scripts/IFCB_modeling/IFCB_feature_data"
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
ifcb_meta <- read.csv("/Users/syrenawhitner/Desktop/KML_Scripts/IFCB_modeling/IFCB_meta.csv") %>%
  mutate(across(c(treatment, sample_name), trimws),
         tank = substr(sample_name, 1, 1),
         timepoint = substr(sample_name, 2, 3))

# Load all CSVs and add sampleID
all_data_list <- lapply(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df <- df[df$Biovolume != 0, ]
  
  # Extract sampleID as first 16 characters to match IFCB_bin_ID
  df$sampleID <- substr(tools::file_path_sans_ext(basename(file)), 1, 16)
  
  return(df)
})

# Combine into single dataframe
all_data <- bind_rows(all_data_list)

# Map sampleID -> sample_name -> treatment/timepoint
meta_join <- ifcb_meta %>%
  select(IFCB_bin_ID, sample_name, treatment, timepoint, volume_analyzed_ml, tank)

all_data <- all_data %>%
  left_join(meta_join, by = c("sampleID" = "IFCB_bin_ID")) %>%
  filter(!is.na(sample_name)) %>%
  mutate(log_bio = log10(Biovolume),
         tank = substr(sample_name, 1, 1),
         timepoint = factor(timepoint),
         treatment = factor(treatment))

# loading up other relevant dataframes i.e parasite counts and density/infection information

# Subset all_data for the following models 
subset_data <- all_data %>%
  select(roi_number, Biovolume, sample_name, EquivDiameter, sampleID, treatment, timepoint, log_bio, numBlobs, tank, volume_analyzed_ml) %>%
  rename(sample = sample_name)

# relative to cell density 
total_cells <- subset_data %>%
  group_by(sample, treatment, timepoint, tank) %>%
  summarise(total_cell = sum(numBlobs, na.rm = TRUE), .groups = "drop") %>% # counts number of cells 
  left_join(ifcb_meta %>% select(sample_name, volume_analyzed_ml), by = c("sample" = "sample_name")) %>%
  mutate(cells_per_ml = total_cell / volume_analyzed_ml)

# relative to parasite counts
parasite_counts <- read.csv("/Users/syrenawhitner/Desktop/KML_Scripts/IFCB_modeling/file_counts.csv") 
parasite_counts$sampleID <- substr(parasite_counts$sample, 1, 16)

parasite_counts <- parasite_counts %>%
  left_join(meta_join, by = c("sampleID" = "IFCB_bin_ID")) %>%
  left_join(total_cells %>% select(sample, total_cell), by = c("sample_name" = "sample")) %>%
  mutate(
    prop_hosts = target / total_hosts,
    prop_total = target / total_cell,
    total_non_target = total_cell - target,
    rowid = row_number(),
    tank = substr(sample_name, 1, 1),
    treatment = factor(treatment),
    timepoint = factor(timepoint)
  )

#################################################################################################
##                                                                                             ##
##                                  cell size                                                  ##
##                                                                                             ##
#################################################################################################

# running the model 
biovolume_full_model <- glmmTMB(
  EquivDiameter ~ treatment * timepoint + 
    (1 | sample) +
    (1 | tank),
  family = gaussian, 
  data = subset_data
)

biovolume_reduced_model <- glmmTMB(
  EquivDiameter ~ treatment + timepoint + 
    (1 | sample) +
    (1 | tank),
  family = gaussian, 
  data = subset_data
)


lrt_biovol <- anova(biovolume_full_model, biovolume_reduced_model)
print(lrt_biovol)

em_biovol <- emmeans(biovolume_full_model, ~ treatment | timepoint)
biovol_contrasts <- as.data.frame(contrast(em_biovol, method = "pairwise", adjust = "fdr"))
write.csv(biovol_contrasts, "biovolume_contrasts.csv")
print(biovol_contrasts)

# plotting the voilins
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(data, ..., draw_quantiles = NULL, trim = TRUE) {
                             data <- transform(data,
                                               xminv = x - violinwidth * (x - xmin),
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- if (grp %% 2 == 1) {
                               transform(data, x = xminv)
                             } else {
                               transform(data, x = xmaxv)
                             }
                             newdata <- newdata[order(newdata$y), ]
                             GeomPolygon$draw_panel(newdata, ...)
                           }
)

# Wrapper function that accepts typical geom_violin arguments
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., draw_quantiles = NULL,
                              trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles,
                  na.rm = na.rm, ...)
  )
}

biovolume_violin <- ggplot(subset_data, aes(x = timepoint, y = EquivDiameter, fill = treatment)) +
  geom_split_violin(scale = "width", adjust = 1.2, color = "gray30", trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  scale_y_log10() +
  labs(
    x = "Timepoint",
    y = "Cell Diameter µm",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  ) +
  scale_fill_manual(values = c("ambient" = "#719da5", "stress" = "#dfd0c4"))

print(biovolume_violin)

#################################################################################################
##                                                                                             ##
##                                  Infection rates                                            ##
##                                                                                             ##
#################################################################################################

# infection model and comparing with LRT
# this tests prevalence of infection within the entire population 
infection_full_model <- glmmTMB(
  cbind(target, total_non_target) ~ treatment * timepoint + 
    (1 | tank) + 
    (1 | rowid),
  family = binomial, 
  data = parasite_counts
)

infection_reduced_model <- glmmTMB(
  cbind(target, total_non_target) ~ treatment + timepoint + 
    (1 | tank) + 
    (1 | rowid),
  family = binomial,  
  data = parasite_counts
)


lrt_infection <- anova(infection_full_model, infection_reduced_model)
print(lrt_infection)

# Post hoc emms comparisons 
em_infection <- emmeans(infection_full_model, ~ treatment | timepoint, type = "response")
infection_contrasts <- contrast(em_infection, method = "pairwise") %>%
  summary(infer = TRUE, adjust = "fdr")
print(infection_contrasts)

# plotting infection rates (total)
infection_summary_total <- parasite_counts %>%
  group_by(treatment, timepoint) %>%
  summarise(mean = mean(prop_total), se = sd(prop_total)/sqrt(n()), .groups = "drop") %>%
  mutate(time_numeric = as.numeric(gsub("T", "", timepoint)))

infection_summary_total <- infection_summary_total %>%
  mutate(
    ci_lower = mean - qt(0.975, df = 4 - 1) * se,  # 95% CI lower
    ci_upper = mean + qt(0.975, df = 4 - 1) * se   # 95% CI upper
  )


total_infection_rates <- ggplot(infection_summary_total, aes(x = time_numeric, y = mean, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean+ se),
                width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_color_manual(values = c("ambient" = "#719da5", "stress" = "#c9b6a6")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Prevalence of Infection (%)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  )

print(total_infection_rates)


# modeling for rate of infection within available hosts only
# running the models
host_infection_full_model <- glmmTMB(
  cbind(target, not_target) ~ treatment * timepoint + 
    (1 | tank) + 
    (1 | rowid),
  family = betabinomial,
  data = parasite_counts
)

host_infection_reduced_model <- glmmTMB(
  cbind(target, not_target) ~ treatment + timepoint + 
    (1 | tank) + 
    (1 | rowid),
  family = betabinomial,
  data = parasite_counts
)

lrt_host_infection <- anova(host_infection_full_model, host_infection_reduced_model)
print(lrt_host_infection)

# Post hoc emms comparisons 
em_infection_hosts <- emmeans(host_infection_full_model, ~ treatment | timepoint, type = "response")
infection_contrasts_hosts <- contrast(em_infection_hosts, method = "pairwise") %>%
  summary(infer = TRUE, adjust = "fdr")
print(infection_contrasts_hosts)

# plotting infection within hosts only 
parasite_summary_hosts <- parasite_counts %>%
  group_by(treatment, timepoint) %>%
  summarise(
    mean_infection = mean(prop_hosts, na.rm = TRUE),
    sd_infection = sd(prop_hosts, na.rm = TRUE),
    n = n(),
    se_infection = sd_infection / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    time_numeric = as.numeric(gsub("T", "", timepoint))  # extract numeric part of timepoint
  )

host_infection_rates <- ggplot(parasite_summary_hosts, aes(x = time_numeric, y = mean_infection, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_infection - se_infection,
                    ymax = mean_infection + se_infection),
                width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_color_manual(values = c("ambient" = "#719da5", "stress" = "#c9b6a6")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Prevalence of Infection (%)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  )

print(host_infection_rates)


# making a model for parasite cell density 
parasite_counts <- parasite_counts %>%
  mutate(parasites_per_ml = target / volume_analyzed_ml)

# model, LRT comparisons 
parasite_density_full_model <- glmmTMB(
  target ~ treatment * timepoint + 
    (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = parasite_counts
)

parasite_density_reduced_model <- glmmTMB(
  target ~ treatment + timepoint + 
    (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = parasite_counts
)

parasite_density_lrt <- anova(parasite_density_full_model, parasite_density_reduced_model)
print(parasite_density_lrt)

# Post hoc
em_parasite_density <- emmeans(parasite_density_full_model, ~ treatment | timepoint)
parasite_density_contrasts <- as.data.frame(contrast(em_parasite_density, method = "pairwise", adjust = "fdr"))
print(parasite_density_contrasts)

# plotting parasite cell density 
parasite_density_summary <- parasite_counts %>%
  group_by(treatment, timepoint) %>%
  summarise(
    mean_parasite_density = mean(parasites_per_ml, na.rm = TRUE),
    sd_parasite_density = sd(parasites_per_ml, na.rm = TRUE),
    n = n(),
    se_parasite_density = sd_parasite_density / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    time_numeric = as.numeric(gsub("T", "", timepoint))  # extract numeric part of timepoint
  )

parasite_density <- ggplot(parasite_density_summary, aes(x = time_numeric, y = mean_parasite_density, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_parasite_density - se_parasite_density,
                    ymax = mean_parasite_density + se_parasite_density),
                width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_color_manual(values = c("ambient" = "#719da5", "stress" = "#c9b6a6")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Parasite Density (cells/mL-1)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  )

print(parasite_density)


#################################################################################################
##                                                                                             ##
##                                  Cell density                                               ##
##                                                                                             ##
#################################################################################################

# running the models and comparing with LRT 
# this adjusts for volume analyzed per sample 
density_full_model <- glmmTMB(
  total_cell ~ treatment * timepoint + 
  (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = total_cells
)

density_reduced_model <- glmmTMB(
  total_cell ~ treatment + timepoint + 
    (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = total_cells
)

density_lrt <- anova(density_full_model, density_reduced_model)
print(density_lrt)

# Post hoc emms comparisons 
em_density <- emmeans(density_full_model, ~ treatment | timepoint)
density_contrasts <- as.data.frame(contrast(em_density, method = "pairwise", adjust = "fdr"))
print(density_contrasts)

# Density plot (total cells)
density_summary <- total_cells %>%
  group_by(treatment, timepoint) %>%
  summarise(mean = mean(cells_per_ml), se = sd(cells_per_ml)/sqrt(n()), .groups = "drop") %>%
  mutate(time_numeric = as.numeric(gsub("T", "", timepoint)))

total_density_plot <- ggplot(density_summary, aes(x = time_numeric, y = mean, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_y_log10() +
  scale_color_manual(values = c("ambient" = "#719da5", "stress" = "#c9b6a6")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Population Density (cells/mL-1)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  )

print(total_density_plot)


###### just looking at host density 
# Density plot (total cells)
host_counts <- parasite_counts %>%
  group_by(treatment, timepoint) %>%
  mutate(hosts_per_ml = total_hosts / volume_analyzed_ml) %>%
  summarise(mean = mean(hosts_per_ml), se = sd(hosts_per_ml)/sqrt(n()), .groups = "drop") %>%
  mutate(time_numeric = as.numeric(gsub("T", "", timepoint)))

# model for host density
host_full_model <- glmmTMB(
  total_hosts ~ treatment * timepoint + 
    (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = parasite_counts
)

host_reduced_model <- glmmTMB(
  total_hosts ~ treatment + timepoint + 
    (1 | tank),
  family = nbinom2,
  offset = log(volume_analyzed_ml),
  data = parasite_counts
)

host_density_lrt <- anova(host_full_model, host_reduced_model)
print(host_density_lrt)

# Post hoc emms comparisons 
em_hosts <- emmeans(host_full_model, ~ treatment | timepoint)
hosts_contrasts <- as.data.frame(contrast(em_hosts, method = "pairwise", adjust = "fdr"))
print(hosts_contrasts)

host_density_plot <- ggplot(host_counts, aes(x = time_numeric, y = mean, color = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_color_manual(values = c("ambient" = "#719da5", "stress" = "#c9b6a6")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Population Density (cells/mL-1)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()  
  )

print(host_density_plot)







