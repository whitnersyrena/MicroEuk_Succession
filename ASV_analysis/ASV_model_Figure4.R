##### This script runs the entire treatment:timepoint:ASV model, timepoint contrasts, and then generates Figure 4. Important outputs: 
# top_12_labeled_ASV_plot was generated and indivudal panals were selected and modified in inkscape to create figure 4. 
# most_responsive2 is the most responsive ASVs when averaging over time 

library(tidyverse)  
library(glmmTMB)     
library(reshape2)    

#########################################################################################################
#
#        FULL MODEL VERSION
#
#########################################################################################################

# if running from Mac 
setwd("~/Desktop/full_run/modeling")

# -load data
feature_table <- read.csv("feature-table.csv")
metadata <- read.csv("kml_meta.csv")
metadata$treatment <- trimws(metadata$treatment)
taxonomy <- read.csv("taxonomy.csv")

feature_table_OG <- feature_table
metadata_OG <- metadata
taxonomy_OG <- taxonomy

# Rename ASV column
feature_table <- feature_table %>% dplyr::rename(ASV = `OUT_ID`)

feature_table2 <- feature_table %>%
mutate(
ASV_ID = paste0("ASV", row_number())
)

# Convert feature table to long format
feature_table_long <- feature_table %>%
pivot_longer(cols = -"ASV", names_to = "sample_name", values_to = "reads")

# Merge feature table with metadata
merged_data <- feature_table_long %>%
left_join(metadata, by = "sample_name")

# Merge with taxonomy information
merged_data <- merged_data %>%
left_join(taxonomy, by = c("ASV" = "Feature.ID"))

# Split Taxon column into separate taxonomic levels
taxonomy_split <- colsplit(merged_data$Taxon, ";",
                 names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Add taxonomy split back to merged data
merged_data <- cbind(merged_data, taxonomy_split) %>%
mutate(across(13:19, ~ ifelse(. == "", NA, .))) %>%     # Convert empty taxonomic ranks to NA
filter(Kingdom != "Unassigned" & Kingdom != "_Mitochondrion" & Kingdom != "k__" & Phylum != "Nematoda", )    # Remove unassigned or mitochondrial classifications

merged_data <- merged_data %>%
mutate(Class = ifelse(Class == "" | Class == "c__" | Class == "c__.", NA, Class)) %>%
mutate(Class = coalesce(Class, Phylum, Kingdom))

# just making a backup
merged_data_OG <- merged_data

# Remove control samples
merged_data <- merged_data %>%
filter(!treatment %in% c("control_ambient", "control_stress")) %>%
filter(!sample_name %in% c("P1.EXT.NEG", "P1.NEG.1", "P1.NEG.2", "P1.POS.1"))

# Remove low-abundance ASVs (present in fewer than threshold samples)
threshold <- 6  

merged_data <- merged_data %>%
group_by(Class) %>%
filter(sum(reads > 0) > threshold) %>%
ungroup()

# Compute Read Proportions #
# Calculate total reads per sample
merged_data <- merged_data %>%
group_by(sample_name) %>%
mutate(total_reads = sum(reads)) %>%
ungroup()

# Compute reads proportion per ASV
merged_data <- merged_data %>%
mutate(reads_proportion = reads / total_reads)

#set reference levels
merged_data <- merged_data %>%
mutate(
treatment = relevel(factor(treatment), ref = "ambient"),
timepoint = relevel(factor(timepoint), ref = "T0"),
ASV = factor(ASV),
tank = factor(tank)
)

# fit/run model
final_model.nodisp <- glmmTMB(
reads_proportion ~ treatment * timepoint +
(1 | tank) +
(1 | ASV) +
(1 | tank:ASV) +
(1 | timepoint:ASV) +
(1 | treatment:ASV) +
(1 | treatment:timepoint:ASV),
family = betabinomial,
data = merged_data,
weights = total_reads,
dispformula = ~ (1 | ASV))


final_model.nodisp2 <- glmmTMB(
reads_proportion ~ treatment * timepoint +
(1 | tank) +
(1 | ASV) +
(1 | tank:ASV) +
(1 | timepoint:ASV) +
(1 | treatment:ASV),
family = betabinomial,
data = merged_data,
weights = total_reads,
dispformula = ~ (1 | ASV))

final_anova = anova(final_model.nodisp, final_model.nodisp2)
print(final_anova)
summary(final_anova)

# extract random effects  
random_effects <- ranef(final_model.nodisp, condVar = TRUE) %>% as.data.frame()

random_effects_OG <- random_effects

#########################################################################################################
#
#        LOOKING AT RANEF VALUES OVER TIME
#
#########################################################################################################
# SECTION #1 - two way interaction - which ASVs are most impacted when averging over time? 
# filter for treatment:ASV interaction
treatment_asv_effects <- random_effects_OG %>%
filter(grpvar == "treatment:ASV")  

# extract treatment labels
# treatment names in the grp column are labeled as "ambient:ASV1" and "stress:ASV1"
treatment_asv_effects <- treatment_asv_effects %>%
separate(grp, into = c("Treatment", "ASV"), sep = ":", extra = "merge") %>%
filter(component == "cond")  

# compute difference in ranef values between treatments (not including timepoints here)
asv_changes2 <- treatment_asv_effects %>%
select(ASV, Treatment, condval) %>%
pivot_wider(names_from = Treatment, values_from = condval) %>%
mutate(effect_change = (stress - ambient))  # true values

asv_changes2_sorted <- asv_changes2 %>%
arrange(desc(effect_change))

most_responsive <- asv_changes2_sorted %>%
slice(c(1:50, 325:375, 699:749))

most_responsive <- most_responsive %>%
left_join(
merged_data %>% select(ASV, Taxon) %>% distinct(),
by = "ASV"
)

response_tax_split <- colsplit(most_responsive$Taxon, ";",
                     names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

most_responsive2 <- cbind(most_responsive, response_tax_split)

most_responsive2 <- most_responsive2 %>%
mutate(across(c(Genus, Family, Order, Class, Phylum, Kingdom),
      ~ ifelse(.x == "" | .x == ".", NA, .x))) %>%
mutate(Genus = coalesce(Genus, Family, Order, Class, Phylum, Kingdom))

most_responsive2 <- most_responsive2 %>%
left_join(feature_table2 %>% select(ASV, ASV_ID), by = "ASV") %>%
mutate(label = ifelse(is.na(Genus) | Genus == "", ASV_ID, paste0(Genus, " (", ASV_ID, ")")))

write.csv(most_responsive2, "most_responsive_asv.csv")

# identify ASVs with largest changes -- just doing all of them
top_asv_changes2 <- asv_changes2 %>%
arrange(desc(effect_change)) %>%
head(100)  

print(top_asv_changes2)

# Section #2 - three way interaction - which ASVs have a variable response to treatment, over time 
# filter for treatment:timepoint:ASV interaction
three_way_effects <- random_effects %>%
filter(grpvar == "treatment:timepoint:ASV")

# separate treatment, timepoint, ASV
three_way_effects <- three_way_effects %>%
separate(grp, into = c("Treatment", "Timepoint", "ASV"), sep = ":", extra = "merge")  

# calculate difference between treatments at each timepoint
treatment_differences <- three_way_effects %>%
select(ASV, Timepoint, Treatment, condval) %>%
pivot_wider(names_from = Treatment, values_from = condval) %>%
mutate(treatment_diff = stress - ambient)  # fifference in condval between treatments for each ASV

# treatment_differences
asv_variance <- treatment_differences %>%
group_by(ASV) %>%
summarize(variance_in_diff = var(treatment_diff, na.rm = TRUE)) %>%
arrange(desc(variance_in_diff))  # sort by highest variance

print(asv_variance)

# subset for ASvs with highest variance
top_asvs <- asv_variance %>%
arrange(desc(variance_in_diff)) %>%
slice(1:75) %>%  # select the top 25 ASVs with the most treatment variance
pull(ASV)

# filter data for these ASVs
plot_data <- three_way_effects %>%
filter(ASV %in% top_asvs)

# PLOT - just getting a sense of whatthe data looks like, can skip
plot1 <- ggplot(plot_data, aes(x = Timepoint, y = condval, color = Treatment, group = Treatment)) +
geom_line() +
geom_point() +
facet_wrap(~ ASV) +
labs(title = " top 75 ASV condval",
x = "Timepoint",
y = "Random Effect Value") +
scale_color_manual(values = c("ambient" = "#58bdd8", "stress" = "#f8564c")) +  # Set colors
theme_minimal()


# Section #3 -  adding the time averaged values from step 1 to the values in step 2         
### ---- adding time averaged ranef values to specific timepoints combinations ---- ###
mod_diff <- treatment_differences
mod_diff$treatment_diff = NULL

mod_diff_updated <- mod_diff %>%
  left_join(asv_changes2 %>% select(ASV, ambient), by = "ASV", suffix = c("", "_top")) %>%
  mutate(ambient_new = ambient + ambient_top)  

mod_diff_updated2 <- mod_diff_updated %>%
  left_join(asv_changes2 %>% select(ASV, stress), by = "ASV", suffix = c("", "_top")) %>%
  mutate(stress_new = stress + stress_top)  

mod_diff2 = mod_diff_updated2[, -c(3, 4, 5, 7)]

mod_diff2 <- mod_diff2 %>%
  rename(
    ambient = ambient_new,
    stress = stress_new
  ) %>%
  mutate(treatment_diff = stress - ambient)

mod_diff2_long <- mod_diff2 %>%
  pivot_longer(cols = c(ambient, stress),
               names_to = "Treatment",    
               values_to = "condval")      

## redoing previous steps with the new time averaged dataset
asv_variance2 <- mod_diff2_long %>%
  group_by(ASV) %>%
  summarize(variance_in_diff = var(treatment_diff, na.rm = TRUE)) %>%
  arrange(desc(variance_in_diff))  # Sort by highest variance

# subset for ASvs with highest variance
top_asvs2 <- asv_variance2 %>%
  arrange(desc(variance_in_diff)) %>%
  slice(1:75) %>%  
  pull(ASV)

# filter data for specific ASVs
plot_data2 <- mod_diff2_long %>%
  filter(ASV %in% top_asvs2)

# plot - just generally looking at all 75 of the top most variable ASVs with the timepoint averaged values added 
plot2 <- ggplot(plot_data2, aes(x = Timepoint, y = condval, color = Treatment, group = Treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ ASV) +
  labs(title = " ",
       x = "Timepoint",
       y = "Random Effect Value") +
  scale_color_manual(values = c("ambient" = "#58bdd8", "stress" = "#f8564c")) +  # Set colors
  theme_minimal()

plot_data_tax <- plot_data2 %>%
  left_join(
    merged_data %>%
      select(ASV, Taxon) %>%
      distinct(ASV, .keep_all = TRUE),
    by = "ASV"
  )

plot_tax_split <- colsplit(plot_data_tax$Taxon, ";",
                           names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

plot_data_tax2 <- cbind(plot_data_tax, plot_tax_split)

plot_data_tax3 <- plot_data_tax2 %>%
  mutate(across(c(Genus, Family, Order, Class, Phylum, Kingdom),
                ~ ifelse(.x == "" | .x == ".", NA, .x))) %>%
  mutate(Genus = coalesce(Genus, Family, Order, Class, Phylum, Kingdom))

plot_data_tax3 <- plot_data_tax3 %>%
  left_join(feature_table2 %>% select(ASV, ASV_ID), by = "ASV")

# create a fallback label if Genus is missing
plot_data_tax3 <- plot_data_tax3 %>%
  mutate(label = ifelse(is.na(Genus) | Genus == "", ASV_ID, paste0(Genus, " (", ASV_ID, ")")))

# Plot using new label for facets -- same as previous plot just with the labels 
plot3 <- ggplot(plot_data_tax3, aes(x = Timepoint, y = condval, color = Treatment, group = Treatment)) +
  geom_hline(yintercept = 0, color = "lightgrey", linewidth = 0.4) +  # Add horizontal reference line
  geom_line() +
  geom_point() +
  facet_wrap(~ label) +
  labs(title = " ",
       x = "Timepoint",
       y = "Random Effect Value") +
  scale_color_manual(values = c("ambient" = "#00acc1", "stress" = "#ff8a65")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )

write.csv(plot_data_tax3, "responsive_asv_time.csv")

##################################################################################################
# looking at RANEF of different trophic modes
plot_tax_4 <- read.csv("plot_tax_4.csv")

# Split based on trophic mode
autotrophs <- plot_tax_4 %>% filter(mode == "autotroph")
heterotrophs <- plot_tax_4 %>% filter(mode == "heterotroph")
mixotrophs  <- plot_tax_4 %>% filter(mode == "mixotroph")

###### variance based ordering 
#mixotrophs
var_order_mixotrophs <- mixotrophs %>%
  group_by(label) %>%
  summarise(condval_var = var(condval, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(condval_var))

top12_mixotrophs <- var_order_mixotrophs %>%
  slice(1:12) %>%
  pull(label)

mixotrophs_top12_var <- mixotrophs %>%
  filter(label %in% top12_mixotrophs) %>%
  left_join(var_order_mixotrophs, by = "label") %>%
  mutate(label = forcats::fct_reorder(label, condval_var))

#autotrophs 
var_order_autotrophs <- autotrophs %>%
  group_by(label) %>%
  summarise(condval_var = var(condval, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(condval_var))

top12_autotrophs <- var_order_autotrophs %>%
  slice(1:12) %>%
  pull(label)

autotrophs_top12_var <- autotrophs %>%
  filter(label %in% top12_autotrophs) %>%
  left_join(var_order_autotrophs, by = "label") %>%
  mutate(label = forcats::fct_reorder(label, condval_var))

# heterotrophs 
var_order_heterotrophs <- heterotrophs %>%
  group_by(label) %>%
  summarise(condval_var = var(condval, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(condval_var))

top12_heterotrophs <- var_order_heterotrophs %>%
  slice(1:12) %>%
  pull(label)

heterotrophs_top12_var <- heterotrophs %>%
  filter(label %in% top12_heterotrophs) %>%
  left_join(var_order_heterotrophs, by = "label") %>%
  mutate(label = forcats::fct_reorder(label, condval_var))


plot_tax_5 <- bind_rows(
  heterotrophs_top12_var,
  mixotrophs_top12_var,
  autotrophs_top12_var
)

top_12_labeled_ASV_plot <- ggplot(plot_tax_5, aes(x = Timepoint, y = condval, color = Treatment, group = Treatment)) +
  geom_hline(yintercept = 0, color = "lightgrey", linewidth = 0.4) +  # Add horizontal reference line
  geom_line() +
  geom_point() +
  facet_wrap(~ label) +
  labs(title = "Random Effect Changes Over Time for Selected ASVs",
       x = "Timepoint",
       y = "Random Effect Value") +
  scale_color_manual(values = c("ambient" = "#008c9d", "stress" = "#9f4403")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )

print(top_12_labeled_ASV_plot)

##################################################################################################
#
#             LOOPING THROUGH EACH TIMEPOINT TO DO CONTRASTS
#
##################################################################################################

#### recalculating read proportions for each model first.... do I need to do this?
# subset data for timepoint 0
subset_data0 <- merged_data %>% filter(timepoint == "T0")

# model for timepoint1
model_timepoint0v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data0,
  weights = total_reads,  
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint0v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data0,
  weights = total_reads,  
  dispformula = ~ treatment +
    (1 | ASV)
)

t0_anova = anova(model_timepoint0v1, model_timepoint0v2)
print(t0_anova)

# Subset data for timepoint 1
subset_data1 <- merged_data %>% filter(timepoint == "T1")

# Fit the model for timepoint1
model_timepoint1v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data1,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint1v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data1,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

t1_anova = anova(model_timepoint1v1, model_timepoint1v2)
print(t1_anova)

# Subset data for timepoint 2
subset_data2 <- merged_data %>% filter(timepoint == "T2")

# Fit the model for timepoint2
model_timepoint2v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data2,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint2v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data2,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

t2_anova = anova(model_timepoint2v1, model_timepoint2v2)
print(t2_anova)

# Subset data for timepoint 3
subset_data3 <- merged_data %>% filter(timepoint == "T3")

# Fit the model for timepoint3
model_timepoint3v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data3,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint3v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data3,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

t3_anova = anova(model_timepoint3v1, model_timepoint3v2)
print(t3_anova)

# Subset data for timepoint 4
subset_data4 <- merged_data %>% filter(timepoint == "T4")

# Fit the model for timepoint4
model_timepoint4v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data4,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint4v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data4,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

t4_anova = anova(model_timepoint4v1, model_timepoint4v2)
print(t4_anova)

# Subset data for timepoint 5
subset_data5 <- merged_data %>% filter(timepoint == "T5")

# Fit the model for timepoint5
model_timepoint5v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data5,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint5v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data5,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

t5_anova = anova(model_timepoint5v1, model_timepoint5v2)
print(t5_anova)

# Subset data for timepoint 6
subset_data6 <- merged_data %>% filter(timepoint == "T6")

# Fit the model for timepoint6
model_timepoint6v1 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV) +
    (1 | treatment:ASV),
  family = betabinomial,
  data = subset_data6,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)

model_timepoint6v2 <- glmmTMB(
  reads_proportion ~ treatment +
    (1 | ASV),
  family = betabinomial,
  data = subset_data6,
  weights = total_reads,
  dispformula = ~ treatment +
    (1 | ASV)
)


t6_anova = anova(model_timepoint6v1, model_timepoint6v2)
print(t6_anova)


contrast_results <- list()
contrast_results[[1]] <- data.frame(timepoint = 0, t0_anova)
contrast_results[[2]] <- data.frame(timepoint = 1, t1_anova)
contrast_results[[3]] <- data.frame(timepoint = 2, t2_anova)
contrast_results[[4]] <- data.frame(timepoint = 3, t3_anova)
contrast_results[[5]] <- data.frame(timepoint = 4, t4_anova)
contrast_results[[6]] <- data.frame(timepoint = 5, t5_anova)
contrast_results[[7]] <- data.frame(timepoint = 6, t6_anova)
time_contrasts <- do.call(rbind, contrast_results)

write.csv(time_contrasts, "time_contrasts.csv")



