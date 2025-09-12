# Load required libraries
library(tidyverse)   
library(microeco)    
library(reshape2)    
library(magrittr)    
library(patchwork)   

# Set working directory
setwd("~/Desktop/full_run/modeling")

# Load OTU table and metadata
otu_table <- read.csv("feature-table.csv", row.names = 1)

sample_info <- read.csv("kml_meta.csv", row.names = 1, check.names = FALSE)
sample_info <- cbind(SampleID = rownames(sample_info), sample_info)

taxonomy_table <- read.csv("taxonomy.csv", row.names = 1, header = TRUE)

# Load and process taxonomy table
tax_split <- colsplit(taxonomy_table$Taxon, ";", 
                      names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

tax_split[tax_split == ""] <- NA

tax <- cbind(taxonomy_table, tax_split)

tax$Confidence <- NULL

tax <- tax %>% 
  filter(!Kingdom %in% c("Unassigned", "_Mitochondrion", "_Nucleomorph", "k__") & Phylum != "Nematoda") %>%
  mutate(Class = ifelse(Class == "" | Class == "c__" | Class == "c__." | Class == ".", NA, Class)) %>%
  mutate(Class = coalesce(Class, Phylum, Kingdom))

tax %<>% tidy_taxonomy

tax2 <- tax %>%
  filter(!is.na(Kingdom))
tax2 <- rownames_to_column(tax2, var = "ASV")


# Align taxonomy with OTU table
original_tax <- tax
tax <- tax[rownames(otu_table), , drop = FALSE]

# Create microtable object
mt <- microtable$new(otu_table = otu_table, sample_table = sample_info, tax_table = tax)
mt$tidy_dataset()
mt

unique(mt$sample_table$treatment)
mt$sample_table$treatment <- trimws(mt$sample_table$treatment)

# Subset for stress and ambient treatments
mt$sample_table <- mt$sample_table %>% filter(treatment %in% c("stress", "ambient"))
mt$tidy_dataset()
mt

# Calculate relative abundance
mt$cal_abund()
mt$filter_taxa()
mt$save_abund(dirpath = "taxa_abund")

# plot taxa abundance 
my_colors2 <- c(
  "#a4b954",
  "#768b21",
  "#1e497c",
  "#497bac",
  "#8eb8d0",
  "#eab676",
  "#de7b48",
  "#70a5d7",
  "#4782cb",
  "#bc5d34",
  "#93472d",
  "#7f4438"
)


t1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 12)

taxa_plot <- t1$plot_bar(
  color_values = my_colors2, 
  others_color = "#c2c4c7", 
  facet = c("treatment", "timepoint"), 
  xtext_keep = FALSE, 
  legend_text_italic = FALSE)


print(taxa_plot)


#### making the helper plot based on metabolic strategy
# Load + reshape 
avg_long <- read.csv("top_12_taxa_avg.csv") %>%
  pivot_longer(
    cols = -c(treatment, timepoint),
    names_to = "Taxa",
    values_to = "Abundance"
  )

avg_long <- avg_long %>%
  mutate(treatment = factor(treatment, levels = c("ambient", "stress")))

# Metabolism mapping 
avg_long <- avg_long %>%
  mutate(Metabolism = case_when(
    Taxa %in% c("Bicoecea", "Ciliophora", "Flabellinia",
                "Labyrinthulomycetes", "Spirotrichea") ~ "heterotroph",
    Taxa %in% c("Chlorarachniophyceae", "Chlorophyta", "Chrysophyceae",
                "Dinophyceae", "Trebouxiophyceae") ~ "mixotroph",
    Taxa %in% c("Bacillariophyceae", "Mediophyceae") ~ "autotroph",
    TRUE ~ NA_character_
  ))

# Sum to metabolism × treatment × timepoint 
meta_sum <- avg_long %>%
  filter(!is.na(Metabolism)) %>%
  group_by(Metabolism, treatment, timepoint) %>%
  summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop")


# ensure ordered timepoints + numeric index for LOESS 
tp_levels <- meta_sum %>%
  mutate(tp_num = as.numeric(gsub("^T", "", as.character(timepoint)))) %>%
  arrange(tp_num) %>% pull(timepoint) %>% unique()

meta_sum <- meta_sum %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels),
         tp_idx    = as.numeric(timepoint))

# plot 
plot_loess_area <- function(dat, title, palette, span = 0.8, alpha_fill = 0.30, lw = 1.2) {
  ggplot(dat, aes(tp_idx, Total, color = treatment, fill = treatment)) +
    # shaded area under the LOESS curve (to y = 0)
    geom_ribbon(stat = "smooth", method = "loess", se = FALSE, span = span,
                aes(ymin = 0, ymax = after_stat(y)), alpha = alpha_fill) +
    # LOESS line
    geom_smooth(method = "loess", se = FALSE, span = span, linewidth = lw) +
    scale_x_continuous(breaks = seq_along(tp_levels), labels = tp_levels) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(title = title, x = "Timepoint", y = "Total abundance",
         color = "Treatment", fill = "Treatment") +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      panel.grid = element_blank()
    )
}

# colors for different groups 
pal_hetero <- c(ambient = "#9f4403", stress = "#cfa281")  # original you liked
pal_mixo   <- c(ambient = "#719da5", stress = "#B8CED2")  # example alt
pal_auto   <- c(ambient = "#516125", stress = "#a4b954")  # example alt

# build each plot SEPARATELY 
p_hetero <- meta_sum %>% filter(Metabolism == "heterotroph") %>%
  plot_loess_area("Heterotrophs", palette = pal_hetero)

p_mixo <- meta_sum %>% filter(Metabolism == "mixotroph") %>%
  plot_loess_area("Mixotrophs", palette = pal_mixo)

p_auto <- meta_sum %>% filter(Metabolism == "autotroph") %>%
  plot_loess_area("Autotrophs", palette = pal_auto)


all_metabolism <- p_hetero + p_mixo + p_auto
print(all_metabolism)
