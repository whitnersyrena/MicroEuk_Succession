#################################################################################################
##                                                                                             ##
##                            qPCR analysis - Figure 5A                                        ##
##                                                                                             ##
#################################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# set working directory and read data 
setwd("/Users/syrenawhitner/Desktop/full_run/qPCR_results")
df <- read.csv("quantities.csv", stringsAsFactors = FALSE)

unique(df$Sample)

# remove controls and negatives 
controls <- c("P1_EXT_NEG","C1","C2","C3","C4","C5","C6","C7","NTC","DC1","DC2")
filt_df <- df %>% filter(!Sample %in% controls)

# keep only the averaged rows every third observation 
third <- filt_df %>% slice(seq(3, n(), by = 3))

# parse tank id and timepoint 
third <- third %>%
  mutate(
    tank      = substr(Sample, 1, 1),                      # first character
    timepoint = stringr::str_extract(Sample, "T\\d+")      # e.g. T0, T1, ...
  )

# assign tanks to groups 
group_mapping <- list(
  "1" = c("A","B","H","K","L"), # stress group
  "2" = c("D","E","F","I","J"), # ambient grop 
  "3" = c("C"),
  "4" = c("G")
)

third <- third %>%
  mutate(
    group = case_when(
      tank %in% group_mapping[["1"]] ~ "1",
      tank %in% group_mapping[["2"]] ~ "2",
      tank %in% group_mapping[["3"]] ~ "3",
      tank %in% group_mapping[["4"]] ~ "4",
      TRUE                           ~ NA_character_
    )
  ) %>%
  filter(group %in% c("1","2"))  # keep only groups 1 and 2

# summarize means and error by group timepoint and target 
averages3 <- third %>%
  group_by(group, timepoint, Target) %>%
  summarise(
    mean_gc = mean(Quantity.Mean, na.rm = TRUE),
    sd_gc   = sd(Quantity.Mean,  na.rm = TRUE),  # preserves your original approach
    n       = n(),
    se_gc   = sd_gc / sqrt(n),
    .groups = "drop"
  )

# create numeric time for consistent spacing on x axis 
averages3_fmt <- averages3 %>%
  mutate(
    time_numeric = readr::parse_number(timepoint),               
    timepoint    = factor(timepoint, levels = paste0("T", 0:9))  
  )

# plot gene copies over time by group/target 
gc_plot <- ggplot(
  averages3_fmt,
  aes(x = time_numeric, y = mean_gc, color = group, group = interaction(Target, group))
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_gc - se_gc, ymax = mean_gc + se_gc), width = 0.2) +
  scale_x_continuous(breaks = 0:6, labels = paste0("T", 0:6)) +
  scale_color_manual(values = c("1" = "#c9b6a6", "2" = "#719da5")) +
  labs(
    title = " ",
    x = "Timepoint",
    y = "Gene Copies",
    color = "Group"
  ) +
  facet_wrap(~ Target, ncol = 2, scales = "free_y") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gc_plot

#################################################################################################
##                                                                                             ##
##             dual axis plot for cell density and bact. abundance - Figure 5b                 ##
##                                                                                             ##
#################################################################################################

# prep data --> this uses the density_summary dataframe from the cell density modeling
avg16s <- averages3 %>%
  mutate(
    treatment = recode(group, "1" = "stress", "2" = "ambient", .default = NA_character_),
    treatment = factor(treatment, levels = c("stress","ambient"))
  ) %>%
  filter(Target == "16S", !is.na(treatment))

dens_tbl <- density_summary %>%
  mutate(
    time_numeric = if ("time_numeric" %in% names(density_summary)) time_numeric else readr::parse_number(as.character(timepoint)),
    treatment    = factor(treatment, levels = c("stress","ambient"))
  ) %>%
  arrange(time_numeric)

tp_levels <- dens_tbl %>% pull(timepoint) %>% unique()

dens_tbl <- dens_tbl %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels),
         tp_idx    = as.numeric(timepoint))

avg16s <- avg16s %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels),
         tp_idx    = as.numeric(timepoint))

# dual-axis mapping 
rng1 <- range(dens_tbl$mean, na.rm = TRUE)         
rng2 <- range(avg16s$mean_gc, na.rm = TRUE)        
stopifnot(diff(rng2) != 0)
b <- (rng1[2] - rng1[1]) / (rng2[2] - rng2[1])
a <- rng1[1] - b * rng2[1]

# combine for plotting 
plot_df <- bind_rows(
  dens_tbl %>% transmute(tp_idx, timepoint, treatment, series = "Density", value = mean),
  avg16s   %>% transmute(tp_idx, timepoint, treatment, series = "16S",     value = a + b * mean_gc)
) %>%
  mutate(group = paste(series, treatment, sep = "_"))

# aesthetics 
col_map <- c(
  "Density_ambient" = "#444279",
  "Density_stress"  = "#969acb",
  "16S_ambient"     = "#005e1c",
  "16S_stress"      = "#92c37d"
)
linetype_map <- c(
  "Density_ambient" = "solid",
  "Density_stress"  = "dashed",
  "16S_ambient"     = "solid",
  "16S_stress"      = "dashed"
)
legend_order <- c("Density_ambient","Density_stress","16S_ambient","16S_stress")
label_map <- c(
  "Density_ambient" = "Density (ambient)",
  "Density_stress"  = "Density (stress)",
  "16S_ambient"     = "16S (ambient)",
  "16S_stress"      = "16S (stress)"
)

# plot 
y_top <- ceiling(max(dens_tbl$mean, na.rm = TRUE) / 1000) * 1000

p_dual <- ggplot(plot_df, aes(tp_idx, value, color = group, linetype = group)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq_along(tp_levels), labels = tp_levels) +
  scale_y_continuous(
    name   = "Cell density",
    limits = c(0, y_top),
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::label_number(big.mark = ","),
    sec.axis = sec_axis(~ (. - a) / b,
                        name   = "16S abundance",
                        labels = scales::label_number(big.mark = ","))
  ) +
  scale_color_manual(values = col_map, breaks = legend_order, labels = label_map[legend_order]) +
  scale_linetype_manual(values = linetype_map, breaks = legend_order, labels = label_map[legend_order]) +
  labs(x = "Timepoint", title = "") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

p_dual






