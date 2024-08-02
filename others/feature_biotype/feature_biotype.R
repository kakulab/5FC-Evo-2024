library(ggplot2)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("others/feature_assignment")

# Feature + Biotype count
## Biotype counts
bio_count = read.csv("featureCounts_biotype_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_bio_count = bio_count |>
  rename(strain = colnames(bio_count)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "type"
  ) |>
  mutate(
    type = factor(str_replace(type, "_", " ")),
    type = fct_relevel(
      type,
      c("Unassigned Ambiguity", "Unassigned NoFeatures", "tRNA pseudogene", "rRNA", "tRNA", "protein coding")
    ),
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
    )
  )
p_biotype = ggplot(
    data = l_bio_count,
    aes(
        x = strain,
        y = counts,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("#A91D3A", "#7469B6", "orange", "green", "black", "#7BC9FF")
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = function(x){paste0(x * 100, "%")}) +
  labs(
    y = "Percentages"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
## Feature assignment
feature_count = read.csv("featureCounts_assignment_plot.csv", header = T, stringsAsFactors = F, check.names = F)
l_feature_count = feature_count |>
  rename(strain = colnames(feature_count)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "type"
  ) |>
  mutate(
    type = factor(str_replace(type, "_", " ")),
    type = fct_relevel(
      type,
      c("Unassigned Ambiguity", "Unassigned NoFeatures", "Unassigned MultiMapping", "Assigned")),
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  )
p_feature = ggplot(
    data = l_feature_count,
    aes(
        x = strain,
        y = counts,
        fill = type
    )
) +
  geom_bar(
    position = "fill",
    stat = "identity"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = function(x){paste0(x*100, "%")}) +
  coord_flip() +
  scale_fill_manual(
    values = c("orange", "#90ed7d", "#322C2B", "#7BC9FF")
  ) +
  labs(
    y = "Percentages"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
feature_fig = p_biotype / p_feature
feature_fig = feature_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 18))
# ggsave("feature_biotype_assignment.png", plot = feature_fig, units = "in", width = 10, height = 6, dpi = 300, device = "png")
ggsave("feature_biotype_assignment.pdf", plot = feature_fig, units = "in", width = 10, height = 6, dpi = 300, device = "pdf")