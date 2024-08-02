library(ggplot2)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("others/library_complex")

# Complexity curve and DupRadar
## Complexity curve
preseq_stat = read.csv("preseq_data.csv", header = T, stringsAsFactors = F, check.names = F)
l_preseq_stat = preseq_stat |>
  rename(total_mol = colnames(preseq_stat)[1]) |>
  pivot_longer(
    cols = c(-1),
    values_to = "counts",
    names_to = "strain"
  ) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  )
complex_curve = ggplot(
    data = l_preseq_stat,
    aes(
        x = total_mol,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain),
    linewidth = 1
  ) +
  geom_segment(
    aes(
      x = 0, xend = 6.88,
      y = 0, yend = 6.88
    ),
    linetype = 2
  ) +
  annotate(
    geom = "richtext",
    size = 6,
    x = 40, y = 6.80,
    label = "_a perfect library where each read is unique_<br>**6.88 M unique molecules**"
  ) +
  scale_x_continuous(breaks = seq(0, 140, 20),
                     labels = function(x){paste0(x, "M")}) +
  scale_y_continuous(breaks = seq(0, 8, 1),
                     labels = function(x){paste0(x, "M")}) +
  scale_color_npg() +
  labs(
    x = "Total molecules (including duplicates)",
    y = "Unique molecules"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
## DupRadar
dup_stat = read.csv("mqc_hcplot_bsqlomaevp.csv", header = T, stringsAsFactors = F)
l_dup_stat = dup_stat |>
  rename(expression = colnames(dup_stat)[1]) |>
  pivot_longer(
    cols = starts_with("strain"),
    values_to = "counts",
    names_to = "strain"
  ) |>
  mutate(
    strain = case_when(
      strain == "strain_185" ~ "185-R",
      strain == "strain_188" ~ "188-S",
      strain == "strain_189" ~ "189-R",
      strain == "strain_191" ~ "191-I",
      TRUE ~ "WT-S"
      )
  ) |>
  filter(!is.na(counts))
dupradar_p = ggplot(
    data = l_dup_stat,
    aes(
        x = expression,
        y = counts,
        group = strain
    )
) +
  geom_line(
    aes(color = strain),
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(0.5, 1000),
    linetype = "dashed",
    color = c("#367E18", "#DF2E38")
    ) +
  # log10 scaling and label with exponential format
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_color_npg() +
  labs(
    x = "expression (reads/kbp)",
    y = "% duplicate reads"
  ) +
  theme_pubclean() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  )
complex_dupradar = complex_curve / dupradar_p
complex_dupradar = complex_dupradar +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 18))
# ggsave("complex_dupradar_figs.png", plot = complex_dupradar, units = "in", width = 14, height = 14, dpi = 300, device = "png")
ggsave("complex_dupradar_figs.pdf", plot = complex_dupradar, units = "in", width = 14, height = 14, dpi = 300, device = "pdf")
