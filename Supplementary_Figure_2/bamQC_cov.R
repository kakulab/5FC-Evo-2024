library(ggplot2)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggpubr)
library(patchwork)

setwd("Supplementary_Figure_2/Figure_S2_B_C")
# BAM coverage
## Cumulative coverage plot
cum_cov = read.delim("mosdepth-cumcoverage-dist-id.txt", sep = "\t", stringsAsFactors = F, check.names = F)
l_cum_cov = cum_cov |>
  pivot_longer(
    cols = !Sample,
    names_to = "coverage",
    values_to = "cov_pct"
  ) |>
  as.data.frame() |>
  mutate(
    pct = as.numeric(str_extract(cov_pct, "\\d+\\.\\d+")),
    coverage = as.numeric(coverage),
    Sample = case_when(
      Sample == "strain_185" ~ "R1",
      Sample == "strain_188" ~ "S1",
      Sample == "strain_189" ~ "R2",
      Sample == "strain_191" ~ "T1",
      TRUE ~ "WT"
    )
  )
cum_cov_p = ggplot(
  data = l_cum_cov,
  aes(
    x = coverage,
    y = pct,
    color = Sample
  )
) +
  geom_line(size = 1) +
  geom_vline(
    xintercept = c(100, 400),
    linetype = "dashed",
    alpha = 0.4
    ) + 
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = function(x){paste0(x, "%")}
  ) +
  scale_x_continuous(
    breaks = seq(0, 600, 100),
    labels = function(x){paste0(x, "X")}
  ) +
  theme_minimal() +
  scale_color_npg() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 24, vjust = 10),
    axis.title.y = element_text(size = 20),
    legend.position = "bottom",
    legend.justification = c(1.4,1.4),
    legend.text = element_text(size = 22)
  ) +
  labs(
    x = "Cumulative Coverage (X)",
    y = "% bases in genome\ncovered by at least X reads"
  )

## Average coverage per contig
avg_cov_contig = read.delim("mosdepth-coverage-per-contig-multi.txt", sep = "\t", header = TRUE, stringsAsFactors = F, check.names = F)
l_avg_cov_contig = avg_cov_contig |>
  pivot_longer(
    cols = !Sample,
    names_to = "no_scaffold",
    values_to = "scaff_cov"
  ) |>
  mutate(
    scaff_cov = str_replace_all(scaff_cov, "[()']", "")
  ) |>
  separate(
    scaff_cov,
    into = c("scaffold", "coverage"),
    sep = ", ",
    convert = TRUE
  ) |>
  mutate(
    scaffold = str_replace(scaffold, "\\.1", ""),
    # Refactor level using natural sorting
    scaffold = fct_relevel(scaffold, unique(str_sort(scaffold, numeric = T))),
    Sample = case_when(
      Sample == "strain_185" ~ "R1",
      Sample == "strain_188" ~ "S1",
      Sample == "strain_189" ~ "R2",
      Sample == "strain_191" ~ "T1",
      TRUE ~ "WT"
    )
  ) |>
  as.data.frame()
contig_cov_p = ggplot(
  data = l_avg_cov_contig,
  aes(
    x = scaffold,
    y = coverage,
    group = Sample,
    color = Sample,
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(
    breaks = seq(50, 350, 50),
    labels = function(x){paste0(x, "X")}
  ) +
  theme_minimal() +
  scale_color_npg() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
    legend.position = "none",
    axis.ticks.y = element_blank()
  ) +
  labs(
    y = "Average Coverage"
  )
cov_fig = (cum_cov_p + contig_cov_p)
cov_fig = cov_fig +
    plot_annotation(
        tag_levels = list(c("B", "C"))
) &
    theme(plot.tag = element_text(size = 24))
cov_fig
ggsave("coverage_summary_figs.pdf", plot = cov_fig, units = "in", width = 20, height = 10, dpi = 300, device = "pdf")
# ggsave("coverage_summary_figs.png", plot = cov_fig, units = "in", width = 20, height = 10, dpi = 300, device = "png")
