library(vcfR)
library(tidyverse)
library(ggsci)
library(ggtext)
library(gghalves)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

setwd("./others/vcfQC")

unfiltered_vcf = read.vcfR("cauris.snp.filtered.PEKT.ann.vcf", verbose = FALSE)

# Create a full chromR object
unfiltered_chrom = create.chromR(name = "", vcf = unfiltered_vcf, verbose = TRUE)
unfiltered_chrom = masker(unfiltered_chrom, min_QUAL = 1, min_DP = 30)
unfiltered_chrom = proc.chromR(unfiltered_chrom, verbose = TRUE)

# Summarize chromR object
head(unfiltered_chrom)

# Extract DP stats in FORMAT column (extract.gt)
dp <- extract.gt(unfiltered_chrom, element = "DP", as.numeric = TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)
is.na(dp[na.omit(dp == 0)]) <- TRUE

# The sequence variant depth for each variant per sample
variant_dp_per_sample_df = dp |>
    as.data.frame() |>
    pivot_longer(
        cols = starts_with("strain"),
        values_to = "DP",
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

# High variant depths between all 5 strains
variant_DP_p = ggplot(
    data = variant_dp_per_sample_df,
    aes(
        x = strain,
        y = DP,
        fill = strain
    )
) +
    geom_half_violin(side = "l") +
    geom_half_boxplot(
        side = "r", width = 0.2,
        position = position_nudge(x = 0.05)
    ) +
    geom_half_point_panel(
        aes(color = strain),
        position = position_nudge(x = 0.05)
    ) +
    scale_y_continuous(
        transform = scales::log2_trans(),
        breaks = c(1, 10, 100, 1000, 8000),
        minor_breaks = c(1:10, 2:10*10, 2:10*100, 2:8*1000)
    ) +
    stat_compare_means(
        aes(group = strain),
        label = "p.signif",
        comparisons = list(
            c("185-R", "WT-S"),
            c("188-S", "WT-S"),
            c("189-R", "WT-S"),
            c("191-I", "WT-S")
        )
    ) +
    labs(
        x = "strain", y = "Variant depth"
    ) +
    theme_minimal() +
    scale_fill_npg() +
    scale_color_npg() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# ggsave("img/variant_depth_across_samples.png", unit = "in", width = 20, height = 8, device = "png", dpi = 300)

# Missingness across all samples
missingness_per_sample = dp |>
    as_tibble() |>
    summarise(
        across(everything(), ~ sum(is.na(.x)))
    ) |>
    pivot_longer(
        cols = starts_with("strain"),
        values_to = "missingness",
        names_to = "strain"
    ) |>
    mutate(
        miss_pct = missingness/nrow(unfiltered_vcf) *100,
        strain = case_when(
            strain == "strain_185" ~ "185-R",
            strain == "strain_188" ~ "188-S",
            strain == "strain_189" ~ "189-R",
            strain == "strain_191" ~ "191-I",
            TRUE ~ "WT-S"
        )
    )
missingness_per_sample_p = ggplot(
    data = missingness_per_sample,
    aes(
        x = strain, y = miss_pct,
        fill = strain,
        label = miss_pct
    )
) +
    geom_col() +
    geom_text(
        aes(label = round(miss_pct, 2)),
        nudge_y = 0.05,
        size = 8
    ) +
    labs(
        y = "Missingness percent (%)"
    ) +
    theme_classic() +
    scale_fill_npg() +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5, face ="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# ggsave("img/missingness_per_sample.png", unit = "in", width = 10, height = 6, device = "png", dpi = 300)

# Missingness across variants
missingness_per_variant = dp |>
    as_tibble() |>
    rowwise() |>
    mutate(
        missingness = sum(is.na(c_across(everything())))
    ) |>
    ungroup() |>
    mutate(
        miss_pct = missingness/ncol(unfiltered_vcf@gt[, -1]) * 100
    ) |>
    select(!contains("strain"))
missingness_per_variant_dp = ggplot(
    data = missingness_per_variant,
    aes(miss_pct)
) +
    geom_histogram(
        aes(y = ..density..),
        fill = "steelblue",
        color = "black", alpha = 0.6
    ) +
    geom_density(
        alpha = 0.6,
        color = "#a063b7",
        fill = "#a063b7"
    ) +
    labs(
        x = "Missingness (%)"
    ) +
    theme_classic() +
    scale_fill_npg() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(hjust = 1, size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, vjust = 10),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
# ggsave("img/missingness_across_variants.png", unit = "in", width = 10, height = 6, device = "png", dpi = 300)

# Combined plots
varDP_missingness_fig = variant_DP_p / (missingness_per_sample_p | missingness_per_variant_dp)
varDP_missingness_fig = varDP_missingness_fig +
    plot_annotation(
        tag_levels = "A"
) &
    theme(plot.tag = element_text(size = 24))
varDP_missingness_fig
# ggsave("varDP_missingness_figs.png", plot = varDP_missingness_fig, units = "in", width = 16, height = 10, dpi = 300, device = "png")
ggsave("varDP_missingness_figs.pdf", plot = varDP_missingness_fig, units = "in", width = 16, height = 10, dpi = 300, device = "pdf")
