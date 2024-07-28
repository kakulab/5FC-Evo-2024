library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ggforce)
library(ggrepel)
library(ggplotify)

setwd("./Figure_2B/")
# Input multiple samples vcf file
## Using "PASS" hard filtering and overall DP >= 30
snp_filtered_vcf = "cauris.snps.PASS.DP30.PEKT.ann.vcf"

## Reformat from vcf to gds object
snpgdsVCF2GDS(snp_filtered_vcf, "cauris.snps.PASS.DP30.PEKT.ann.gds", method = "biallelic.only")

## Summary
snpgdsSummary("cauris.snps.PASS.DP30.PEKT.ann.gds")

# Open a GDS file
snp_filtered_gds = snpgdsOpen("cauris.snps.PASS.DP30.PEKT.ann.gds")

# ---LD-based SNP pruning---

## Different LD thresholds for sensitivity analysis
set.seed(134)   # For reproducibility
filtered_snpset = snpgdsLDpruning(snp_filtered_gds, ld.threshold = 0.2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)
str(filtered_snpset)
names(filtered_snpset)

## Get all selected snp id
filtered_snpset_id = unlist(unname(filtered_snpset))
filtered_snpset_id

# ---Principal Component Analysis (PCA)---

## Pruning SNPs with LD threshold = 0.2 to exclude highly correlated SNPs
filtered_pca <- snpgdsPCA(snp_filtered_gds, snp.id = filtered_snpset_id, num.thread = 2, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN)

filtered_pca_df <- data.frame(sample.id = filtered_pca$sample.id,
    EV1 = filtered_pca$eigenvect[,1],
    EV2 = filtered_pca$eigenvect[,2],
    stringsAsFactors = FALSE)
filtered_pca_df = filtered_pca_df |> mutate(
    phenotype = case_when(
        sample.id == "strain_188" ~ "188-S",
        sample.id == "strain_185" ~ "185-R",
        sample.id == "strain_189" ~ "189-R",
        sample.id == "strain_191" ~ "191-I",
        TRUE ~ "WT-S"
    )
) |> mutate(
    sample.id = str_replace(sample.id, "_", " ")
)
head(filtered_pca_df)

## Variance proportion (%)
filtered_pc_percent <- filtered_pca$varprop*100
head(round(filtered_pc_percent, 2))

## PCA plot using pruned SNPs
set.seed(123)
PCA_vcf_p = ggplot(
    data = filtered_pca_df,
    aes(
        x = EV1, y = EV2,
    )
) +
    geom_point(
        aes(color = factor(phenotype)),
        position = position_jitter(width = 0.15, height = 0.15),
        size = 6, alpha = 1
    ) +
    # xlim(c(-0.8, 0.9)) +
    scale_x_continuous(limits = c(-0.75, 0.8),
                       breaks = seq(-0.75, 0.8, by = 0.5)) +
    labs(
        x = "PC1 (62.5%)", y = "PC2 (37.5%)",
    ) +
    theme_classic() +
    scale_color_npg() +
    scale_fill_npg() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "top",
        # legend.justification = c(1,1)
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)  # Adjust the size of legend point label
        )
    )
# ggsave("./PCA_prunedSNPs.pdf", plot = PCA_vcf_p, units = "in", width = 6, height = 5, dpi = 300, device = "pdf")
ggsave("./PCA_prunedSNPs.png", plot = PCA_vcf_p, units = "in", width = 6, height = 5, dpi = 300, device = "png")
