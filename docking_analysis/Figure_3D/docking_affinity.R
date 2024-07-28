library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggpol)
library(ggtext)
library(rstatix)

setwd("docking_analysis/Figure_3D")

affi = read.table("docking_analysis/autodock_vina/binding_affinity.csv", sep = ",", header = T, stringsAsFactors = F, check.names = F)
affi$affinity = abs(affi$affinity)
affi$receptor = factor(affi$receptor, levels = c("WT", "R1", "R4"))

# Barplot for best binding affinity comparison between WT, R1 and R4
affi |>
    group_by(receptor, ligand) |>
    filter(affinity == max(affinity)) |>
    mutate(
        ligand = case_when(
            ligand == "5FU" ~ "5FC",
            ligand == "PRPP" ~ "PRPP",
            ligand == "PRPP_5FU" ~ "PRPP + 5FC"
        )
    ) |>
    ggplot(
        aes(x = receptor, y = affinity, fill = ligand)
    ) +
        geom_bar(
            position = "dodge",
            stat = "identity"
        ) +
        geom_text(
            aes(label = round(affinity,2)),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3,
        ) +
        facet_wrap(~ligand) +
        scale_fill_manual(values = c("#4361cf", "#FFAF00", "#346d33")) +
        scale_y_continuous(
            labels = function(x){paste0("-", x)}
        ) +
        labs(
            y = "Best binding affinity (kcal/mol)"
        ) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, color = "black"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 12)
        )
ggsave("binding_affinity_score.pdf", unit = "in", width = 6, height = 5, device = "pdf", dpi = 300)
# ggsave("binding_affinity_score.png", unit = "in", width = 6, height = 5, device = "png", dpi = 300)
