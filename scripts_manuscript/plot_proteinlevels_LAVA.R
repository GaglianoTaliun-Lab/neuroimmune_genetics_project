# Description: summary plots for LAVA results with pQTLs

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(cowplot)
library(gridExtra)
library(grid)
library(ggraph)

# Arguments ------------------------------------------------------

project_dir = "Documents/research-projects/neuro_immune"

date=str_remove_all(Sys.Date(), "-")

cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# pQTLs
# "Number of total bivariate tests: 1863. The bonferroni corrected p-value = 2.68384326355341e-05."
# "Number of total bivariate tests: 1863. The FDR corrected p-value = 0.0016264."

FDR_thres = 0.00164396

# 0.05 divided by the number of total bivariate tests performed
bonferroni_thres = 0.05/1863

# Load files -----------------------------------------------------------

sign_results_fdr <- read.table(here(project_dir, "lava_results", "pqtls", "FDR_significant_pQTLs_bivar_window_100000_20220803.tsv"), sep = "\t", header = T)
sign_results_bonf <- read.table(here(project_dir, "lava_results","pqtls", "bonferroni_significant_pQTLs_bivar_window_100000_20220803.tsv"), sep = "\t", header = T)
all_results <- read.table(here(project_dir, "lava_results", "pqtls", "pQTLs_bivar_window_100000_20220803.tsv"), sep = "\t", header = T)

sample_sizes <- read.table(here(project_dir, "lava_results", "N_GWAS_traits.txt"), sep = "\t", header = F)
colnames(sample_sizes) <- c("phen1_clean","N")

# Get list of genes for each trait -------------------------------------

list_trait <- sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  group_split(., phen1_clean)

lapply(list_trait, function(x) {
  
  write.table(x$gene_locus, here(project_dir, "lava_results", "pqtls", "genes_per_trait",str_c(x$phen1_clean[1],"_list_of_genes.txt")),
              sep = "\t", row.names = F, quote = F, col.names = F)
  
})

# Main plots -----------------------------------------------------------

### plot: bar plot of summary results

summary_counts <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", "gene"), sep = "_") %>%
  group_by(., phen1_clean) %>%
  count(., phen2_clean) %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","LBD","PD","ALS","SCZ","MS","UC","CD"))
  ) %>%
  left_join(., sample_sizes, by = "phen1_clean") %>%
  mutate(phen1_plot = str_c(order_phen1," \n(",N,")"))

# plot code:
ggplot(data = summary_counts, aes(x = reorder(phen1_plot, N), y = n)) +
  geom_col(fill = "#CCCCCC", colour = "black") +
  theme_bw() +
  labs(x = "Trait (GWAS sample size)", y = "Number of significant \n local genetic correlations") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 20)) +
  geom_text(aes(label = n, vjust = -0.2), size= 8) +
  ylim(0, 85)

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "pqtls", stringr::str_c("summary_of_significant_local_correlations_pqtls_FDR_",date,".jpg")),
       width = 40, height = 30, units = "cm")
# --------------------------------------------------------------------------------------------------------

### plot: count plot for FDR significant results

# - count genes
bivar_to_plot <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", "gene"), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","LBD","PD","ALS","SCZ","MS","UC","CD")),
         fill_rho = case_when(
           p <= FDR_thres ~ round(rho, 2)
         ),
         chr_gene = str_c("chr ",chr,": ", gene_locus)
  )

# - create and save plot:
ggplot(bivar_to_plot, aes(x=phen1_clean, y=reorder(gene_locus, start))) + 
  geom_count(colour = "black", shape = 15, size = 8) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "Trait", y = "Expressed Genes", title = "") +
  guides(colour = guide_colourbar(title = "Local rg")) +
  theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size =20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20)) +
  coord_flip() +
  facet_wrap(~chr, nrow = 4, scales = "free")

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "pqtls", stringr::str_c("FDR_local_genetic_correlations_pqtls_",date,".jpg")),
       width = 80, height = 60, units = "cm")

# --------------------------------------------------------------------------------------------------------

### plot: count plot for FDR significant results, only those significant for a test trait and at least one control trait

# - count genes
bivar_to_plot <-
  all_results %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", "gene"), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","LBD","PD","ALS","SCZ","MS","UC","CD")),
         fill_rho = case_when(
           p <= 0.05 ~ round(rho, 2)
         ),
         chr_gene = str_c("chr ",chr,": ", gene_locus),
         sig = case_when(
           p < FDR_thres ~ "**"
         ),
         nom_sig = case_when(
           p < 0.05 & p >= FDR_thres ~ ".",
         ),
         chr_gene = str_c("chr ",chr,": ", gene_locus),
         not_sig = case_when(
           p >= 0.05 ~ "ns"
         )
  )

# list of significant genes for ND
list_genes_ND <- sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  filter(., phen1_clean %in% c("AD","PD","LBD","ALS","FTD","SCZ")) %>%
  dplyr::select(gene_id)

# list of significant genes for immune-related diseases
list_genes_imm <- sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  filter(., phen1_clean %in% c("MS","UC","CD")) %>%
  dplyr::select(gene_id)

# overlap across lists
list_overlap <- list_genes_ND %>%
  filter(., gene_id %in% list_genes_imm$gene_id)

# plot_genes
bivar_to_plot2 <- bivar_to_plot %>%
  filter(., gene_id %in% list_overlap$gene_id)

# reorder based on chr, start position:
bivar_to_plot2 <- bivar_to_plot2[with(bivar_to_plot2, order(chr, start)),]
bivar_to_plot2$chr_gene <- factor(bivar_to_plot2$chr_gene, levels = unique(bivar_to_plot2$chr_gene))

# - create and save plot:
ggplot(bivar_to_plot2, aes(x=order_phen1, y=chr_gene)) + 
  geom_count(colour = "black", shape = 15, size = 8) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "GWAS Trait", y = "Protein gene symbol", title = "") +
  guides(colour = guide_colourbar(title = "Local rg")) +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size =20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20)) +
  geom_text(aes(label=sig), size = 4, vjust = 0.75) +
  geom_text(aes(label=nom_sig), size = 5, vjust = 0.1) +
  geom_text(aes(label=not_sig), size = 4, vjust = 0.4) +
  coord_flip()

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "pqtls", stringr::str_c("FDR_shared_local_genetic_correlations_pqtls_",date,".jpg")),
       width = 20, height = 15, units = "cm")
