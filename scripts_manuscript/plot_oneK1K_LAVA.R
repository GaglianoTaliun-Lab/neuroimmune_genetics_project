# Description: summary plots for LAVA results with OneK1Ké eQTLs

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(cowplot)
library(gridExtra)
library(grid)
library(ggraph)
library(ggpubr)
library(RColorBrewer)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/neuro_immune"

date=str_remove_all(Sys.Date(), "-")

cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ONEK1K
# "Number of total bivariate tests: 1628. The bonferroni corrected p-value = 3.07125307125307e-05."
# "Number of total bivariate tests: 1628. The FDR corrected p-value = 0.0038654."

FDR_thres = 0.0038654

# 0.05 divided by the number of total bivariate tests performed
bonferroni_thres = 0.05/1628

# chr6:28,477,797-33,448,354 (GRCh37)
hla_start <- 28477797
hla_stop <- 33448354

# Load files -----------------------------------------------------

sign_results_fdr <- read.table(here(project_dir, "lava_results", "onek1k", "FDR_significant_onek1k_bivar_window_100000_20220711.tsv"), sep = "\t", header = T)
sign_results_bonf <- read.table(here(project_dir, "lava_results","onek1k", "bonferroni_significant_onek1k_bivar_window_100000_20220711.tsv"), sep = "\t", header = T)
all_results <- read.table(here(project_dir, "lava_results", "onek1k", "onek1k_bivar_window_100000_20220711.tsv"), sep = "\t", header = T)

sample_sizes <- read.table(here(project_dir, "lava_results", "N_GWAS_traits.txt"), sep = "\t", header = F)
colnames(sample_sizes) <- c("phen1_clean","N")

eQTL_significant <- read.table(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", "significant_genes_per_celltype.txt"), sep = "\t", header = T)

# Main plots -----------------------------------------------------------

### plot 1: bar plot

summary_counts <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean1", NA, "gene"), sep = "_") %>%
  group_by(., phen1_clean) %>%
  count(., phen2_clean1) %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         phen2_clean = case_when(
           phen2_clean1 == "MonoC" ~ "Monocytes",
           phen2_clean1 == "BIN" ~ "Immature and naïve B cells",
           phen2_clean1 == "BMem" ~ "Memory B cells",
           phen2_clean1 == "CD4ET" ~ "Effector memory CD4+ T cells",
           phen2_clean1 == "CD8ET" ~ "Effector memory CD8+ T cells",
           phen2_clean1 == "CD4NC" ~ "Naïve/central memory CD4+ T cells",
           phen2_clean1 == "CD8NC" ~ "Naïve/central memory CD8+ T cells"
         )
  ) %>%
  left_join(., sample_sizes, by = "phen1_clean") %>%
  mutate(phen1_plot = str_c(order_phen1," \n(",N,")"))

# plot (Proportions):
ggplot(data = summary_counts, aes(x = reorder(order_phen1, N), y = n, fill = phen2_clean)) +
  geom_col(colour = "black", position = "fill") +
  theme_bw() +
  scale_fill_manual(values=cbbPalette, name = "QTL cell type (OneK1K)") +
  labs(x = "GWAS trait", y = "Proportion of significant \n local genetic correlations") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 20))

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("proportion_of_significant_local_correlations_onek1k_FDR_",date,".jpg")),
       width = 40, height = 30, units = "cm")

# --------------------------------------------------------------------------------------------------------

### plot: significant correlations including the number of cell types in which the correlation is significant

bivar_to_plot <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA, "gene"), sep = "_") %>%
  mutate(., phen2_legend = case_when(
    phen2_clean == "MonoC" ~ "Monocytes",
    phen2_clean == "BIN" ~ "Immature and naïve B cells",
    phen2_clean == "BMem" ~ "Memory B cells",
    phen2_clean == "CD4ET" ~ "Effector memory CD4+ T cells",
    phen2_clean == "CD8ET" ~ "Effector memory CD8+ T cells",
    phen2_clean == "CD4NC" ~ "Naïve/central memory CD4+ T cells",
    phen2_clean == "CD8NC" ~ "Naïve/central memory CD8+ T cells"
  ),
  fill_rho = case_when(
    p < 0.05 ~ round(rho, 2)
  ),
  chr_gene = str_c("chr ",chr,": ", gene_name),
  cor_pair = str_c(gene_locus,"_",phen1_clean)
  )

# remove HLA region
HLA <- bivar_to_plot %>%
  dplyr::filter(., chr == 6 & start >= hla_start) %>%
  dplyr::filter(., chr == 6 & stop <= hla_stop) %>%
  select(gene_locus)

noHLA = setdiff(bivar_to_plot$gene_locus, HLA$gene_locus)

bivar_noHLA <- bivar_to_plot %>%
  filter(., gene_locus %in% noHLA)

n_celltypes <- count(bivar_noHLA, cor_pair)

bivar_to_plot <- bivar_noHLA %>%
  left_join(., n_celltypes, by = "cor_pair")

suppfig2B <- ggplot(bivar_to_plot, aes(x=phen1_clean, y=reorder(gene_name, start))) + 
  geom_count(colour = "black", shape = 15, size = 12) +
  geom_count(aes(colour = as.factor(n)), shape = 15, size = 11) +
  coord_fixed() +
  theme_bw() +
  scale_colour_brewer() +
  labs(x = "Traits", y = "Expressed Genes", title = "") +
  guides(colour = guide_colourbar(title = "Local rg")) +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 0.8, vjust = 0),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24)) +
  geom_text(aes(label=n), size = 8, vjust = 0.5) +
  coord_flip() +
  facet_wrap(~chr, scales = "free", nrow = 3)

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("all_FDR_significant_shared_across_cell_types_onek1k_noHLA",date,".jpg")),
       width = 110, height = 50, units = "cm")
# --------------------------------------------------------------------------------------------------------

### plot: number of significantly expressed genes 

# plot code:
ggplot(data = eQTL_significant, aes(x = reorder(cell_type, significant_genes), y = significant_genes)) +
  geom_col(fill = "#CCCCCC", colour = "black") +
  theme_bw() +
  labs(x = "Cell type", y = "Number of significantly \n expressed genes") +
  theme(axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 20)) + 
  geom_text(aes(label = significant_genes, vjust = 2), size= 8)

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("summary_of_significant_QTL_genes_expressed_",date,".jpg")),
       width = 30, height = 30, units = "cm")
# ----------------------------------------------------------------------------------------------------

### plot: significant correlations for a specific gene across cell types
# plot modified for paper including BIN1, RAB7L1 and SCIMP

bivar_to_plot <-
  all_results %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA, "gene"), sep = "_") %>%
  mutate(., phen2_legend = case_when(
    phen2_clean == "MonoC" ~ "Monocytes",
    phen2_clean == "BIN" ~ "Immature and naïve B cells",
    phen2_clean == "BMem" ~ "Memory B cells",
    phen2_clean == "CD4ET" ~ "Effector memory CD4+ T cells",
    phen2_clean == "CD8ET" ~ "Effector memory CD8+ T cells",
    phen2_clean == "CD4NC" ~ "Naïve/central memory CD4+ T cells",
    phen2_clean == "CD8NC" ~ "Naïve/central memory CD8+ T cells"),
    fill_rho = case_when(
      p < 0.05 ~ round(rho, 2)),
    chr_gene = str_c("chr ",chr,": ", gene_name),
    cor_pair = str_c(gene_locus,"_",phen1_clean),
    sig = case_when(
      p < FDR_thres ~ "**"),
    nom_sig = case_when(
      p < 0.05 & p >= FDR_thres ~ "."),
    not_sig = case_when(
      p >= 0.05 ~ "ns")
  )

# remove HLA region
HLA <- bivar_to_plot %>%
  dplyr::filter(., chr == 6 & start >= hla_start) %>%
  dplyr::filter(., chr == 6 & stop <= hla_stop) %>%
  select(gene_locus)

noHLA = setdiff(bivar_to_plot$gene_locus, HLA$gene_locus)

bivar_noHLA <- bivar_to_plot %>%
  filter(., gene_locus %in% noHLA)

n_celltypes <- count(bivar_noHLA, cor_pair)

bivar_to_plot <- bivar_noHLA %>%
  left_join(., n_celltypes, by = "cor_pair")

# plot only one/subset of genes:
bivar_to_plot <- bivar_to_plot %>% filter(., chr == 2 | chr == 17 | chr == 1) %>% 
  filter(., gene_name == "BIN1" | gene_name == "SCIMP" | gene_name == "RAB7L1")

ggplot(bivar_to_plot, aes(x=phen1_clean, y=phen2_clean)) + 
  geom_count(colour = "black", shape = 15, size = 12) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 11) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(palette = "RdYlBu", limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "Traits", y = "Cell types", title = "") +
  guides(colour = guide_colourbar(title = "Local rg")) +
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 0.8, vjust = 0),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_text(size = 24)) +
  geom_text(aes(label=sig), size = 8, vjust = 0.75) +
  geom_text(aes(label=nom_sig), size = 8, vjust = 0.1) +
  geom_text(aes(label=not_sig), size = 5, vjust = 0.4) +
  coord_flip() +
  facet_wrap(~gene_name, scales = "free", nrow = 1)

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("ORMDL3_FDR_significant_shared_across_cell_types_onek1k",date,".jpg")),
       width = 12, height = 10, units = "cm")
# ----------------------------------------------------------------------------------------------------

# plot: bar plot of chromosomes stratified by disease

# palettes:
fishy = c("#6388b4", "#ffae34", "#ef6f6a", "#8cc2ca", "#c3bc3f", "#55ad89", "#bb7693", "#baa094", "#767676")
fishy2 =c("#6388b4", "#ffae34", "#8cc2ca", "#55ad89", "#bb7693", "#baa094", "#767676")

# remove HLA region
HLA <- sign_results_fdr %>%
  dplyr::filter(., chr == 6 & start >= hla_start) %>%
  dplyr::filter(., chr == 6 & stop <= hla_stop) %>%
  select(gene_locus)

noHLA = setdiff(sign_results_fdr$gene_locus, HLA$gene_locus)

bivar_noHLA <- sign_results_fdr %>%
  filter(., gene_locus %in% noHLA)

summary_counts_bis <-
  bivar_noHLA %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean1", NA, "gene"), sep = "_") %>%
  group_by(., chr) %>%
  count(., phen1_clean)

# add chr 9 with zeros for visual purposes:
chr9 <- data.frame("chr" = rep(9, 7), "phen1_clean" = c("AD","PD","ALS","SCZ","MS","UC","CD"), "n" = rep(0, 7))

summary_counts_bis <- rbind(summary_counts_bis, chr9) %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD"))
  )

suppfig2A <- ggplot(data = summary_counts_bis, aes(x = reorder(as.character(chr), chr), y = n, fill = order_phen1)) +
  geom_col(colour = "black") +  
  theme_bw() +
  scale_fill_manual(values=fishy2, name = "Traits") +
  labs(x = "Chromosomes", y = "Number of significant \n local genetic correlations") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 20)) 

suppfig2A
ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("barplot_signals_per_chromosome",date,".jpg")),
       width = 30, height = 20, units = "cm")

###### merge two plots in one figure for supplementary figure of paper:
ggarrange(suppfig2A, suppfig2B,
          labels = c("A", "B"),
          ncol = 1, nrow = 2, font.label = 30)

ggsave(here(project_dir, "lava_plots", "lava_QTLs", "onek1k", stringr::str_c("Supp_Figure_2_paper_",date,".jpg")),
       width = 120, height = 80, units = "cm")

