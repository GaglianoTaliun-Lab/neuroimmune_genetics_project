# Description: summary plots for LAVA results with QTLs

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(cowplot)
library(gridExtra)
library(grid)

# Arguments ------------------------------------------------------

project_dir = "Documents/research-projects/neuro_immune"

date=str_remove_all(Sys.Date(), "-")

cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

FDR_thres = 0.00437756

# 0.05 divided by the number of total bivariate tests performed
bonferroni_thres = 0.05/358

# Load files -----------------------------------------------------

sign_results_fdr <- read.table(here(project_dir, "lava_results", "FDR_significant_results_bivar_window_100000_20220425.tsv"), sep = "\t", header = T)
sign_results_bonf <- read.table(here(project_dir, "lava_results", "bonferroni_significant_results_bivar_window_100000_20220425.tsv"), sep = "\t", header = T)
all_results <- read.table(here(project_dir, "lava_results", "all_results_bivar_window_100000_20220425.tsv"), sep = "\t", header = T)

# Main plots -----------------------------------------------------------

### plot 1: bar plot

# wrangle data for plot:
# - remove author,year from phen1 and "_1MscBloodNL" from phen2
# - group by phen1 (GWAS) and count incidences per phen2 (QTL):

summary_counts <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA, "gene"), sep = "_") %>%
  group_by(., phen1_clean) %>%
  count(., phen2_clean) %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         phen2_clean = case_when(
           phen2_clean == "monocyte" ~ "Monocytes",
           phen2_clean == "Bcells" ~ "B cells",
           phen2_clean == "CD4Tcells" ~ "CD4+ T cells",
           phen2_clean == "CD8Tcells" ~ "CD8+ T cells"
         )
  )

# plot code:
ggplot(data = summary_counts, aes(x = order_phen1, y = n, fill = phen2_clean)) +
  geom_col(colour = "black") +
  theme_bw() +
  scale_fill_manual(values=cbbPalette, name = "QTL cell type") +
  labs(x = "Trait", y = "Number of significant \n local genetic correlations") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(2,40,2)) +
  geom_text(aes(label = "FDR < 0.01", x = "CD", y = 40))

# Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", stringr::str_c("summary_of_significant_local_correlations_FDR_",date,"_.pdf")),
       width = 40, height = 30, units = "cm")
# --------------------------------------------------------------------------------------------------------

### plot 2: count plot per cell type for significant results

# wrangle data for plot:

# - count genes grouped by cell type
bivar_to_plot <-
  sign_results_fdr %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA, "gene"), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         phen2_clean = case_when(
           phen2_clean == "monocyte" ~ "Monocytes",
           phen2_clean == "Bcells" ~ "B cells",
           phen2_clean == "CD4Tcells" ~ "CD4+ T cells",
           phen2_clean == "CD8Tcells" ~ "CD8+ T cells"
         )
  ) %>%
  mutate(chr_gene = str_c("chr ",chr,": ", gene_name))

# - plot code:
plots = list()
for (i in 1:nrow(count(bivar_to_plot,phen2_clean))){
  
  bivar_per_celltype <- filter(bivar_to_plot, phen2_clean == count(bivar_to_plot,phen2_clean) %>% pull(1) %>% sort() %>% .[i])
  
  plots[[i]] <- ggplot(bivar_per_celltype, aes(x=phen1_clean, y=reorder(chr_gene, chr))) + 
    geom_count(colour = "black", shape = 15, size = 8) +
    geom_count(aes(colour = rho), shape = 15, size = 7) +
    coord_fixed() +
    theme_minimal() +
    scale_colour_distiller(
      palette = "RdYlBu",
      limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
    labs(x = "Trait", y = "Expressed Genes") +
    guides(colour = guide_colourbar(title = "Local rg")) +
    theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size =16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    coord_flip()
  
}

# - display 4 figures in one plot
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
          ncol = 2, 
           labels = unique(bivar_to_plot$phen2_clean) %>% sort(), 
           rel_heights = c(2,2), rel_widths = c(2.5,3))

# - Save plot --------------------------------------------------------------------------------------------
ggsave(here(project_dir, "lava_plots", "lava_QTLs", stringr::str_c("local_genetic_correlations_per_cell_type_",date,".pdf")),
       width = 80, height = 35, units = "cm")
# --------------------------------------------------------------------------------------------------------

### plot 3: count plot per cell type for all results

# - count genes grouped by cell type
bivar_to_plot <-
  all_results %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA, "gene"), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         phen2_clean = case_when(
           phen2_clean == "monocyte" ~ "Monocytes",
           phen2_clean == "Bcells" ~ "B cells",
           phen2_clean == "CD4Tcells" ~ "CD4+ T cells",
           phen2_clean == "CD8Tcells" ~ "CD8+ T cells"
         ),
         fill_rho = case_when(
           p < 0.05 ~ round(rho, 2)
         ),
         sig = case_when(
           p < FDR_thres ~ "**"
         ),
         nom_sig = case_when(
           p < 0.05 & p >= FDR_thres ~ ".",
         ),
         chr_gene = str_c("chr ",chr,": ", gene_name)
         )

# - create and save one plot per cell type:
for (i in 1:nrow(count(bivar_to_plot,phen2_clean))){
  
  cell_type = count(bivar_to_plot,phen2_clean) %>% pull(1) %>% sort() %>% .[i]
  
  bivar_per_celltype <- filter(bivar_to_plot, phen2_clean == cell_type)
  
  ggplot(bivar_per_celltype, aes(x=phen1_clean, y=reorder(chr_gene, chr))) + 
    geom_count(colour = "black", shape = 15, size = 8) +
    geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
    coord_fixed() +
    theme_minimal() +
    scale_colour_distiller(
      palette = "RdYlBu",
      limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
    labs(x = "Trait", y = "Expressed Genes", title = stringr::str_c("Cell type: ", cell_type)) +
    guides(colour = guide_colourbar(title = "Local rg")) +
    theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size =16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    geom_text(aes(label=sig), size = 4, vjust = 0.75) +
    geom_text(aes(label=nom_sig), size = 5, vjust = 0.1) +
    coord_flip()
  
  ggsave(here(project_dir, "lava_plots", "lava_QTLs", stringr::str_c("all_local_genetic_correlations_",cell_type,"_",date,".pdf")),
         width = 60, height = 25, units = "cm")
  
}
