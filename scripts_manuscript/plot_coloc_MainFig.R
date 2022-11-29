# Description: summary plot for COLOC results

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(cowplot)
library(gridExtra)
library(grid)
library(ggraph)
library(RColorBrewer)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/neuro_immune"

date=str_remove_all(Sys.Date(), "-")

# okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbbPalette <- c("purple", "#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette_man <- c("purple", "#E69F00")

coloc_susie <- list()

# chr6:28,477,797-33,448,354 (GRCh37)
hla_start <- 28477797
hla_stop <- 33448354

# Load files -----------------------------------------------------

coloc_res_H4 <- read.table(here(project_dir, "colocalization", "coloc_results_FDR", "summary_H4_0.8", "summary_of_results_H4_0.8.tsv"), sep = "\t", header = T) %>%
  rename(PP = H4) %>%
  mutate(PP_type = "H4")
coloc_res_H3 <- read.table(here(project_dir, "colocalization", "coloc_results_FDR", "summary_H3_0.8", "summary_of_results_H3_0.8.tsv"), sep = "\t", header = T) %>%
  rename(PP = H3) %>%
  mutate(PP_type = "H3")
coloc_res <- rbind(coloc_res_H4, coloc_res_H3)

coloc_susie_files <- list.files(here(project_dir, "colocalization", "coloc_results_FDR", "coloc_with_susie", "summary_all"), full.names = T)

coloc_susie_names <- list.files(here(project_dir, "colocalization", "coloc_results_FDR", "coloc_with_susie", "summary_all"), full.names = F) %>%
  stringr::str_remove(., ".tsv") %>%
  stringr::str_remove(., "coloc.susie_summary_")

for (i in 1:length(coloc_susie_files)) {
  
  coloc_susie[[i]] <- read.table(coloc_susie_files[i], sep = "\t", header = F)
  
}

for (i in 1:length(coloc_susie_names)){
  
  coloc_susie[[i]] <- coloc_susie[[i]] %>% 
    mutate(coloc_name = coloc_susie_names[i])
  
}

names(coloc_susie) <- coloc_susie_names
coloc_susie_res <- data.table::rbindlist(coloc_susie) %>%
  dplyr::select(coloc_name,
                n_snps = V2,
                hit1 = V3,
                hit2 = V4,
                H0 = V5,
                H1 = V6,
                H2 = V7,
                H3 = V8,
                H4 = V9,
                idx1 = V10,
                idx2 = V11) %>%
  separate(., coloc_name,  c("trait", "cell_type", NA, "ensembl_id"), sep = "_")

# Plot -----------------------------------------------------------

### Figure 5: count plot per cell type for FDR significant results (H3 or H4 >= 0.8) - no HLA

coloc_plot <-
  coloc_res %>%
  filter(., PP >= 0.8) %>%
  filter(., PP_type == "H4" | PP_type == "H3" | PP_type == "H2") %>%
  mutate(., order_phen1 = factor(trait, levels = c("AD","PD","LBD","ALS","SCZ","MS","UC","CD")),
         cell_type_legend = case_when(
           cell_type == "MonoC" ~ "Monocytes",
           cell_type == "BIN" ~ "Immature and naïve B cells",
           cell_type == "BMem" ~ "Memory B cells",
           cell_type == "CD4ET" ~ "Effector memory CD4+ T cells",
           cell_type == "CD8ET" ~ "Effector memory CD8+ T cells",
           cell_type == "CD4NC" ~ "Naïve/central memory CD4+ T cells",
           cell_type == "CD8NC" ~ "Naïve/central memory CD8+ T cells"
         ),
         fill_PP = case_when(
           PP > 0.8 ~ round(PP,1)
         ),
         chr_gene = str_c(chr, ": ", gene),
         cell_type_trait = str_c(ensembl_id, "_", trait),
         PP_type_desc = case_when(
           PP_type == "H3" ~ "distinct causal variants",
           PP_type == "H4" ~ "same causal variant",
           PP_type == "H2" ~ "no signal for one trait"
         )
  )

# remove HLA genes
HLA <- coloc_plot %>%
  dplyr::filter(., chr == 6 & start >= hla_start) %>%
  dplyr::filter(., chr == 6 & end <= hla_stop) %>%
  select(ensembl_id)
noHLA = setdiff(coloc_plot$ensembl_id, HLA$ensembl_id)
coloc_noHLA <- coloc_plot %>%
  filter(., ensembl_id %in% noHLA)

# reorder based on chr, start position:
coloc_noHLA <- coloc_noHLA[with(coloc_noHLA, order(chr, start)),]
coloc_noHLA$chr_gene <- factor(coloc_noHLA$chr_gene, levels = unique(coloc_noHLA$chr_gene))

# number of cell types for significant colocalization
n_celltypes <- coloc_noHLA %>%
  filter(., PP_type == "H4") %>%
  count(., cell_type_trait)

coloc_noHLA <- coloc_noHLA %>%
  left_join(., n_celltypes, by = "cell_type_trait") %>%
  mutate(n_H4 = case_when(
    PP_type == "H4" ~ n
  ))

# plot:
ggplot(coloc_noHLA, aes(x=order_phen1, y=chr_gene, chr)) + 
  geom_count(colour = "black", shape = 16, size = 12) +
  geom_count(aes(fill = PP_type_desc), shape = 21, size = 11) +
  coord_fixed() +
  theme_bw() +
  scale_fill_manual(values=palette_man, name = "PP >= 0.8") +
  labs(x = "GWAS Trait", y = "Expressed Genes", title = "") +
  # guides(colour = guide_colourbar(title = "Cell type")) +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 0.95, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size =22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20)) +
  # geom_text(aes(label=fill_PP), size = 4, vjust = 0.5) +
  geom_text(aes(label=n_H4), size = 6, vjust = 0.5) +
  coord_flip() #+
#facet_wrap(~chr, scale = "free", nrow = 4)

ggsave(here(project_dir, "colocalization", "plots", stringr::str_c("FDR_coloc_H4_H3_noHLA_nocells_",date,".jpg")),
       width = 50, height = 20, units = "cm")
