# Description: analysis of shared significant genes across cell types

# Load packages -----------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

list_eqtls <- list.files(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","QTLs"), 
                         recursive = T, full.names = F, pattern = "*.lava.gz") %>%
  gsub(".lava.gz","",.) %>%
  gsub("_1MscBloodNL","",.) %>%
  #str_remove(., ".lava.gz") %>%
  #str_remove(., "_1MscBloodNL") %>%
  as.data.frame() %>%
  rename(., file_name = .) %>%
  separate(file_name, c("cell_type","gene"), sep = "_")

# Main --------------------------------------------------------------------

# counts all
count_per_celltype <- count(list_eqtls, cell_type)
count_per_gene <- count(list_eqtls, gene) %>%
  count(., n, name = "freq")

# count only 4 cell types
count_per_celltype_filt <- list_eqtls %>%
  filter(., cell_type == "Bcells" | cell_type == "CD4Tcells" | cell_type == "CD8Tcells" | cell_type == "monocyte") %>%
  count(., cell_type)

count_per_gene_filt <- list_eqtls %>%
  filter(., cell_type == "Bcells" | cell_type == "CD4Tcells" | cell_type == "CD8Tcells" | cell_type == "monocyte") %>%
  count(., gene) %>%
  count(., n, name = "freq")

# Plots --------------------------------------------------------------------

# histogram of counts per gene
ggplot(count_per_gene, aes(x = n, y = freq)) + 
  geom_col(fill = "darkolivegreen3", colour = "black") +
  theme_classic() +
  labs(x = "Number of cell types across expressed genes", y = "Frequency") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7)) +
  scale_y_continuous(breaks = seq(0,250,50)) +
  theme(axis.title.x = element_text(face="bold", size=26),
        axis.text.x  = element_text(size=22),
        axis.title.y = element_text(face="bold", size=26),
        axis.text.y  = element_text(size=22)) +
  geom_text(aes(label = freq), vjust = 1.5, size = 12, colour = "black")
ggsave(here(project_dir,"test_loci","figures","counts_per_gene_expressed_across_cell_types.pdf"),
       width = 40, height = 30, units = "cm")

# histogram of counts per gene filtered
ggplot(count_per_gene_filt, aes(x = n, y = freq)) + 
  geom_col(fill = "darkolivegreen3", colour = "black") +
  theme_classic() +
  labs(x = "Number of cell types across expressed genes", y = "Frequency") +
  scale_x_continuous(breaks=c(1,2,3,4)) +
  scale_y_continuous(breaks = seq(0,250,50)) +
  theme(axis.title.x = element_text(face="bold", size=26),
        axis.text.x  = element_text(size=22),
        axis.title.y = element_text(face="bold", size=26),
        axis.text.y  = element_text(size=22)) +
  geom_text(aes(label = freq), vjust = 1.5, size = 12, colour = "black")
ggsave(here(project_dir,"test_loci","figures","counts_per_gene_expressed_across_cell_types_filtered.pdf"),
       width = 40, height = 30, units = "cm")

# histogram of counts per cell type
ggplot(count_per_celltype, aes(x = cell_type, y = n)) + 
  geom_col(fill = "darkmagenta", colour = "black") +
  theme_classic() +
  labs(x = "Cell types", y = "Number of genes") +
  scale_y_continuous(breaks = seq(0,200,50)) +
  theme(axis.title.x = element_text(face="bold", size=26),
        axis.text.x  = element_text(size=22),
        axis.title.y = element_text(face="bold", size=26),
        axis.text.y  = element_text(size=22)) +
  geom_text(aes(label = n), vjust = 1.5, size = 12, colour = "white")
ggsave(here(project_dir,"test_loci","figures","counts_of_genes_expressed_per_cell_type.pdf"),
       width = 40, height = 30, units = "cm")

# histogram of counts per cell type filtered
ggplot(count_per_celltype_filt, aes(x = cell_type, y = n)) + 
  geom_col(fill = "darkmagenta", colour = "black") +
  theme_classic() +
  labs(x = "Cell types", y = "Number of genes") +
  scale_y_continuous(breaks = seq(0,200,50)) +
  theme(axis.title.x = element_text(face="bold", size=26),
        axis.text.x  = element_text(size=22),
        axis.title.y = element_text(face="bold", size=26),
        axis.text.y  = element_text(size=22)) +
  geom_text(aes(label = n), vjust = 1.5, size = 12, colour = "white")
ggsave(here(project_dir,"test_loci","figures","counts_of_genes_expressed_per_cell_type_filtered.pdf"),
       width = 40, height = 30, units = "cm")

