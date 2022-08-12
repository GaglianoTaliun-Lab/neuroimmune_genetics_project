# Description: compare cell types effect sizes from OneK1K data

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

celltype1 = "CD8NC"
celltype2 = "CD8ET"

# Load files -----------------------------------------------------

cell1_files <- list.files(path = here(project_dir,"GWAS_summary_statistics","onek1k","OneK1K_matrix_eQTL_results"),
                          pattern = str_c(celltype1,"_chr*"), full.names = T, all.files = T)
cell2_files <- list.files(path = here(project_dir,"GWAS_summary_statistics","onek1k","OneK1K_matrix_eQTL_results"),
                          pattern = str_c(celltype2,"_chr*"), full.names = T, all.files = T) 

# Read files and prepare for plot --------------------------------

cell1 <- lapply(cell1_files, function(x) read.table(x, sep = "\t", header  = T)) %>%
  data.table::rbindlist(.) %>%
  # mutate(SE_cell1 = beta/t.stat) %>%
  rename(beta_cell1 = beta, pvalue_cell1 = p.value)

cell2 <- lapply(cell2_files, function(x) read.table(x, sep = "\t", header  = T)) %>%
  data.table::rbindlist(.) %>%
  rename(beta_cell2 = beta, pvalue_cell2 = p.value)

cell1_genes <- cell1 %>%
  filter(., pvalue_cell1 < 5e-08) %>%
  dplyr::select(gene) %>%
  distinct(., gene)

cell2_genes <- cell2 %>%
  filter(., pvalue_cell2 < 5e-08) %>%
  dplyr::select(gene) %>%
  distinct(., gene)

unique_genes <- rbind(cell1_genes, cell2_genes) %>%
  distinct(., gene)

cells_plot <- inner_join(cell1, cell2, by = c("gene", "SNP")) %>%
# filter(., pvalue_cell1 < 5e-08 | pvalue_cell2 < 5e-08)
  filter(., gene %in% unique_genes$gene) %>%
  mutate(sign = case_when(
    pvalue_cell1 < 5e-08 & pvalue_cell2 < 5e-08 ~ "both significant",
    pvalue_cell1 < 5e-08 & pvalue_cell2 >= 5e-08 ~ "significant cell 1",
    pvalue_cell1 >= 5e-08 & pvalue_cell2 < 5e-08 ~ "significant cell 2",
    pvalue_cell1 >= 5e-08 & pvalue_cell2 >= 5e-08 ~ "not significant"
    )
  )
 
cat("Rows present in merged file: ", nrow(cells_plot), ".\n")

cells_plot_ss <- sample_n(cells_plot, 1000000)

head(cells_plot_ss)

# Plot -----------------------------------------------------------

ggplot(data = cells_plot_ss, aes(x = beta_cell1, y = beta_cell2, colour = sign)) +
  geom_point() + 
  # geom_errorbar(aes(ymin = beta_cell2 - SE_cell2, ymax = beta_cell2 + SE_cell2)) +
  theme_minimal()
ggsave(here(project_dir, "GWAS_summary_statistics", "onek1k", str_c("cell_comparison_plot_",celltype1,"_",celltype2,".pdf")),
       width = 40, height = 40, units = "cm")
