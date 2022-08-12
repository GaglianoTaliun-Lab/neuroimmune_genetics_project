# Description: obtain the significant genes from the sc-eQTL results (with FDR correction) to filter the locus file used for the bulk QTL LAVA analysis.

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments -----------------------------------------------------

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# Read files ----------------------------------------------------

locus_file <- read.table(here(project_dir, "test_loci", "window_100000", "gene_filtered.loci"), sep = "\t", header = T) %>%
  rename(gene_locus = LOC)

fdr_results <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_results_bivar_window_100000_20220521.tsv"), sep = "\t", header = T)

# Join -----------------------------------------------------------

locus_filtered <- right_join(locus_file, fdr_results, by = "gene_locus") %>%
  select(LOC = gene_locus,
         CHR, 
         START,
         STOP) %>%
  unique(.)

# save new locus file --------------------------------------------

write.table(locus_filtered, here(project_dir, "test_loci", "window_100000", "gene_filtered_bulkQTL.loci"), row.names = F, quote = F, sep = "\t")
