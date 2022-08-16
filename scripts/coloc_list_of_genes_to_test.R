# create a table with the genic regions to test for coloc
# to be used to create regional files for GWAS and QTL summary statistics

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# Load files -----------------------------------------------------------

test_list <- read.table(here(project_dir, "colocalization", "test_list_coloc.txt"),
                        sep = "\t", header = T)

gtf_file <- read.table(here(project_dir, "reference_data", "Homo_sapiens.GRCh37.87.gtf"),
                       sep = "\t", header = F, skip = 5) %>%
  .[,c(1:5,9)]

# Format GTF file ------------------------------------------------------

colnames(gtf_file) = c("chr","source","type","start","end","info")

gtf_file <- gtf_file %>%
  mutate(
    ensembl_id = stringr::str_extract(info, "ENS[:alpha:][:alnum:]+"),
    gene_name = stringr::str_extract(info, "gene_name [:graph:]+") %>%
      stringr::str_remove_all(., "gene_name ") %>%
      stringr::str_remove_all(., ";") %>%
      stringr::str_remove_all(., '"')
  ) %>%
    select(.,
           source, type, start, end, ensembl_id, gene_name)

# Merge both dataframes ------------------------------------------------

gene_table <- left_join(test_list, gtf_file, by = c("gene" = "gene_name")) %>%
filter(., type == "gene") %>%
filter(., source %in% c("ensembl","ensembl_havana"))

# Save new table -------------------------------------------------------

write.table(gene_table, here(project_dir, "colocalization", "gene_table_coloc.txt"),
            sep = "\t", row.names = F, quote = F)
