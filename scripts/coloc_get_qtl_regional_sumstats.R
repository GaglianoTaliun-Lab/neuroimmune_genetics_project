# create regional files for onek1k sc-eQTLs to be used for coloc

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(rutils)
library(R.utils)
library(MafDb.1Kgenomes.phase3.hs37d5)
library(GenomicScores)

# Arguments ------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

mafdb <- MafDb.1Kgenomes.phase3.hs37d5

# Load files -----------------------------------------------------------

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc_FDR.txt"),
                         sep = "\t", header = T)

lava_qtl_files <- list.files(here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "onek1k"), 
                             pattern = ".lava.gz", full.names = T) %>% as.data.frame()
colnames(lava_qtl_files) <- "lava_path"

mafs <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_eur.frq")) %>% dplyr::select(SNP, MAF)
# colnames(mafs) <- c("CHR","SNP","A1","A2","MAF","NCHROBS")

head(mafs)

# Format dataframes and save regional sumstats -----------------------------------------------------------

lava_qtl_files <- lava_qtl_files %>%
  mutate(lava_filename = basename(lava_path))

head(lava_qtl_files)

gene_table <- gene_table %>%
  mutate(lava_filename = stringr::str_c(gene_table$cell_type, "_onek1k_", gene_table$ensembl_id, ".lava.gz")) %>%
  left_join(., lava_qtl_files, by = "lava_filename")

nrow(gene_table)
gene_table

for (i in 1:nrow(gene_table)) {
  
    lava_qtl_sumstats <- fread(gene_table$lava_path[i])
      
      chr_pos = gene_table$chr[i]
      start_pos = gene_table$start[i]
      end_pos = gene_table$end[i]
      
        lava_sumstats_region <- filter(lava_qtl_sumstats, CHR == chr_pos & BP >= start_pos & BP <= end_pos) %>%
          mutate(eQTL_dataset = str_c(gene_table$cell_type[i], "_onek1k_", gene_table$ensembl_id[i]),
                 gene = gene_table$ensembl_id[i]) %>%
          dplyr::select(., eQTL_dataset, gene, SNP, CHR, BP, beta = BETA, se = SE, pvalues = P, A1, A2, N) %>%
          colochelpR::get_varbeta(.)
    
#      mafs <- GenomicScores::gscores(x = mafdb, ranges = unique(lava_sumstats_region$SNP) %>% as.character(), pop = "EUR_AF")
#      mafs <- mafs %>%
#        as.data.frame() %>%
#        tibble::rownames_to_column(var = "SNP") %>%
#        dplyr::rename(MAF = EUR_AF) %>%
#        dplyr::select(SNP, MAF)
      
      lava_sumstats_region <- lava_sumstats_region %>%
        inner_join(mafs, by = "SNP")

      write.table(lava_sumstats_region, here(project_dir, "colocalization", "QTL_regional_sumstats",
                                             stringr::str_c("onek1k_", gene_table$cell_type[i], "_", gene_table$gene[i], ".tsv")),
                  sep = "\t", row.names = F, quote = F)

}

