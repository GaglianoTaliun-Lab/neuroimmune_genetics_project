# Description: load and wrangle Onek1k beta summary statistics for coloc per gene/celltype/region

# Load packages -----------------------------------------------------------

library(here)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(rutils)
library(R.utils)
library(utils)
library(MafDb.1Kgenomes.phase3.hs37d5)
library(GenomicScores)
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(colochelpR)

# Arguments ---------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

mafdb <- MafDb.1Kgenomes.phase3.hs37d5
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# Main ---------------------------------------------------------------

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc_FDR.txt"),
                         sep = "\t", header = T)

# chr	gene	trait	cell_type	start	end	ensembl_id
# 1	RP11_108M9_4	PD	CD8ET	17115033	17316144	ENSG00000238142
# 1	AK5	PD	CD8ET	77647736	78125651	ENSG00000154027

mafs <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_eur.frq")) %>% dplyr::select(SNP, MAF)
# colnames(mafs) <- c("CHR","SNP","A1","A2","MAF","NCHROBS")

for (row in 1:nrow(gene_table)) {
  
  cat("Start processing row # ", row, ".\n")
  
  curr_chr = gene_table[row,1]
  curr_gene = gene_table[row,2]
  curr_celltype = gene_table[row,4]
  curr_start = gene_table[row,5]
  curr_end = gene_table[row,6]
  curr_ensemblid = gene_table[row,7]
  
  curr_gene2 <- str_replace_all(curr_gene, "_", "-")
  
  eQTL_name = str_c(curr_celltype, "_onek1k_", curr_ensemblid)
  
  # read a file with chosen chromosome and specific cell type (from gene_table) 
  # onek1k columns: SNP	gene	beta	t-stat	p-value	FDR
  
  chr_file <- read.table(here(project_dir, "GWAS_summary_statistics", "onek1k", "OneK1K_matrix_eQTL_results", str_c(curr_celltype,"_chr",curr_chr,"_cis_eqtls_210211.tsv")), sep = "\t", header = TRUE) %>%
    filter(., gene == curr_gene2)
  head(chr_file)
  
  colnames(chr_file) <- c("SNP","gene","beta","tstat","pvalues","FDR")
  head(chr_file)
  
  if (nrow(chr_file) > 0) {
    
    chr_file <- chr_file %>%
      filter(., SNP %like% "^rs") %>%
      colochelpR::convert_rs_to_loc(df = ., SNP_column = "SNP", dbSNP = dbsnp_144) %>%
      separate(loc, c("CHR", "BP"), sep = ":") %>%
      filter(., BP >= curr_start & BP <= curr_end) %>%
      mutate(gene = curr_gene, eQTL_dataset = eQTL_name, N = 982, se = beta/tstat) %>%
      colochelpR::get_varbeta(.) %>%
      dplyr::select(., eQTL_dataset, 
                    gene,
                    SNP,
                    CHR,
                    BP,
                    BETA = beta,
                    varbeta,
                    SE = se,
                    pvalues,
                    N)
  }
  
    coloc_out <- chr_file %>%
      inner_join(., mafs, by = "SNP")
    
    out_name=here(project_dir, "colocalization", "QTL_regional_sumstats", str_c("coloc_onek1k_",curr_celltype,"_",curr_gene,".tsv"))
    fwrite(coloc_out, file = out_name, sep = "\t")
    
    cat("---Finished processing row # ", row, ".\n")
  
}

