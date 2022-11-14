# Description: formatting of pQTLs summary statistics for LAVA.

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(devtools)
library(stringr)
library(BiocManager)
library(dplyr)
library(data.table)
library(biomaRt)
library(BSgenome)
library(colochelpR)
library(rutils)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"
path_hg38_hg37 <- here(project_dir,"reference_data","hg38ToHg19.over.chain")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                             pQTLs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

files <- setNames(list.files(path = here(project_dir,"GWAS_summary_statistics","pQTLs_Zhang2022","EA"), 
                             full.names = T, pattern = "*PHENO1.glm.linear"),
                  nm = list.files(path = here(project_dir,"GWAS_summary_statistics","pQTLs_Zhang2022","EA"), 
                                  full.names = F, pattern = "*PHENO1.glm.linear") %>%
                    stringr::str_remove_all(., ".PHENO1.glm.linear"))

seq_ids <- read.table(here(project_dir,"GWAS_summary_statistics","pQTLs_Zhang2022", "seqid.txt"),
                      sep = "\t", header = T) %>%
  mutate(gene_name = stringr::str_extract(entrezgenesymbol, "[:alnum:]+"))

pqtls_all <- lapply(files, function(x) fread(x, sep = "\t", header= T) %>%
  mutate(P = as.numeric(P))) %>%
  dplyr::bind_rows(., .id = "seqid_in_sample") %>%
  left_join(., seq_ids, by = "seqid_in_sample")

significant_genes <- pqtls_all %>%
  dplyr::filter(P < 5e-08) %>%
  dplyr::select(gene_name) %>%
  distinct(., gene_name)

cat("There are",nrow(significant_genes),"significant genes across all pQTLs.\n")

pqtls_significant <- pqtls_all %>%
  dplyr::mutate(CHR = as.factor(CHROM), N = 7213,
                A2 = case_when(A1 != REF ~ REF, TRUE ~ ALT)) %>%
  dplyr::filter(gene_name %in% significant_genes$gene_name) %>%
  dplyr::select(
    SNP = ID,
    GENE = gene_name,
    CHR,
    BP = POS,
    A1,
    A2,
    BETA,
    SE,
    P,
    N
  ) %>%
  rutils::liftover_coord(df = ., path_to_chain = path_hg38_hg37)

pqtls_list <- setNames(pqtls_significant %>%
                         dplyr::group_split(GENE),
                       nm = str_c("pqtls_",
                                  pqtls_significant$GENE %>%
                                    unique() %>%
                                    sort()
                       )
)

for(i in 1:length(pqtls_list)){
  
  fwrite(
    pqtls_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","pqtls",
           stringr::str_c(names(pqtls_list)[i],lava_ext)
      ),
    sep = "\t"
  )
  
}

