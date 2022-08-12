# Description: preprocess onek1k sc-eQTL summary statistics into LAVA format, 
# to output one file per gene and per cell type

# Load libraries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(devtools)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(here)
library(purrr)
library(glue)
library(utils)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
# project_dir = here("Documents/research-projects/neuro_immune")
lava_ext = ".lava.gz"

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read file downloaded an unzipped from https://onek1k.org/
onek1k <- read_tsv(unzip(here(project_dir,"GWAS_summary_statistics","onek1k","onek1k_eqtl_dataset.zip")))
# onek1k <- read.table(here(project_dir, "GWAS_summary_statistics/onek1k/onek1k_eqtl_dataset.tsv"), sep = "\t", header = T)
# onek1k <- read.table(here(project_dir, "test.tsv"), sep = "\t", header = T)

# create list of tibbles by cell type
onek1k_list <- setNames(onek1k %>%
                           dplyr::group_split(CELL_ID) %>%
                          as.list(),
                         nm = str_c(onek1k$CELL_ID %>%
                                      unique() %>%
                                      sort(),
                                    "_onek1k"
                         )
)

# extract unique eGenes with significant eQTLs for each cell type
unique_sign_genes_list <- lapply(onek1k_list, function(x){
  # filter(x, P_VALUE <= 1) %>%
  filter(x, P_VALUE <= 5e-08) %>%  
    select(GENE_ID) %>%
    distinct(., GENE_ID)
})

# print number of significant genes per cell type
for (i in 1:length(unique_sign_genes_list)){
  cat("The number of significant genes for cell type ID <", names(unique_sign_genes_list[i]), "> is = ", nrow(unique_sign_genes_list[[i]]),".\n")
}

# print total number of unique genes across all cell types
total_unique_genes <- rbindlist(unique_sign_genes_list) %>%
  distinct(., GENE_ID) %>%
  nrow(.)
cat("The total number of unique genes across all cell types is = ", total_unique_genes, ".\n")

# print number of cell types
cat("The total number of cell types 'onek1k_list' is =",length(onek1k_list),".\n")

# keep only genes that have significant eQTLs and format table into LAVA format
for (i in 1:length(onek1k_list)){
  onek1k_list[[i]] <- onek1k_list[[i]] %>%
    dplyr::mutate(., N = 982, SE = "NA") %>%
    dplyr::filter(GENE_ID %in% unique_sign_genes_list[[i]]$GENE_ID) %>%
    dplyr::select(.,
      CELL_ID, # will remove this column after
      SNP = RSID,
      GENE = GENE_ID,
      CHR,
      BP = POS,
      A1 = A2,
      A2 = A1, # A2 is the effect allele (confirmed by looking at the allele count in the paper plots and the effect size in the database)
      BETA = SPEARMANS_RHO,
      SE,
      P = P_VALUE,
      N) %>%
    dplyr::mutate(., variant_ID = paste(CHR,BP,A1,A2,sep=":"))
}

# split per gene and name tibbles in list as: cell id + onek1k + gene id
onek1k_list_per_gene <- list()

for (i in 1:length(onek1k_list)){
  
  if (nrow(onek1k_list[[i]]) > 0){
    
    onek1k_list_per_gene[[i]] <- onek1k_list[[i]] %>%
      dplyr::group_split(GENE) %>%
      as.list() %>%
      set_names(., nm = map(.x = .,
                            ~glue("{first(.x$CELL_ID)}_onek1k_{first(.x$GENE)}")
                            )
                )
  }
  
}

# save one file per gene and per cell type
for(i in 1:length(onek1k_list_per_gene)){ # one list per cell type composed of lists
  for(j in 1:length(onek1k_list_per_gene[[i]])){ # each cell-type list composed of one tibble per gene
    onek1k_list_per_gene[[i]][[j]] %>%
      select(-CELL_ID) %>%
   fwrite(.,
      file =
        # here(project_dir,
        here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "onek1k",
             stringr::str_c(names(onek1k_list_per_gene[[i]][j]),lava_ext)), sep = "\t"
    )
  }
}
