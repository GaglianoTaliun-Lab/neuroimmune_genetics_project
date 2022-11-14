# Description: run hyper-colocalization analysis after LAVA

# Packages -------------------------------------------------------

library(hyprcoloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(colochelpR)

# Arguments ----------------------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# default options:
p1 = 1e-04
p2 = 1e-04
p12 = 1e-05

# empty outputs:
dfs <- list()
GWAS_file_path <- vector()
QTL_file_path <- vector()

# Read files -------------------------------------------------------

test_gene = "SCIMP"

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc.txt"),
                         sep = "\t", header = T) %>%
  filter(., gene == test_gene)

# read file paths:
for (i in 1:nrow(gene_table)){
  
  GWAS_file_path[i] = stringr::str_c(here(project_dir, "colocalization", "GWAS_regional_sumstats"), "/", gene_table$trait[i], "_", gene_table$gene[i], ".tsv")
  QTL_file_path[i] = stringr::str_c(here(project_dir, "colocalization", "QTL_regional_sumstats"), "/onek1k_", gene_table$cell_type[i], "_", gene_table$gene[i], ".tsv")
  
}

# extract unique file paths (in case cell type or trait is repeated across rows in gene_table):
GWAS_file_path <- GWAS_file_path %>% unique()
QTL_file_path <- QTL_file_path %>% unique()
all_file_path <- c(GWAS_file_path, QTL_file_path)

# obtain vector of trait names for hyprcoloc input:
all_file_names <- basename(all_file_path) %>%
  stringr::str_remove(., ".tsv")

# read coloc data frames for each path and store in list:
for (i in 1:length(all_file_path)) {
  
  dfs[[i]] <- read.table(all_file_path[i], sep = "\t", header = T) %>%
    filter(., MAF > 0) %>%
    distinct(., SNP, .keep_all = TRUE)
  
}

# loop to get string of common SNPs across all datasets:
common_rsids <- as.data.frame(dfs[[1]]$SNP)
colnames(common_rsids) <- "SNP"

for (i in 2:length(dfs)) {
  
  common_rsids <- inner_join(common_rsids, dfs[[i]], by = "SNP") %>%
    dplyr::select(SNP)
  
}

cat("The number of common SNPs is:", nrow(common_rsids),".\n")

# create a list of datasets, keeping only common SNPs across all datasets:
dfs_subset_rsids <- lapply(dfs, function(x) {
  
  filter(x, SNP %in% common_rsids$SNP)
  
})

# create the vector for binary.outcomes needed for hyprcoloc (1 = binary; 0 = quant):
binary_outcomes = replicate(length(dfs_subset_rsids), 0)

for (i in 1:length(dfs_subset_rsids)) {
  
  if ("GWAS" %in% names(dfs_subset_rsids[[i]])) { # because all GWAS here are binary
    binary_outcomes[i] = 1
  }
  
}

cat("The binary outcomes vector is the following:\n")
binary_outcomes

# extract betas:
betas_coloc <- matrix(nrow = nrow(common_rsids), ncol = length(dfs_subset_rsids))

for (i in 1:length(dfs_subset_rsids)) {
  
  betas_coloc[,i] <- dfs_subset_rsids[[i]]$beta
  
}

colnames(betas_coloc) <- all_file_names

# extract SE:
se_coloc <- matrix(nrow = nrow(common_rsids), ncol = length(dfs_subset_rsids))

for (i in 1:length(dfs_subset_rsids)) {
  
  se_coloc[,i] <- dfs_subset_rsids[[i]]$se
  
}

colnames(se_coloc) <- all_file_names

# main hyprcoloc function:

hyprcoloc_res <- hyprcoloc(
  effect.est = betas_coloc,
  effect.se = se_coloc,
  binary.outcomes = binary_outcomes,
  trait.names = all_file_names,
  snp.id = common_rsids$SNP,
  sensitivity = FALSE,
  prior.12 = 1e-05,
  snpscores = TRUE
)

cat("Analysis done!\n")
cat("Hyprcoloc output looks like this:\n")
head(hyprcoloc_res)

hyprcoloc_df <- hyprcoloc_res$results %>%
  mutate(n_SNPs = nrow(common_rsids))

snp_scores <- data.frame("SNP" = common_rsids$SNP, "PP" = hyprcoloc_res$snpscores[[1]]) %>%
  arrange(., desc("PP"))

write.table(hyprcoloc_df, here(project_dir, "colocalization", "hyprcoloc_results", str_c("hyprcoloc_", test_gene,".tsv")), sep = "\t", row.names = F, quote = F)
write.table(snp_scores, here(project_dir, "colocalization", "hyprcoloc_results", str_c("hyprcoloc_snpscores_", test_gene,".tsv")), sep = "\t", row.names = F, quote = F)


