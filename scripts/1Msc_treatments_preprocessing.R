library(devtools)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(here)

# devtools::install_github("RHReynolds/colochelpR")
# devtools::install_github("RHReynolds/rutils")
library(colochelpR)
library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args <- commandArgs(TRUE)
treatment <- args[1]
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | treatment eQTLs for LAVA |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- here(project_dir, "GWAS_summary_statistics", "1M-scBloodNL", "treatments_eqtls", str_c(treatment,".txt.gz"))
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

cat("The number of unique genes for '",treatment,"' is:",nrow(unique_significant_genes),".\n")

one_cell_type <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(SNPChr), EA = AlleleAssessed) %>%
 tidyr::separate(., col = SNPType, into = c("A1", "A2"), sep = "/") %>%
 dplyr::mutate(NEA = case_when(EA != A1 ~ A1, TRUE ~ A2)) %>%
 dplyr::filter(ProbeName %in% unique_significant_genes$ProbeName) %>%
 dplyr::mutate(N = 120) %>%
 dplyr::select(
   SNP = SNPName,
   GENE = ProbeName,
   CHR,
   BP = SNPChrPos,
   A1 = EA,
   A2 = NEA,
   Z = OverallZScore,
   P = PValue,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

treatment_name = str_remove_all(treatment, "_eQTLs")

cell_type_list <- setNames(one_cell_type %>%
			 dplyr::group_split(GENE),
			 nm = str_c(treatment_name,"_1MscBloodNL_",
				   one_cell_type$GENE %>%
				    unique() %>%
				    sort()
				    )
			 )

for(i in 1:length(cell_type_list)){

  fwrite(
    cell_type_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL_treatments",
	stringr::str_c(names(cell_type_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
