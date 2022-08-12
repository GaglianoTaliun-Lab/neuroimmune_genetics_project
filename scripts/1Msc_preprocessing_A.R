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

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | B cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/B_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

B_cells <- fread(file_path) %>%
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

B_cells_list <- setNames(B_cells %>%
			 dplyr::group_split(GENE),
			 nm = str_c("Bcells_1MscBloodNL_",
				   B_cells$GENE %>%
				    unique() %>%
				    sort()
				    )
			 )

for(i in 1:length(B_cells_list)){

  fwrite(
    B_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
	stringr::str_c(names(B_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CD4 T cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/CD4T_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

CD4T_cells <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(SNPChr), EA = AlleleAssessed) %>%
 tidyr::separate(., col = SNPType, into = c("A1", "A2"), sep = "/") %>%
 dplyr::mutate(NEA = case_when(EA != A1 ~ A1, TRUE ~ A2)) %>%
 dplyr::filter(ProbeName %in% unique_significant_genes$ProbeName)	%>%
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

CD4T_cells_list <- setNames(CD4T_cells %>%
                         dplyr::group_split(GENE),
                         nm = str_c("CD4Tcells_1MscBloodNL_",
                                   CD4T_cells$GENE %>%
                                    unique() %>%
                                    sort()
                                    )
                         )

for(i in 1:length(CD4T_cells_list)){

  fwrite(
    CD4T_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(CD4T_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CD8 T cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/CD8T_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

CD8T_cells <- fread(file_path) %>%
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

CD8T_cells_list <- setNames(CD8T_cells %>%
                         dplyr::group_split(GENE),
                         nm = str_c("CD8Tcells_1MscBloodNL_",
                                   CD8T_cells$GENE %>%
                                    unique() %>%
                                    sort()
                                    )
                         )

for(i in 1:length(CD8T_cells_list)){

  fwrite(
    CD8T_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(CD8T_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Dendrytic cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/DC_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

DC_cells <- fread(file_path) %>%
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

DC_cells_list <- setNames(DC_cells %>%
                         dplyr::group_split(GENE),
                         nm = str_c("DCcells_1MscBloodNL_",
                                   DC_cells$GENE %>%
                                    unique() %>%
                                    sort()
                                    )
                         )

for(i in 1:length(DC_cells_list)){

  fwrite(
    DC_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(DC_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
