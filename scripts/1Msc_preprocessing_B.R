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
#                                               | Megakaryocytes cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/megakaryocyte_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

mk_cells <- fread(file_path) %>%
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

MK_cells_list <- setNames(mk_cells %>%
			 dplyr::group_split(GENE),
			 nm = str_c("MK_1MscBloodNL_",
				   mk_cells$GENE %>%
				    unique() %>%
				    sort()
				    )
			 )

for(i in 1:length(MK_cells_list)){

  fwrite(
    MK_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(MK_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Monocyte cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/monocyte_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

monocyte_cells <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(SNPChr), EA = AlleleAssessed) %>%
 tidyr::separate(., col = SNPType, into = c("A1", "A2"), sep = "/") %>%
 dplyr::mutate(NEA = case_when(EA != A1 ~ A1, TRUE ~ A2)) %>%
 dplyr::mutate(N = 120) %>%
 dplyr::filter(ProbeName %in% unique_significant_genes$ProbeName) %>%
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

monocyte_cells_list <- setNames(monocyte_cells %>%
                         dplyr::group_split(GENE),
                         nm = str_c("monocyte_1MscBloodNL_",
                                   monocyte_cells$GENE %>%
                                    unique() %>%
                                    sort()
                                    )
                         )

for(i in 1:length(monocyte_cells_list)){

  fwrite(
    monocyte_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(monocyte_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Natural Killer cells |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/1M-scBloodNL/NK_eQTLs.txt")
fread(file_path, nrows = 6)

unique_significant_genes <- fread(file_path) %>%
  dplyr::filter(., PValue <= 5e-08) %>%
  select(ProbeName) %>%
  distinct(., ProbeName)

NK_cells <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(SNPChr), EA = AlleleAssessed) %>%
 tidyr::separate(., col = SNPType, into = c("A1", "A2"), sep = "/") %>%
 dplyr::mutate(NEA = case_when(EA != A1 ~ A1, TRUE ~ A2)) %>%
 dplyr::mutate(N = 120) %>%
 dplyr::filter(ProbeName %in% unique_significant_genes$ProbeName) %>%
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

NK_cells_list <- setNames(NK_cells %>%
                         dplyr::group_split(GENE),
                         nm = str_c("NKcells_1MscBloodNL_",
                                   NK_cells$GENE %>%
                                    unique() %>%
                                    sort()
                                    )
                         )

for(i in 1:length(NK_cells_list)){

  fwrite(
    NK_cells_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL",
        stringr::str_c(names(NK_cells_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
