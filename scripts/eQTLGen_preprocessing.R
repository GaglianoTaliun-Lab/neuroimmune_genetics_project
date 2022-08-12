library(devtools)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(here)

library(colochelpR)
library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Bulk eQTLs |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/eQTLGenConsortium/eQTLGen_all.txt")
fread(file_path, nrows = 6)

# obtain significant genes from 1MscBlood sc-eQTLs
significant_genes <- read.table(here(project_dir, "test_loci", "window_100000", "unique_GWS_genes.txt"), sep = "\t", header = T)

eqtls <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(SNPChr)) %>%
 dplyr::filter(Gene %in% significant_genes$GENE) %>%
 dplyr::select(
   SNP = SNP,
   GENE = Gene,
   CHR,
   BP = SNPPos,
   A1 = AssessedAllele,
   A2 = OtherAllele,
   Z = Zscore,
   P = BonferroniP,
   N = NrSamples) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

eqtls_list <- setNames(eqtls %>%
			 dplyr::group_split(GENE),
			 nm = str_c("eQTLGen_",
				   eqtls$GENE %>%
				    unique() %>%
				    sort()
				    )
			 )

for(i in 1:length(eqtls_list)){

  fwrite(
    eqtls_list[[i]],
    file =
      here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","eQTLGen",
	stringr::str_c(names(eqtls_list)[i],lava_ext)
      ),
    sep = "\t"
  )

}
