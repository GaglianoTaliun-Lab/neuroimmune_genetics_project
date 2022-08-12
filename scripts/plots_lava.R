library(here)
library(stringr)
library(tidyverse)
library(circlize)
library(tidygraph)
library(ggraph)
library(rtracklayer)
library(gghighlight)

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# arguments
args <- commandArgs(TRUE)

# GWAS phenotypes that were included in the LAVA run:
datasets <- array(args[1:length(args)]) %>% sort()
lava_output_name <- str_c(datasets, collapse = ":")

bivar_all <- read.table(here(project_dir,"lava_results", "tsv_univ_bivar_GWAS","phenotypes_all_results_bivar.tsv"), sep = "\t", header = T)
bivar_bonf <- read.table(here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_significant_results_bivar.tsv"), sep = "\t", header = T)

source(here(project_dir, "plots_RHR.R"))

# get number of bivariate tests
bivar_rds <- readRDS(here(project_dir,"lava_results","lava_GWAS",stringr::str_c(lava_output_name,".bivar.lava.rds")))

# get number of bivariate tests
ntests = 0
for (i in 1:length(bivar_rds)) {
  condition = nrow(bivar_rds[[i]]) > 0
  if (condition == "TRUE" && length(condition) != 0) {
    ntests <- ntests + nrow(bivar_rds[[i]])
  }
}

# get pvalue bonferroni threshold
pvalue_bivar = 0.05/ntests

fct_disease = factor(datasets)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### scatter plot of n snps and sample size

bivar_bonf %>% 
  dplyr::count(phen1) %>% 
  dplyr::bind_rows(
    bivar_bonf %>% 
      dplyr::count(phen2) %>% 
      dplyr::rename(phen1 = phen2)
  ) %>% 
  dplyr::group_by(phen1) %>% 
  dplyr::summarise(
    n_bivar = sum(n)
  ) %>% 
  dplyr::inner_join(
    read.table(
      here(project_dir,"info_files","input.info.txt"), sep = "\t", header = T) %>% 
      dplyr::mutate(
        n_individ = cases + controls
      ),
    by = c("phen1" = "phenotype")
  ) %>% 
  ggplot(aes(x = n_bivar, y = n_individ, colour = phen1)) +
  geom_point(size = 4) +
  scale_colour_brewer(
    palette = "BrBG",
    type = "div",
    name = "GWAS"
  ) +
  labs(
    x = "Number of significant local genetic correlations",
    y = "Sample size"
  ) + 
  theme(axis.title.x = element_text(size = 9),
	axis.title.y = element_text(size = 9),
	axis.text.x = element_text(size = 6),
	axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11)) + theme_rhr

ggsave(here(project_dir, "lava_results", "figures", "phenotypes_scatter_plot.pdf"), width = 40, height = 30, units = "cm")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Chord diagram:

# plot_bivar_chord_diagram(
#  bivar_corr = bivar_bonf,
#  fct_phen = fct_disease,
#  palette = RColorBrewer::brewer.pal(n = 6, name = "BrBG")
#)
# ggsave(here(project_dir, "lava_results", "figures",stringr::str_c("phenotypes_chord_diagram.",run_version,".pdf")), width = 40, height = 30, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Heat map per LD block

bivar_all %>% 
  dplyr::filter(locus %in% unique(bivar_bonf$locus)) %>% 
  dplyr::mutate(
    rho_fill = 
      case_when(
        p < pvalue_bivar ~ round(rho, 2)
      )
  ) %>% 
  ggplot(
    aes(
      x = phen1,
      y = phen2,
      fill = rho_fill,
      label = round(rho, 2)
    )
  ) +
  geom_tile(colour = "black") +
  geom_text(
    size = 2
  ) +
  facet_wrap(vars(locus), ncol = 5) +
  scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "#cccccc", limits = c(-1, 1)) + theme_rhr +
  theme(axis.text.x = element_text(angle = 90))

ggsave(here(project_dir, "lava_results", "figures", "phenotypes_heatmaps_per_LDblock.pdf"), width = 35, height = 65, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Edge Diagram

plot_edge_diagram(
  bivar_corr = 
    bivar_bonf,
  phen = fct_disease, 
  ncol = 3,
)

ggsave(here(project_dir, "lava_results", "figures", "phenotypes_edge_diagram.pdf"), width = 40, height = 50, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Annotation Plots per LD block

ref <- rtracklayer::import(here(project_dir,"reference_data","Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse") 
ref <- ref[ref$type == "gene"]

loci_gr <-
  bivar_bonf %>%
  dplyr::count(locus, chr, start, stop, n_snps) %>%
  dplyr::arrange(locus) %>% 
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "stop"
  )

fig_list = vector(mode = "list", length = length(loci_gr))

for(i in 1:length(loci_gr)){
  
  fig_list[[i]] <- 
    plot_locus(
    locus_gr = loci_gr[i], 
    ref = ref
    )
  
  names(fig_list)[i] <- str_c("locus_", loci_gr[i]$locus)
  
}

pdf(here(project_dir, "lava_results", "figures", "phenotypes_LDblock_annotations.pdf"))
fig_list
dev.off()

