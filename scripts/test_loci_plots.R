library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(Homo.sapiens)

run_version = "v3"

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# source("functions_granges.R")

loci <- read.table(here(project_dir,"test_loci",stringr::str_c("gwas_filtered.loci.",run_version)), sep = "\t", header = T)

loci %>% 
  dplyr::mutate(width = STOP - START) %>% 
  ggplot(
    aes(
      x = as.factor(CHR),
      y = width)
  ) +
  geom_boxplot() +
  scale_y_sqrt(n.breaks = 10) +
  labs(x = "Chromosome",
       y = "Locus width (bp, square root scale)")
ggsave(here(project_dir,"test_loci","figures",stringr::str_c("width_gws_LD_blocks.",run_version,".pdf")), width = 40, height = 30, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

loci %>% 
  dplyr::count(CHR) %>% 
  ggplot(
    aes(
      x = forcats::fct_reorder(.f = as.factor(CHR), 
                               .x = n, 
                               .fun = median, 
                               .desc = TRUE),
      y = n)
  ) +
  geom_col() +
  labs(x = "Chromosome",
       y = "Number of LD blocks") +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
ggsave(here(project_dir,"test_loci","figures",stringr::str_c("number_gws_LD_blocks.",run_version,".pdf")), width = 40, height = 30, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

gwas_loci <- 
  readRDS(here(project_dir,"test_loci",stringr::str_c("gwas_filtered_loci.",run_version,".rds")))
levels_gwas <- count(gwas_loci, GWAS) %>% .[,1]

gwas_loci <- gwas_loci %>%
  dplyr::mutate(
    GWAS = 
      fct_relevel(
        GWAS,
        levels_gwas
      )
  )

##### annotation of LD blocks (https://support.bioconductor.org/p/67118/)

gwas_loci <- gwas_loci %>% mutate(chr_LDblock = paste0("chr", LD_CHR))
unique_loci <- gwas_loci[, c(18, 14, 15, 17)] %>% unique(.)
write.table(unique_loci, "unique_LDblocks.bed", sep = "\t", row.names = F, quote = F)

LD_loci <- makeGRangesFromDataFrame(unique_loci,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="chr_LDblock",
                         start.field="LD_start",
                         end.field="LD_end")

####################################################################

gwas_loci %>% 
  dplyr::count(LD_CHR, LD_LOC, GWAS) %>% 
  dplyr::inner_join(
    gwas_loci %>% 
      dplyr::count(LD_CHR, LD_LOC, name = "n_total")
  ) %>% 
  dplyr::mutate(
    LD_LOC = as.factor(LD_LOC)
  ) %>% 
  dplyr::slice_max(order_by = n_total, n = 25) %>% 
  ggplot(
    aes(
      x = n,
      y = forcats::fct_reorder(.f = LD_LOC, 
                               .x = n, 
                               .fun = sum, 
                               .desc = FALSE),
      fill = GWAS) 
  ) +
  geom_col(colour = "black") +
  scale_fill_brewer(
    palette = "BrBG",
    type = "div",
    name = "GWAS"
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  labs(x = "Number of SNPs",
       y = "LD block")
ggsave(here(project_dir,"test_loci","figures",stringr::str_c("SNPs_per_LD_block.",run_version,".pdf")), width = 40, height = 30, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

gwas_loci %>% 
  dplyr::count(LD_CHR, GWAS) %>% 
  dplyr::mutate(
    LD_CHR = as.factor(LD_CHR)
  ) %>% 
  ggplot(
    aes(
      x = fct_rev(LD_CHR),
      y = n,
      fill = GWAS
    )
  ) +
  geom_col(
    position = position_dodge2(preserve = "single"),
    colour = "black"
  ) +
  scale_fill_brewer(
    palette = "BrBG",
    type = "div"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Chromosome",
       y = "Number of SNPs") + 
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="top") +
  coord_flip()
ggsave(here(project_dir,"test_loci","figures",stringr::str_c("SNPs_per_chromosome.",run_version,".pdf")), width = 40, height = 30, units = "cm")
