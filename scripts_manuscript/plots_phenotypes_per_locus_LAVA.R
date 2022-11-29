library(tidyverse)
library(here)
library(gghighlight)
library(stringr)
library(reshape2)
library(corrplot)
library(gridExtra)
library(ggpubr)

# Arguments ----------------------------------------------------------------------------------------

project_dir <- "Documents/research-projects/neuro_immune"

bivar_all <- read.table(here(project_dir, "lava_results", "phenotypes", "phenotypes_all_results_bivar.tsv"), sep = "\t", header = T)
bivar_sign <- read.table(here(project_dir, "lava_results", "phenotypes", "phenotypes_significant_results_bivar.tsv"), sep = "\t", header = T)

date=str_remove_all(Sys.Date(), "-")

bonferroni = 2.629e-05

######### -----------------------------------------------------

# locus to plot:

# - 118
bivar_to_plot <-
  bivar_all %>%
  dplyr::filter(., locus == 118) %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         fill_rho = case_when(
           p < 0.05 ~ round(rho, 2)
         ),
         sig = case_when(
           p < bonferroni ~ "**"
         ),
         nom_sig = case_when(
           p < 0.05 & p >= bonferroni ~ ".",
         ),
         not_sig = case_when(
           p >= 0.05 ~ "ns"
         ),
         chr_pos = str_c("chr ",chr,": ", start, "-", stop)
  )

plot_118 <- ggplot(bivar_to_plot, aes(x=phen2_clean, y=phen1_clean)) + 
  geom_count(colour = "black", shape = 15, size = 8) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "", y = "", title = "") +
  guides(colour = guide_colourbar(title = "Regional rg")) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14)) +
  geom_text(aes(label=sig), size = 4, vjust = 0.75) +
  geom_text(aes(label=nom_sig), size = 5, vjust = 0.1) +
  geom_text(aes(label=not_sig), size = 4, vjust = 0.4) +
  coord_flip() +
  facet_wrap(~chr_pos)

# - 752
bivar_to_plot <-
  bivar_all %>%
  dplyr::filter(., locus == 752) %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         fill_rho = case_when(
           p < 0.05 ~ round(rho, 2)
         ),
         sig = case_when(
           p < bonferroni ~ "**"
         ),
         nom_sig = case_when(
           p < 0.05 & p >= bonferroni ~ ".",
         ),
         not_sig = case_when(
           p >= 0.05 ~ "ns"
         ),
         chr_pos = str_c("chr ",chr,": ", start, "-", stop)
  )

plot_752 <- ggplot(bivar_to_plot, aes(x=phen2_clean, y=phen1_clean)) + 
  geom_count(colour = "black", shape = 15, size = 8) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "", y = "", title = "") +
  guides(colour = guide_colourbar(title = "Regional rg")) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14)) +
  geom_text(aes(label=sig), size = 4, vjust = 0.75) +
  geom_text(aes(label=nom_sig), size = 5, vjust = 0.1) +
  geom_text(aes(label=not_sig), size = 4, vjust = 0.4) +
  coord_flip() +
  facet_wrap(~chr_pos)

# - 1146
bivar_to_plot <-
  bivar_all %>%
  dplyr::filter(., locus == 1146) %>%
  tidyr::separate(phen1, c("phen1_clean", NA), sep = "_") %>%
  tidyr::separate(phen2, c("phen2_clean", NA), sep = "_") %>%
  mutate(., order_phen1 = factor(phen1_clean, levels = c("AD","PD","ALS","SCZ","MS","UC","CD")),
         fill_rho = case_when(
           p < 0.05 ~ round(rho, 2)
         ),
         sig = case_when(
           p < bonferroni ~ "**"
         ),
         nom_sig = case_when(
           p < 0.05 & p >= bonferroni ~ ".",
         ),
         not_sig = case_when(
           p >= 0.05 ~ "ns"
         ),
         chr_pos = str_c("chr ",chr,": ", start, "-", stop)
  )

plot_1146 <- ggplot(bivar_to_plot, aes(x=phen2_clean, y=phen1_clean)) + 
  geom_count(colour = "black", shape = 15, size = 8) +
  geom_count(aes(colour = fill_rho), shape = 15, size = 7) +
  coord_fixed() +
  theme_bw() +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1,0,1), na.value = "grey") +
  labs(x = "", y = "", title = "") +
  guides(colour = guide_colourbar(title = "Regional rg")) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14)) +
  geom_text(aes(label=sig), size = 4, vjust = 0.75) +
  geom_text(aes(label=nom_sig), size = 5, vjust = 0.1) +
  geom_text(aes(label=not_sig), size = 4, vjust = 0.4) +
  coord_flip() +
  facet_wrap(~chr_pos)

ggarrange(plot_118, plot_752, plot_1146, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

ggsave(here(project_dir, "lava_plots", "lava_GWAS", stringr::str_c("phenotypes_heatmap_locus_118_752_1146_",date,".jpg")),
       width = 40, height = 10, units = "cm")
