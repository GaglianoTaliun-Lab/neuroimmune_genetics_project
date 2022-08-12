# neuro_immune_project
Scripts to compute genetic correlations and other analyses between neurodegenerative disorders and immune system

**PART A: correlations across GWAS datasets**
1. "prepare_input_info.R" = script to obtain the 'input.info.txt' file needed as input in LAVA, with the following columns: phenotype_name, N_cases, N_controls, path_to_lava.gz_file. Outputs a txt file tab separated.
2. "GWAS_preprocessing_*.R" = scripts (where * = A-D) to wrangle GWAS summary statistics and save them into the LAVA format. Outputs a gz compressed file for each GWAS sumstats.
3. "get_test_loci.R" = script to derive the LD blocks that will be used in LAVA. When changing the summary statistics included, make sure to include all possible column types across all GWAS sumstats (e.g. Z, beta, OR). Outputs a .loci file (LOCUS ID and positions) to use as loc_file in LAVA and an rds file.
4. "get_tsv_sumstats_for_ldsc_munge.R" = script to wrangle LAVA sumstats and prepare for munge sumstats in LDSC. It outputs a tsv table for munge_sumstats.py input.
5. "run_munge_sumstats.sh" = script to wrangle summary statistics for running LDSC across all GWAS summary statistics. Outputs a gz compressed file for each sumstats to use in LDSC.
6. "run_ldsc.sh" = script to run LDSC across all GWAS summary statistics (using the "list_of_GWAS_studies.txt" file) to get global genetic correlations and sample overlaps between datasets. Outputs a log file for each pairwise comparison. 
7. "get_corr_matrix.R" = script to process LDSC log file results and obtain the sample overlap matrix and global genetic correlations matrix. Outputs two text files with matrices.
8. "test_loci_plots.R" = script to generate plots that provide an overview of the tested loci across chromosomes and phenotypes.
9. "global_correlations_plot.R" = script to obtain a heatmap of the LDSC correlations results, highlighting significant correlations.
10. "lava.R" = script to run main univariate and bivariate LAVA analyses together. Outputs two rds files: results for the univariate analysis, and a second with all pairwise correlations performed if passed the univariate analysis.
11. "wrangle_lava_results.R" = script to extract significant results from the bivariate analyses.
12. "plots_lava.R" = generate plots of LAVA results (heatmap for each LD locus, edge diagrams, etc).
13. "prepare_PD_data_blauwendraat2019.R" = script to wrangle raw PD_blauwendraat2019 summary statistics, to include sample size per SNP depending on the number of studies.

>> Partial correlations and multiple regressions as post-hoc analysis of correlations <<
13. "lava_multi_reg.R" = script to run multiple regressions across loci which significant genetic correlation includes 3 or more traits.
14. "lava_partial_corr.R" = script to run partial correlations across loci which significant genetic correlation includes 3 or more traits.
15. "plots_lava_multi_reg.R" = script that plots the multiple regression coefficient across all models and a box plot of the multivariate r2. Uses raw LAVA results as input. 
16. "plots_lava_partial_corr.R" = script that creates a heatmap plot of the genetic correlations of each pair trait conditioned on a third trait, and comparing it with the unconditioned (i.e. initial bivariate) genetic correlations. Uses raw LAVA results as input.

**Given that I have run these scripts in Compute Canada (Cedar), I have also created shell scripts to run the R scripts. These are named as: "run_${Rscript_name}.sh"**

**PART B: correlations between GWAS datasets and eQTLs**
1. "1Msc_preprocessing_*.R" = scripts (where * = A-B) to wrangle QTL summary statistics and save them into the LAVA format. It outputs one file per genome-wide significantly expressed gene with all eQTLs, for each cell type. The +/- 100kb region will be defined by the gene_filtered.loci file.
2. “prepare_input_info_with_qtls.R” = using the file names for all GWAS and eQTL lava summary statistics, outputs a text file will case control sample sizes and NAs for the eQTLs.
3. "get_overlap_GWAS_QTL_per_egene.R" = script to obtain the sample overlap across all GWAS datasets and QTL datasets. However, it is not possible to run LDSC for QTLs, so we will assume that the overlap = 0, which is a reasonable assumption. It outputs a matrix of sample overlap estimates, with zeros in the case of QTLs.
4. "get_genic_regions_GWS.R" = script to obtain a +/-100kb window of the TSS and TES for any genome-wide significantly expressed gene across all cell types. It outputs a *.loci text file to be used as loc_file in LAVA, where instead of LD block ID, it is the ENSEMBL gene name. The LD block IDs are used as arguments to run the LAVA script.
5. "lava_with_qtls.R" = script to run main univariate and bivariate LAVA analyses together, with each QTL as target. It outputs one rds file per gene tested with the bivariate and univariate results. Needs to be run in parallel using chunks of loci of around 100 loci per script. It only tests cell types where there is a significant eQTL for the gene in question, and it performs the tests one GWAS trait at a time (to avoid losing SNPs if all GWAS traits are processed together).
6. "wrangle_lava_qtls.R" = script to format raw LAVA results with QTLs, and outputs one tsv table per test (univariate and bivariate). Additionally, it outputs a table of only the significant genetic correlations with a Bonferroni correction (# of bivariate tests) and another with an FDR < 0.01 correction.
7. “lava_qtl_plots.R” = script to plot the significant and all results from LAVA correlations across all cell types.

**Given that I have run these scripts in Compute Canada (Cedar), I have also created shell scripts to run the R scripts. These are named as: "run_${Rscript_name}.sh"**
