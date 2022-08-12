#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --job-name=get_overlap_GWAS_1MscTreatments
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# arguments:
# arg 1-x = GWAS phenotypes to include in the sample overlap matrix

# Rscript get_overlap_GWAS_QTL_per_egene.R
# Rscript get_overlap_GWAS_bulkQTL_per_egene.R

# Run the following with arrays (1-7):
# project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project/GWAS_summary_statistics/onek1k"
# cell_type=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/list_cell_types_subset.txt)
# Rscript get_overlap_GWAS_onek1k_per_egene.R ${cell_type}

Rscript get_overlap_GWAS_1Msc_treatment_per_egeneB.R
