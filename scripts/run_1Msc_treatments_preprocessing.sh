#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=10:00:00
#SBATCH --array=1-31
#SBATCH --job-name=GWAS_preprocessing
#SBATCH --output=slurm-%x-array-%a.out
#SBATCH --error=slurm-%x-array-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project/GWAS_summary_statistics/1M-scBloodNL/treatments_eqtls"

treatment=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/treatments_list.txt)

Rscript 1Msc_treatments_preprocessing.R ${treatment}
