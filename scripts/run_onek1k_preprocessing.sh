#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=10:00:00
#SBATCH --array=1-4,6,7,10
#SBATCH --job-name=onek1k_preprocessing
#SBATCH --output=slurm-%x-array-%a.out
#SBATCH --error=slurm-%x-array-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
cell_type1=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/onek1k/list_cell_types.txt | awk '{print $1}')
cell_type2=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/onek1k/list_cell_types.txt | awk '{print $2}')

Rscript onek1k_preprocessing.R ${cell_type1} ${cell_type2}
