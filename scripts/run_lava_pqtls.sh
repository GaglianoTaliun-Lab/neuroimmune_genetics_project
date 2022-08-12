#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=4-00:00:00
#SBATCH --array=1-15
#SBATCH --job-name=lava_pqtls
#SBATCH --output=slurm-%x-array-%a.out
#SBATCH --error=slurm-%x-array-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=6G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

chunk_name=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_pqtls.txt | awk '{print $1}')
chunk_start=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_pqtls.txt | awk '{print $2}')
chunk_end=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_pqtls.txt | awk '{print $3}')
Rscript lava_pqtls_target_loci.R ${chunk_name} ${chunk_start} ${chunk_end}
