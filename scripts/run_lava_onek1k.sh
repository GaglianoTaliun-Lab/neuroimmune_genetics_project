#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=3-00:00:00
#SBATCH --array=1-26
#SBATCH --job-name=MonoC_lava_qtls_onek1k
#SBATCH --output=slurm-%x-array-%a.out
#SBATCH --error=slurm-%x-array-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

#################### LAVA with ONEK1K QTLs:

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
chunk_name=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_onek1k.txt | awk '{print $1}')
chunk_start=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_onek1k.txt | awk '{print $2}')
chunk_end=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/parallelisation_list_onek1k.txt | awk '{print $3}')

# cell types: BIN, BMem, CD4ET, CD4NC, CD8ET, CD8NC, MonoC - choose one as first argument, defined below
cell_type="MonoC"

Rscript lava_onek1k_qtls_target_loci.R ${cell_type} ${chunk_name} ${chunk_start} ${chunk_end}
