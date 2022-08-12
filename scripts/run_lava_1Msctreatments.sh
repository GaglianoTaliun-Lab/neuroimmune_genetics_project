#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=6-00:00:00
#SBATCH --array=2
#SBATCH --job-name=lava_qtls_1Msc_treatments
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
treatment=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/1M-scBloodNL/treatments_eqtls/treatment_list_LAVA.txt)

for i in {1..7}
do
	chunk_name=$(sed -n ${i}p ${project_dir}/parallelisation_list_1Msctreatments.txt | awk '{print $1}')
	chunk_start=$(sed -n ${i}p ${project_dir}/parallelisation_list_1Msctreatments.txt | awk '{print $2}')
	chunk_end=$(sed -n ${i}p ${project_dir}/parallelisation_list_1Msctreatments.txt | awk '{print $3}')
	Rscript lava_1Msc_treatments.R ${treatment} ${chunk_name} ${chunk_start} ${chunk_end}
done
