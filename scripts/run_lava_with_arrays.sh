#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2-00:00:00
#SBATCH --array=5
#SBATCH --job-name=lava_qtls_array_B
#SBATCH --output=slurm-%x-array-%a.out
#SBATCH --error=slurm-%x-array-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

#################### LAVA with QTLs:
Rscript lava_qtls_target_QTL_loci_B.R ${SLURM_ARRAY_TASK_ID}
