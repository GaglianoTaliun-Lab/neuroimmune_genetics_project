#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=10:00:00
#SBATCH --job-name=GWAS_preprocessing
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# Rscript GWAS_preprocessing_tmp.R
# Rscript 1Msc_preprocessing_A.R
# Rscript 1Msc_preprocessing_B.R
# Rscript eQTLGen_preprocessing.R
Rscript pqtls_preprocessing.R
