#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=lava_multi_reg
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

##################### LAVA for GWAS traits
# arguments:
# arg 1 = lava run version (e.g. v1, v2, v3)
# arg 2-x = phenotypes to include in the LAVA analysis

Rscript lava_multi_reg.R v3 AD_schwartzentruber2021 CD_delange2017 MS_imsgc2019 UC_delange2017

