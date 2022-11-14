#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=get-maf-g1000eur
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

plink \
	--bfile ${project_dir}/reference_data/g1000_eur/g1000_eur \
	--freq \
	--out ${project_dir}/reference_data/g1000_eur/g1000_eur
