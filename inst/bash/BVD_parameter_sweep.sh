#!/bin/bash
#SBATCH --job-name=BVD
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --time=05:00:00
#SBATCH --output=BVD_parameter_sweep.log

module purge

echo "Running BVD parameter sweep script"

module load R/4.4.2-gfbf-2024a

Rscript inst/scripts/BVD_parameter_sweep.R
