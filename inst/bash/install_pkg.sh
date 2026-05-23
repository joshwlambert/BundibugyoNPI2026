#!/bin/bash
#SBATCH --job-name=install_pkg
#SBATCH --ntasks=1
#SBATCH --mem=5GB
#SBATCH --time=01:00:00
#SBATCH --output=install_pkg.log


echo "Installing BundibugyoNPI2026 and its dependencies"

module load R/4.4.2-gfbf-2024a

Rscript -e "remotes::install_github('joshwlambert/BundibugyoNPI2026')"
