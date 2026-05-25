#!/bin/bash
#SBATCH --job-name=BVD
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --time=05:00:00
#SBATCH --output=BVD_parameter_sweep.log

module purge

echo "Running BVD parameter sweep script on ${SLURM_CPUS_PER_TASK} cores"

module load R/4.4.2-gfbf-2024a

# Stop BLAS/OpenMP from oversubscribing cores when each future worker
# spawns its own linear-algebra threads. FlexiBLAS dispatches BLAS calls
# in the gfbf toolchain; OpenBLAS is its expected backend.
export OMP_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript inst/scripts/BVD_parameter_sweep.R
