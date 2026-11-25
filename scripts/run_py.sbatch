#!/bin/bash
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch --job-name=<name> scripts/run_py.sbatch scripts/<script>.py [args]
export HDF5_USE_FILE_LOCKING=FALSE OPENBLAS_NUM_THREADS=2
cd /vast/djishnu/Palak/GRN_ABM
/ihome/djishnu/pac271/.conda/envs/spatial/bin/python "$@"
