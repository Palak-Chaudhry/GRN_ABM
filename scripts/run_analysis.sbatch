#!/bin/bash
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch -J analysis scripts/run_analysis.sbatch results/runs/<name>
# Plots (spatial env) then RNA velocity (velo_py env, needs the embeddings written by the first step).
export HDF5_USE_FILE_LOCKING=FALSE OPENBLAS_NUM_THREADS=4
cd /vast/djishnu/Palak/GRN_ABM
/ihome/djishnu/pac271/.conda/envs/spatial/bin/python scripts/09_analysis_plots.py --run "$1" ${2:+--plots $2}
/ihome/djishnu/pac271/.conda/envs/velo_py/bin/python scripts/10_velocity_plots.py --run "$1"
