#!/bin/bash
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.log
# Usage: sbatch -J sms_spatial comparison/run_scmultisim.sbatch <spatial|velocity> <cells> <out dir>
source /etc/profile.d/lmod.sh; module load r/4.4.0
export R_LIBS_USER=/vast/djishnu/Palak/Rlibs/4.4
cd /vast/djishnu/Palak/GRN_ABM
/usr/bin/time -v Rscript comparison/run_scmultisim.R "$@" 2>&1 | grep -v "^\s*$"
