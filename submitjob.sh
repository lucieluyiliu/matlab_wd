#!/bin/bas
#SBATCH --account=def-vihang01
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=128000M
#SBATCH -t 4:00:00
#SBATCH --job-name=ds1

module load matlab
matlab -nodesktop -nosplash < jackniferun.m > outputjacknife.log
