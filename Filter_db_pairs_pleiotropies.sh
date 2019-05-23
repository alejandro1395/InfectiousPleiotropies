#!/bin/bash

#set the job name
#SBATCH --job-name=db_pleio
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=23:00
 
# run the application

module load Python
python3 Filter_db_pairs_pleiotropies.py




