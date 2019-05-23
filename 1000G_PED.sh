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
#PATHS
INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/
module load Python
module load VCFtools

### Construct haplotypes with vcftools application
vcftools --gzvcf ${INPUT}ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --plink --out ${INPUT}1000G_integrated_calls
