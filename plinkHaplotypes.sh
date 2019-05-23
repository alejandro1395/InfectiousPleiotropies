#!/bin/bash

#set the job name
#SBATCH --job-name=db_pleio
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=23:00

#run the application
#PATHS
INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/results/ 
VCF=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/
BIN=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/bin/plink-1.07-x86_64/
mkdir -p ${INPUT}InfecImmun_pleios/
mkdir -p ${INPUT}InfecImmun_snps/

module load Python
module load modulepath/noarch

#### Crear haplotipos con Plink. Solo lo hace una vez por valor de R2
### Esta comanda de aqui abajo selecciona pleiotropias por el valor de r2 que tienen en la tabla.
sqlite3 ${INPUT}GWASpleiotropies.sqlite 'SELECT DISTINCT SNPA,DiseaseA,RiskAllA,OnsetA,SNPB,DiseaseB,RiskAllB,onsetB,ID FROM filteredPairs WHERE R2 >= 0.8 AND ((GroupA IN ("infectious disease") AND GroupB IN ("immune system disease")) OR (GroupA IN ("immune system disease") AND GroupB IN ("infectious disease")));' > ${INPUT}InfecImmun_pleios/pleios.txt

### Construct haplotypes with plink
cat ${INPUT}InfecImmun_pleios/pleios.txt | while read line; do
    snpA=$(echo "$line" | cut -f 1);
    echo $snpA;
    snpB=$(echo "$line" | cut -f 5);
    echo -e '*' ${snpA}'\t'${snpB} > ${INPUT}InfecImmun_snps/snps.hlist;
    sed -i -e 's/ /\t/g' ${INPUT}InfecImmun_snps/snps.hlist;
    ${BIN}plink --bfile ${VCF}hapmap_CEU_r23a \
--hap ${INPUT}InfecImmun_snps/snps.hlist \
--hap-freq \
--noweb \
--out ${INPUT}InfecImmun_snps/${snpA}_${snpB};
grep -v LOCUS ${INPUT}InfecImmun_snps/${snpA}_${snpB}.frq.hap | awk '{print $2,$3}' > ${INPUT}InfecImmun_snps/${snpA}_${snpB}.fhtp
done

### 3 #### Create pleiotropies agon/antagon

