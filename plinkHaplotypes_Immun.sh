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
mkdir -p ${INPUT}Immun_pleios/
mkdir -p ${INPUT}Immun_snps/

module load Python

#### Crear haplotipos con Plink. Solo lo hace una vez por valor de R2
### Esta comanda de aqui abajo selecciona pleiotropias por el valor de r2 que tienen en la tabla.

sqlite3 ${INPUT}GWASpleiotropies.sqlite \
'SELECT DISTINCT SNPA,DiseaseA,RiskAllA,OnsetA,SNPB,DiseaseB,RiskAllB,onsetB,ID,CHR FROM filteredPairs WHERE R2 >= 0.8 AND (GroupA IN ("immune system disease") 
 OR GroupB IN ("immune system disease")) AND CHR != "" ;' > ${INPUT}Immun_pleios/pleios.txt

### Construct haplotypes with plink
cat ${INPUT}Immun_pleios/pleios.txt | while read line; do
    snpA=$(echo "$line" | cut -f 1);
    chr=$(echo "$line" | cut -f 10);
    echo $snpA;
    snpB=$(echo "$line" | cut -f 5);
    echo -e '*' ${snpA}'\t'${snpB} > ${INPUT}Immun_snps/snps.hlist;
    sed -i -e 's/ /\t/g' ${INPUT}Immun_snps/snps.hlist;
    ${BIN}plink --file ${VCF}chr${chr}_CEU_genotypes \
--hap ${INPUT}Immun_snps/snps.hlist \
--hap-freq \
--noweb \
--out ${INPUT}Immun_snps/${snpA}_${snpB};
grep -v LOCUS ${INPUT}Immun_snps/${snpA}_${snpB}.frq.hap | awk '{print $2,$3}' > ${INPUT}Immun_snps/${snpA}_${snpB}.fhtp
done

### 3 #### Create pleiotropies agon/antagon
cat ${INPUT}Immun_pleios/pleios.txt | while read line; do    
	snpA=$(echo "$line" | cut -f 1);
	snpB=$(echo "$line" | cut -f 5);
	riskHap=$(echo "$line" | awk -F"\t" '{print $3$7}');
### Python script to count pleiotropies in haplotypes
	python3 ./countPleiotropies.py ${INPUT}Immun_snps/${snpA}_${snpB}.fhtp ${riskHap} ${line} ${INPUT}Immun_snps/; ##### poner ${arx} para hacerlo como el original DONE. 	
    done

