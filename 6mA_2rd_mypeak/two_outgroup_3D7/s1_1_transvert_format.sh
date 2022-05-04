#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
input_file=input_template
export bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/OCE/africa_biallelic_variant
outfile=$(echo ${input_file} | sed 's/.vcf.gz//g')
$bcftools view -H $input_file | awk 'BEGIN{OFS="\t"}{for(i=10;i<=NF;i++){split($i,a,":");split(a[1],b,/[|/]/);printf("%s%s%s%s",b[1],OFS,b[2],OFS)};printf(ORS)}' >$outfile
