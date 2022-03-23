#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6
for i in $(echo "01 02 03 04 05 06 07 08 09 10 11 12 13 14" | tr " " "\n"); do
    echo $i
    before=$(wc -l WSEA_PASS_${i}_polysite.bed | awk '{print $1}')
    after=$(bcftools view -H WSEA_polymor_${i}.vcf.gz | wc -l | awk '{print $1}')
    diff=$(echo "$after-$before" | bc)
    por=$(echo "scale=2;$diff/$before" | bc)
    echo
done
