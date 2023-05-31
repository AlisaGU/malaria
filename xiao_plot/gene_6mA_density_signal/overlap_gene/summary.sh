#!/bin/bash
#SBATCH -n 8
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE

for bam in $(/usr/bin/ls *bam); do
    $samtools flagstat $bam --threads 8 >${bam}.summary
done
