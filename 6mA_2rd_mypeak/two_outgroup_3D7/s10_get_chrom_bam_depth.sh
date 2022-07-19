#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=chrom_template
types=types_template
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
genome_length=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length
chrom_length=$(grep $chrom $genome_length | awk '{print $NF}')
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/bam_depth

echo -e "${chrom}\t0\t${chrom_length}" | $samtools depth -b - ../3D7-T3_${types}.bam -a | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' >${chrom}.bam.${types}.depth.bed
