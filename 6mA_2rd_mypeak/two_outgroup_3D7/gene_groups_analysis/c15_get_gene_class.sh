#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7.noncore.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $working_dir
$bedtools intersect -a genes.bed -b $noncore -wa | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""noncore"}' >noncore.bed
$bedtools intersect -a genes.bed -b $noncore -wa -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""core"}' >core.bed
cat core.bed noncore.bed | sort -k1,1 -k2,2n >genes_class.bed
