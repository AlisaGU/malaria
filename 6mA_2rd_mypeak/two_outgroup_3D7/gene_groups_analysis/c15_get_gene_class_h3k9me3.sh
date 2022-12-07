#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/histone_methy/WT-H3K9me3-30h/BAM_broad_nomodel_dup1/WT-H3K9me3-30h_peaks.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $working_dir
$bedtools intersect -a genes.bed -b $noncore -wa | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""noncore"}' | sort | uniq >noncore_h3k9me3.bed
$bedtools intersect -a genes.bed -b $noncore -wa -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""core"}' | sort | uniq >core_h3k9me3.bed
cat core_h3k9me3.bed noncore_h3k9me3.bed | sort -k1,1 -k2,2n >genes_class_h3k9me3.bed
