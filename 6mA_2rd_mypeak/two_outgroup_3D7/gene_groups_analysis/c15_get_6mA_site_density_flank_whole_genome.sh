#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
third_seq_6mA_site="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/3D7-6mA-9Cell_FRACgt20_COVgt25.bed"
gene_flank="gene_flank.bed"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome
$bedtools intersect -loj -a $gene_flank -b $third_seq_6mA_site >gene_flank_6mA_site_info
