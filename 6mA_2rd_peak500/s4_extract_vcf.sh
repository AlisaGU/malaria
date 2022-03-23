#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:


# VARIABLE NAMING TODO:
pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
peak500_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/top500"
motif_analysis_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif"
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
# chrom="01"
chrom=$1

variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
peak500_and_surr_position=$peak500_dir/"top500_and_surr.bed"
WSEA_samplelist=$motif_analysis_dir/"WSEA.sample.list"
WSEA_peak500_and_surr_QSlt0_variant=$pf6_variant_dir/"WSEA_peak500_surr_${chrom}.vcf.gz"
WSEA_peak500_and_surr_QSlt0_bialle_snp=$pf6_variant_dir/"WSEA_peak500_surr_bisnp_${chrom}.vcf.gz"
bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# extract specific region, specific sample, satisfying variants
$bcftools view --include 'VQSLOD>0.0' -S $WSEA_samplelist -R $peak500_and_surr_position -O z -o $WSEA_peak500_and_surr_QSlt0_variant --threads 4 $variant_file
$bcftools view --include 'TYPE="snp" && N_ALT=1' -O z -o $WSEA_peak500_and_surr_QSlt0_bialle_snp  $WSEA_peak500_and_surr_QSlt0_variant

# 