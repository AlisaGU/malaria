#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=$1

pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
motif_analysis_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif"
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"

WSEA_samplelist=$motif_analysis_dir/"WSEA.sample.list"

variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
WSEA_PASS_variant=$pf6_variant_dir/"WSEA_PASS_${chrom}.vcf.gz"
WSEA_polysite_bed=$pf6_variant_dir/"WSEA_PASS_${chrom}_polysite.bed"
WSEA_snp_bed=$pf6_variant_dir/WSEA_PASS_${chrom}_snp_polysite.bed
WSEA_indel_bed=$pf6_variant_dir/WSEA_PASS_${chrom}_indel_polysite.bed

WSEA_polymorphic_variant=$pf6_variant_dir/"WSEA_polymor_${chrom}.vcf"
WSEA_snp=$pf6_variant_dir/"WSEA_snp_${chrom}.vcf.gz"
WSEA_indel=$pf6_variant_dir/"WSEA_indel_${chrom}.vcf.gz"

bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
tabix="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/tabix"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
bgzip="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/bgzip"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# $bcftools view --include 'FILTER="PASS"' -S $WSEA_samplelist -O z -o $WSEA_PASS_variant --threads 4 $variant_file
# $tabix $WSEA_PASS_variant

# $code_dir/s4_get_heterosite_bed.R $WSEA_PASS_variant

awk '{print $3}' $WSEA_polysite_bed >polysite
rm -rf $WSEA_polymorphic_variant
$bcftools view -h $WSEA_PASS_variant >$WSEA_polymorphic_variant
$bcftools view -H $WSEA_PASS_variant | csvtk grep -H -t -f 2 -P polysite >>$WSEA_polymorphic_variant
$bgzip $WSEA_polymorphic_variant
rm -rf polysite
$tabix ${WSEA_polymorphic_variant}.gz

$bcftools view --types "snps" -O z -o $WSEA_snp --threads 4 ${WSEA_polymorphic_variant}.gz
$tabix $WSEA_snp
$bcftools view -H $WSEA_snp | awk '{print $2}' >polysite
csvtk grep -H -t -f 3 -P polysite $WSEA_polysite_bed >$WSEA_snp_bed
rm -rf polysite

$bcftools view --types "indels" -O z -o $WSEA_indel --threads 4 ${WSEA_polymorphic_variant}.gz
$tabix $WSEA_indel
$bcftools view -H $WSEA_indel | awk '{print $2}' >polysite
csvtk grep -H -t -f 3 -P polysite $WSEA_polysite_bed >$WSEA_indel_bed
rm -rf polysite
