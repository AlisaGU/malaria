#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=$1
popu_output_dir=$2

code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"

popu_symbol=$(basename $popu_output_dir)
popu_samplelist=$popu_output_dir/${popu_symbol}.sample.list
variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
popu_variant_dir=$popu_output_dir/variant

popu_PASS_variant=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"
popu_polysite_bed=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}_polysite.bed"
popu_snp_bed=$popu_variant_dir/${popu_symbol}_PASS_${chrom}_snp_polysite.bed
popu_indel_bed=$popu_variant_dir/${popu_symbol}_PASS_${chrom}_indel_polysite.bed

popu_polymorphic_variant=$popu_variant_dir/"${popu_symbol}_polymor_${chrom}.vcf"
popu_snp=$popu_variant_dir/"${popu_symbol}_snp_${chrom}.vcf.gz"
popu_indel=$popu_variant_dir/"${popu_symbol}_indel_${chrom}.vcf.gz"

bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
tabix="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/tabix"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
bgzip="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/bgzip"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# mkdir -p $popu_variant_dir
# $bcftools view --include 'FILTER="PASS"' -S $popu_samplelist -O z -o $popu_PASS_variant --threads 4 $variant_file
# $tabix $popu_PASS_variant

# $code_dir/sss3_2_1_get_heterosite_bed.R $popu_PASS_variant

# awk '{print $3}' $popu_polysite_bed >polysite
# rm -rf $popu_polymorphic_variant
# $bcftools view -h $popu_PASS_variant >$popu_polymorphic_variant
# $bcftools view -H $popu_PASS_variant | csvtk grep -H -t -f 2 -P polysite >>$popu_polymorphic_variant
$bgzip $popu_polymorphic_variant
rm -rf polysite
$tabix ${popu_polymorphic_variant}.gz

$bcftools view --types "snps" -O z -o $popu_snp --threads 4 ${popu_polymorphic_variant}.gz
$tabix $popu_snp
$bcftools view -H $popu_snp | awk '{print $2}' >polysite
csvtk grep -H -t -f 3 -P polysite $popu_polysite_bed >$popu_snp_bed
rm -rf polysite

$bcftools view --types "indels" -O z -o $popu_indel --threads 4 ${popu_polymorphic_variant}.gz
$tabix $popu_indel
$bcftools view -H $popu_indel | awk '{print $2}' >polysite
csvtk grep -H -t -f 3 -P polysite $popu_polysite_bed >$popu_indel_bed
rm -rf polysite
