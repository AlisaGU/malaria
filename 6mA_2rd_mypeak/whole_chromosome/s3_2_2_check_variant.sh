#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom="chrom_template"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/"
popu_symbol="popu_template"
variant_dir=$global_output_dir/$popu_symbol/variant
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
cd $variant_dir
echo $chrom
echo "all variants:"
echo "bed: "$(wc -l ${popu_symbol}_PASS_${chrom}_polysite.bed)
echo "vcf: "$(bcftools view -H ${popu_symbol}_polymor_${chrom}.vcf.gz | wc -l)
echo "SNPs:"
echo "bed: "$(wc -l ${popu_symbol}_PASS_${chrom}_snp_polysite.bed)
echo "vcf: "$(bcftools view -H ${popu_symbol}_snp_${chrom}.vcf.gz | wc -l)
echo "Indels:"
echo "bed: "$(wc -l ${popu_symbol}_PASS_${chrom}_indel_polysite.bed)
echo "vcf: "$(bcftools view -H ${popu_symbol}_indel_${chrom}.vcf.gz | wc -l)
