#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mem=10G

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=chrom_template
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/across_species"
variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf"
# inputfile=$variant_dir/${chrom}.species3.include_indel.vcf.gz
inputfile=$variant_dir/${chrom}.noFilterBlock.include_indel.vcf.gz

# maf_prefix="P_falciparum.P_reichenowi.P_praefalciparum"
maf_prefix="P_falciparum.P_billcollinsi.P_reichenowi"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# $code_dir/s3_get_mutation_load_for_each_chrom_species2.R $inputfile

# $code_dir/s3_get_3d7_outgroup_A_consistent_region_include_indel.r $variant_dir/$chrom.${maf_prefix}.species3.fasta $chrom
# $code_dir/s3_get_3d7_outgroup_A_consistent_region_include_indel.r $variant_dir/$chrom.${maf_prefix}.noFilterBlock.fasta $chrom
$code_dir/s3_get_3d7_outgroupA_base_exist_region.r $variant_dir/$chrom.${maf_prefix}.noFilterBlock.fasta $chrom
