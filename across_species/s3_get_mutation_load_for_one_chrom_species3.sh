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
inputfile=$variant_dir/${chrom}.species3.vcf.gz
# maf_prefix="P_falciparum.P_reichenowi.P_praefalciparum"
maf_prefix="P_falciparum.P_billcollinsi.P_reichenowi"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
$code_dir/s3_get_mutation_load_for_each_chrom_species3.R $inputfile
$code_dir/s3_get_3d7_outgroup_A_consistent_region.r $variant_dir/$chrom.${maf_prefix}.species3.fasta $chrom
