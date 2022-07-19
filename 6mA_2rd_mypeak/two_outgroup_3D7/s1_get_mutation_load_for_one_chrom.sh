#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mem=150G

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
chrom=chrom_template
export popu_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}"
export popu_variant_dir=$popu_output_dir/variant_two_group
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
# inputfile=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"
inputfile=$popu_variant_dir/"${popu_symbol}_VQSLOD_gt0_${chrom}.vcf.gz"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
$code_dir/s1_2_get_mutation_load_for_each_site.R $inputfile
# $code_dir/s1_2_get_mutation_load_for_each_site_before_correct_indel_site.R $inputfile
