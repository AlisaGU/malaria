#!/bin/bash
#SBATCH -n 4
#SBATCH -p hpc
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mem=150G
# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
chrom=chrom_template
export popu_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}"
export popu_variant_dir=$popu_output_dir/variant_3D7
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/0_get_ancestral"
## 只取核心区域的变异
# inputfile=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"
## 整个基因组的变异都取
inputfile=$popu_variant_dir/"${popu_symbol}_VQSLOD_gt0_${chrom}.vcf.gz"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# $code_dir/s1_2_get_mutation_load_for_each_site.R $inputfile #这个有问题，多次统计了spanning deletion
# $code_dir/s1_2_get_mutation_load_for_each_site_before_remove_spanning_deletion.R $inputfile
$code_dir/s1_2_get_mutation_load_for_each_site_remove_spanning_deletion.R $inputfile
# $code_dir/s1_2_get_mutation_load_for_each_site_remove_spanning_deletion_just_in_population.R $inputfile
# $code_dir/s1_2_get_mutation_load_for_each_site_remove_spanning_deletion_remove_popu_monomorphism.R $inputfile
# $code_dir/s1_2_get_mutation_load_for_each_site_nomatter_ancestor.R $inputfile
