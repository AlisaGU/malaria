#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mem=35G

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
chrom=chrom_template
export popu_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}"
export popu_variant_dir=$popu_output_dir/variant_3D7
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/0_get_ancestral"
inputfile=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
$code_dir/s1_2_get_mutation_load_for_each_site.R $inputfile
