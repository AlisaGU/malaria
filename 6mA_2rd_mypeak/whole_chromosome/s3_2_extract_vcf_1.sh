#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd/whole_chromosome"
popu_output_dir=$1
chroms=$2
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
#01 02 03 04 05 06 07 08 09 10 11 12 13 14
for chrom in $(echo $chroms | tr " " "\n"); do
    $code_dir/ss3_2_extract_vcf_bychrom_specific_population.sh $chrom $popu_output_dir
done
