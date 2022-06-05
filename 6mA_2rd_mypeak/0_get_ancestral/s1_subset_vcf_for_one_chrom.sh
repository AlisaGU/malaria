#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=chrom_template
popu_symbol=popu_symbol_template
export pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
export popu_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}"
export popu_variant_dir=$popu_output_dir/variant_3D7
export popu_samplelist=$popu_output_dir/${popu_symbol}.sample.list

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
popu_PASS_variant=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"

$bcftools view --include 'FILTER="PASS"' -S $popu_samplelist -O z -o $popu_PASS_variant --threads 4 $variant_file 1>$popu_variant_dir/$chrom"_pass_select_sample.log" 2>&1
