#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_sample_list() {
    popu_symbol=$1
    $code_dir/s3_1_select_specific_population_sample.R $popu_symbol
}

extract_popu_vcf() {
    popu_output_dir=$1
    $code_dir/s3_2_extract_vcf.sh $popu_output_dir
}

compute_mutation_density() {
    popu_symbol=$1
    $code_dir/s3_3_compare_chrom.sh $popu_symbol
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd/whole_chromosome"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x

# popu_symbols="WAF CAF EAF WAF_CAF_EAF SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE"
popu_symbols="CAF"
for popu_symbol in $(echo ${popu_symbols} | tr " " "\n"); do
    popu_output_dir=$global_output_dir/$popu_symbol
    mkdir -p $popu_output_dir

    # get_sample_list $popu_symbol
    # extract_popu_vcf $popu_output_dir
    compute_mutation_density $popu_symbol
done
