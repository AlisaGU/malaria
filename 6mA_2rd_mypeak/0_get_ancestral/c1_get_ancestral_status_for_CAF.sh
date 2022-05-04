#!/bin/bash
#SBATCH -n 8
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --mem=30G

# FUNCTIONS TODO:
subset_vcf() {
    chrom1=$1
    mkdir -p $popu_variant_dir
    echo $chrom1 | tr " " "\n" | $parallel -j 2 subset_vcf_for_one_chrom
}

subset_vcf_for_one_chrom() {
    chrom=$1
    variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
    popu_PASS_variant=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"

    $bcftools view --include 'FILTER="PASS"' -S $popu_samplelist -O z -o $popu_PASS_variant --threads 4 $variant_file 1>$popu_variant_dir/$chrom"_pass_select_sample.log" 2>&1
    $tabix $popu_PASS_variant
}

get_popu_and_outgroup_vcf() {
    # $code_dir/s1_1_get_popu_outgroup_sample_list_ref_popuid.R $popu_symbol $outgroup_symbol
    subset_vcf $chroms
}

get_mutation_load_for_all_chroms() {
    # echo "01" | $parallel -j 1 get_mutation_load_for_one_chrom
    echo $chroms | tr " " "\n" | grep -v "01" | $parallel -j 2 get_mutation_load_for_one_chrom
}

get_mutation_load_for_one_chrom() {
    chrom=$1
    inputfile=$popu_variant_dir/"${popu_symbol}_PASS_${chrom}.vcf.gz"
    $code_dir/s1_2_get_mutation_load_for_each_site.R $inputfile
}

export -f subset_vcf_for_one_chrom get_mutation_load_for_one_chrom
# VARIABLE NAMING TODO:
export popu_symbol="CAF"
export outgroup_symbol="OCE"

export code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/0_get_ancestral"
export popu_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}"
export pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
export popu_variant_dir=$popu_output_dir/variant
export popu_samplelist=$popu_output_dir/${popu_symbol}.sample.list

export chroms="01 02 03 04 05 06 07 08 09 10 11 12 13 14"

export bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
export tabix="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/tabix"
export parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# time get_popu_and_outgroup_vcf
# get_mutation_load_for_each_site
time get_mutation_load_for_all_chroms
