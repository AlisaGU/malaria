#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_average_variant_count_global_variant() {
    var_type=$1
    gene_groups=$2

    gene_groups1=$(echo $gene_groups | sed "s/ /_/g")
    mean_variant_count=$mutation_density_dir/${gene_groups1}_${var_type}_mean_variant_count

    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        for gene_group in $(echo $gene_groups | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $specific_gene_global_dir/${gene_group}/${gene_group}_genes.bed $var_type)"
        done

        echo $chrom_var_count >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region() {
    chrom=$1
    region_bed=$2
    var_type=$3

    gene_group=$(basename $region_bed | sed "s/_genes.bed//g")
    intermediate=$mutation_density_dir/${chrom}.${gene_group}_intermediate
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    $bedtools intersect -a $two_outgroup_consistent_region_in_core_and_noncore -b $region_bed | awk -v chrom=$chrom '$1==chrom' >$intermediate

    chrom_region_len=$(awk '{print $3-$2}' $intermediate | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b $intermediate | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    if [ $chrom_region_len -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_region_len=0
    fi

    echo $chrom_mutation_count" "$chrom_region_len

    rm -rf $intermediate
}

# VARIABLE NAMING TODO:
popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
two_outgroup_consistent_region_in_core_and_noncore_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore"
specific_gene_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
mutation_density_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_groups
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
variant_dir=$global_output_dir/$popu_symbol/variant_two_group
two_outgroup_consistent_region_in_core_and_noncore=$two_outgroup_consistent_region_in_core_and_noncore_dir/consistent_region_in_core_and_noncore.bed

gene_groups="VAR HDR invasion"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
get_average_variant_count_global_variant "allvar" "$gene_groups"
get_average_variant_count_global_variant "snps" "$gene_groups"
get_average_variant_count_global_variant "indels" "$gene_groups"
$code_dir/s9_plot.R $(echo $gene_groups | sed 's/ /_/g')
