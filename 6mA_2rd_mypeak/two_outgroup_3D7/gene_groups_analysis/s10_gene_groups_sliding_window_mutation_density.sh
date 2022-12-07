#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_average_variant_count_global_variant() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for window in $(seq 1 50); do
        gene_count=$(grep "\_${window}$" ${gene}_genes.window.bed | wc -l)
        window_gene_mutation_density=""
        for gene_order in $(seq 1 $gene_count); do
            chrom=$(grep "\_${window}$" ${gene}_genes.window.bed | awk -v gene_order=$gene_order 'NR==gene_order{print $1}')
            bed=$(grep "\_${window}$" ${gene}_genes.window.bed | awk -v gene_order=$gene_order 'NR==gene_order{print $1"\t"$2"\t"$3}')
            window_gene_mutation_density+=" $(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")"
        done
        echo $window_gene_mutation_density >>$mean_variant_count
    done
}

get_average_variant_count_global_variant_peak_window() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_peak_window_${window_length}bp
    rm -rf $mean_variant_count

    for window in $(awk '{print $NF}' ${gene}_genes.peak.window.${window_length}bp.bed | sort -n | uniq); do
        gene_count=$(grep -e "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | wc -l)
        window_gene_mutation_density=""
        for gene_order in $(seq 1 $gene_count); do
            chrom=$(grep -e "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | awk -v gene_order=$gene_order 'NR==gene_order{print $1}')
            bed=$(grep -e "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | awk -v gene_order=$gene_order 'NR==gene_order{print $1"\t"$2"\t"$3}')
            window_gene_mutation_density+=" $(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")"
        done
        echo $(echo $"$window_gene_mutation_density" | grep -v "^$" | awk '{s+=$1}END{print s}') $(echo $"$window_gene_mutation_density" | grep -v "^$" | awk '{s+=$2}END{print s}') >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{print $3-$2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    # chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    # chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    # echo $chrom_region_mutation_count_mean1
    echo "$chrom_mutation_count $chrom_region_len"

}
# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
pseudo=pseudo_template
gene=gene_template
window_length=window_length_template
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
mutation_density_global_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_groups
mutation_density_dir=$mutation_density_global_dir/$pseudo/$gene
specific_gene_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups"
variant_dir=$global_output_dir/$popu_symbol/variant_two_group

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $specific_gene_global_dir/$pseudo/$gene
mkdir -p $mutation_density_dir
# get_average_variant_count_global_variant "allvar"
# get_average_variant_count_global_variant "snps"
# get_average_variant_count_global_variant "indels"

get_average_variant_count_global_variant_peak_window "allvar"
get_average_variant_count_global_variant_peak_window "snps"
get_average_variant_count_global_variant_peak_window "indels"
