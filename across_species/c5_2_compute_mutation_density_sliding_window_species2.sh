#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_density_count() {
    # compute_density_count_each_variant "snps"
    # compute_density_count_each_variant "indels"

    compute_density_count_each_variant_include_indel "snps"
    compute_density_count_each_variant_include_indel "indels"
}

compute_density_count_each_variant() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/$chrom.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.bed -b $m6mA_peak_bed >$global_dir/$chrom.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count $chrom $global_dir/$chrom.nopeak.bed $var_type)
        value+=" $nopeak_mutation_count $nopeak_region_len"
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_include_indel() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_include_indel
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/$chrom.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_include_indel $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.bed -b $m6mA_peak_bed >$global_dir/$chrom.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_include_indel $chrom $global_dir/$chrom.nopeak.bed $var_type)
        value+=" $nopeak_mutation_count $nopeak_region_len"
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_include_indel() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species2.include_indel_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}
# VARIABLE NAMING TODO:
sliding_window_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/sliding_window/relative.0.05"

m6mA_peak_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex

# for species2 in $(echo "P_praefalciparum P_reichenowi" | tr " " "\n"); do
for species2 in $(echo "P_reichenowi" | tr " " "\n"); do
    global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/P_falciparum_only/$species2"
    mutation_density_dir=$global_dir/mutation_density
    mkdir -p $mutation_density_dir
    variant_dir=$global_dir/../../maf2vcf/P_falciparum.$species2
    compute_density_count
done
