#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_density_count() {
    # compute_density_count_each_variant "snps"
    # compute_density_count_each_variant "single_del"

    # compute_density_count_each_variant_include_indel "snps"
    # compute_density_count_each_variant_include_indel "indels"

    # compute_density_count_each_variant_noFilter_include_indel "snps"
    # compute_density_count_each_variant_noFilter_include_indel "indels"

    # compute_density_count_each_variant_noFilter_include_indel_bedFirst "snps"
    # compute_density_count_each_variant_noFilter_include_indel_bedFirst "indels"

    # compute_density_count_each_variant_noFilter_include_indel_bedFirst_snp "snps"
    # compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel "indels"

    # compute_density_count_each_variant_noFilter_include_indel_bedFirst_snp "snps"
    # compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap "indels"

    # compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap_both "indels"

    # compute_density_count_indel_LI "indels"
    # compute_density_count_indel_LI_1 "indels"
    compute_density_count_indel_LI_2 "indels"

}

compute_density_count_each_variant() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/$chrom.two_consistent_region.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/peak_nopeak/$chrom.two_consistent_region.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count $chrom $global_dir/peak_nopeak/$chrom.two_consistent_region.nopeak.bed $var_type)
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
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
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
            sort -k1,1 -k2,2n $global_dir/$chrom.two_consistent_region_include_indel.bed | uniq | $bedtools intersect -a stdin -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_include_indel $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.two_consistent_region_include_indel.bed -b $m6mA_peak_bed >$global_dir/$chrom.two_consistent_region_include_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.two_consistent_region_include_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_include_indel $chrom $global_dir/$chrom.two_consistent_region_include_indel.nopeak.bed $var_type)
        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.two_consistent_region_include_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_include_indel() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3.include_indel_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            sort -k1,1 -k2,2n $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed | uniq | $bedtools intersect -a stdin -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed -b $m6mA_peak_bed >$global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel $chrom $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed $var_type)
        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel_bedFirst() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_index
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            sort -k1,1 -k2,2n $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed | uniq | $bedtools intersect -a stdin -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst $chrom $global_dir/$chrom.w${window}.bed $var_type)
            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed -b $m6mA_peak_bed >$global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst $chrom $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed $var_type)
        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel_bedFirst_snp() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_snps
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            sort -k1,1 -k2,2n $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed | uniq | $bedtools intersect -a stdin -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_snp $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/$chrom.noFilter.two_consistent_region_include_indel.bed -b $m6mA_peak_bed >$global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_snp $chrom $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_snp() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_indels
    rm -rf $mean_variant_count
    indel_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/P_falciparum_only/P_reichenowi"

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            # $bedtools intersect -a $indel_bed_dir/${chrom}.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed
            $bedtools intersect -a $indel_bed_dir/${chrom}.bed -b $sliding_window_dir/$chrom/w${window}.bed | $bedtools subtract -a stdin -b /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/prank_speciesB_variants/${chrom}.noBill.3d7.region.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $indel_bed_dir/${chrom}.bed -b $m6mA_peak_bed | $bedtools subtract -a stdin -b /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/prank_speciesB_variants/${chrom}.noBill.3d7.region.bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.two_consistent_region_noFilter_include_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel() {
    chrom=$1
    region_bed=$2
    varType=$3

    # chrom_mutation_count=$($bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/prank_speciesB_variants/${chrom}.scaled_genome_pos.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    chrom_mutation_count=$($bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/prank_speciesB_variants/${chrom}.scaled_genome_pos.except_no_bill_block.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_twoOutgroup_nongap_indels
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $m6mA_peak_bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.only_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.allvariant.noFilterBlock.include_indel.bedFirst.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_each_variant_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap_both() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_twoOutgroup_nongap_indels_both_just_del
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap_both $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $m6mA_peak_bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap_both $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.only_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count_noFilter_include_indel_bedFirst_indel_twoOutgroup_nongap_both() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.twoOutgroup_nogap_only_indel.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_indel_LI() {
    ##从block角度出发，选择有外群的block
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_LI
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $variant_dir/${chrom}.include_P_billcollinsi.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_mutation_density_count_indel_LI $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $variant_dir/${chrom}.include_P_billcollinsi.bed -b $m6mA_peak_bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_mutation_density_count_indel_LI $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.only_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_mutation_density_count_indel_LI() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.include_P_billcollinsi.intermediate_only_indel_polysite.scaled_genome_pos.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_density_count_indel_LI_1() {
    ##相比于block角度，进一步深化，挑选外群出现的位置。舍弃block中有一部分没有外群的地方
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_LI_1
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $variant_dir/$chrom.outgroup_region.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_mutation_density_count_indel_LI $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $variant_dir/$chrom.outgroup_region.bed -b $m6mA_peak_bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_mutation_density_count_indel_LI $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.only_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

compute_density_count_indel_LI_2() {
    ## 相比于compute_density_count_indel_LI_1，排除标准株是gap的区间和变异
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_sliding_window_mean_variant_count_LI_2
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for window in $(seq 1 10); do
            $bedtools intersect -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $sliding_window_dir/$chrom/w${window}.bed >$global_dir/$chrom.w${window}.bed

            window_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.w${window}.bed)
            window_mutation_count=$(get_mutation_density_count_indel_LI_2 $chrom $global_dir/$chrom.w${window}.bed $var_type)

            value+=" $window_mutation_count $window_region_len"
            rm -rf $global_dir/$chrom.w${window}.bed
        done
        $bedtools subtract -a $global_dir/${chrom}.noFilter.3D7_outgroupA_both_exist.bed -b $m6mA_peak_bed >$global_dir/$chrom.only_indel.nopeak.bed
        nopeak_region_len=$(awk '{s+=$3-$2}END{print s}' $global_dir/$chrom.only_indel.nopeak.bed)
        nopeak_mutation_count=$(get_mutation_density_count_indel_LI_2 $chrom $global_dir/$chrom.only_indel.nopeak.bed $var_type)

        value+=" $nopeak_mutation_count $nopeak_region_len"
        rm -rf $global_dir/$chrom.only_indel.nopeak.bed
        echo $value >>$mean_variant_count
    done
}

get_mutation_density_count_indel_LI_2() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.include_P_billcollinsi.intermediate_only_indel_polysite_standardStrainBaseExist.scaled_genome_pos.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}
# VARIABLE NAMING TODO:
# global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_reichenowi"
global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi"

sliding_window_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/sliding_window/relative.0.05"

mutation_density_dir=$global_dir/mutation_density
variant_dir=$global_dir/../../maf2vcf
m6mA_peak_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
compute_density_count
