#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
# FUNCTIONS TODO:

get_average_variant_count_3control_N2N_indel() {
    var_type=$1

    # mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_3control_count
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_3control

    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""
        chrom_var_count="$(get_average_variant_count_for_chrom_base_count $chrom $macs2_two_outgroup_consistent_dir/${chrom}.nopeak.bed $var_type)"
        chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_base_count $chrom $macs2_two_outgroup_consistent_dir/${chrom}.nopeak.bed $var_type)"
        chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_base_count $chrom $macs2_two_outgroup_consistent_dir/$chrom.remove_margin_500_nopeak.bed $var_type)"

        regions=""
        for por in $(echo $regions_type | tr " " "\n"); do
            regions+=" $peak_motif_region_loc_dir/${chrom}_${por}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom_base_count $chrom $region_bed $var_type)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done

}

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    var_type=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    ori_base=$(echo $var_type | awk -F"2|_" '{print $1}')
    chrom_region_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${ori_base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1

    # echo $chrom_mutation_count
}

get_average_variant_count_for_chrom_only_motif_base_count() {
    chrom=$1
    region_bed=$2
    var_type=$3

    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    ori_base=$(echo $var_type | awk -F"2|_" '{print $1}')

    motif_dir=""
    if [ $(basename $region_bed | grep "nopeak" | wc -l) -eq 1 ]; then
        motif_dir=$(echo $nopeak_motif_bed_dir)
    else
        motif_dir=$(echo $peak_motif_dir)
    fi

    awk '{print $0"\t""motif_"NR}' ${motif_dir}/${chrom}.*.bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' | $bedtools intersect -wa -a stdin -b $region_bed >$motif_dir/${chrom}.${popu_symbol}.intermediate

    chrom_region_len=$($seqkit subseq --quiet --bed $motif_dir/${chrom}.${popu_symbol}.intermediate $genome | $seqkit fx2tab -i -H -n -C "${ori_base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
    # echo $chrom_mutation_count

    rm -rf $motif_dir/${chrom}.${popu_symbol}.intermediate
}

get_average_variant_count_3control_global_variant() {
    var_type=$1

    # mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_3control_count
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_3control

    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""
        chrom_var_count="$(get_average_variant_count_for_chrom_whole_region $chrom $macs2_two_outgroup_consistent_dir/${chrom}.nopeak.bed $var_type)"
        chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_whole_region $chrom $macs2_two_outgroup_consistent_dir/${chrom}.nopeak.bed $var_type)"
        chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_whole_region $chrom $macs2_two_outgroup_consistent_dir/$chrom.remove_margin_500_nopeak.bed $var_type)"

        regions=""
        for por in $(echo $regions_type | tr " " "\n"); do
            regions+=" $peak_motif_region_loc_dir/${chrom}_${por}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom_whole_region $chrom $region_bed $var_type)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_whole_region() {
    chrom=$1
    region_bed=$2
    var_type=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1

    # echo $chrom_mutation_count
}

get_average_variant_count_for_chrom_only_motif_whole_region() {
    chrom=$1
    region_bed=$2
    var_type=$3

    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    ori_base=$(echo $var_type | awk -F"2|_" '{print $1}')

    motif_dir=""
    if [ $(basename $region_bed | grep "nopeak" | wc -l) -eq 1 ]; then
        motif_dir=$(echo $nopeak_motif_bed_dir)
    else
        motif_dir=$(echo $peak_motif_dir)
    fi

    awk '{print $0"\t""motif_"NR}' ${motif_dir}/${chrom}.*.bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' | $bedtools intersect -wa -a stdin -b $region_bed >$motif_dir/${chrom}.${popu_symbol}.intermediate

    chrom_region_len=$(awk '{print $3-$2}' $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.${popu_symbol}.intermediate | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
    # echo $chrom_mutation_count

    rm -rf $motif_dir/${chrom}.${popu_symbol}.intermediate
}
# VARIABLE NAMING TODO:
popu_symbol=$1
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
variant_dir=$global_output_dir/$popu_symbol/variant_two_group
mutation_density_dir=$global_output_dir"/${popu_symbol}/motif/mutation_dentisy_two_outgroup_consistent"

motif_id="2"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
nopeak_motif_bed_dir=$macs2_two_outgroup_consistent_dir/motif_pattern_${motif_id}_nopeak/motif_loc
peak_motif_region_loc_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_${motif_id}/region_loc"

regions_type="collapsed_motif region1 region2 region3 region4"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
mkdir -p $mutation_density_dir
get_average_variant_count_3control_global_variant "allvar"
get_average_variant_count_3control_global_variant "snps"
get_average_variant_count_3control_global_variant "indels"
# get_average_variant_count_3control_N2N_indel "A2T"
# get_average_variant_count_3control_N2N_indel "A2G"
# get_average_variant_count_3control_N2N_indel "A2C"
# get_average_variant_count_3control_N2N_indel "T2A"
# get_average_variant_count_3control_N2N_indel "T2G"
# get_average_variant_count_3control_N2N_indel "T2C"
# get_average_variant_count_3control_N2N_indel "C2T"
# get_average_variant_count_3control_N2N_indel "C2G"
# get_average_variant_count_3control_N2N_indel "C2A"
# get_average_variant_count_3control_N2N_indel "G2A"
# get_average_variant_count_3control_N2N_indel "G2C"
# get_average_variant_count_3control_N2N_indel "G2T"
# get_average_variant_count_3control_N2N_indel "A_single_del"
# get_average_variant_count_3control_N2N_indel "T_single_del"
# get_average_variant_count_3control_N2N_indel "C_single_del"
# get_average_variant_count_3control_N2N_indel "G_single_del"

# get_average_variant_count_3control_N2N_indel "A_indel"
# get_average_variant_count_3control_N2N_indel "T_indel"
# get_average_variant_count_3control_N2N_indel "C_indel"
# get_average_variant_count_3control_N2N_indel "G_indel"
