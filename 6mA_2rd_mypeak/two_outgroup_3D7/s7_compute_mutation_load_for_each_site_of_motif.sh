#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_average_variant_count_based_on_motif_each_site_status() {
    site_index=$1
    motif_seq=$2
    base_on_this_site=$(echo ${motif_seq} | grep -o . | awk -v order=${site_index} 'NR==order{print}')

    if [ $base_on_this_site == "A" ]; then
        get_average_variant_count_based_on_site_status "A2G" $site_index
        get_average_variant_count_based_on_site_status "A2C" $site_index
        get_average_variant_count_based_on_site_status "A2T" $site_index
        get_average_variant_count_based_on_site_status "A_single_del" $site_index
    elif [ $base_on_this_site == "T" ]; then
        get_average_variant_count_based_on_site_status "T2C" $site_index
        get_average_variant_count_based_on_site_status "T2A" $site_index
        get_average_variant_count_based_on_site_status "T2G" $site_index
        get_average_variant_count_based_on_site_status "T_single_del" $site_index
    elif [ $base_on_this_site == "C" ]; then
        get_average_variant_count_based_on_site_status "C2T" $site_index
        get_average_variant_count_based_on_site_status "C2A" $site_index
        get_average_variant_count_based_on_site_status "C2G" $site_index
        get_average_variant_count_based_on_site_status "C_single_del" $site_index
    elif [ $base_on_this_site == "G" ]; then
        get_average_variant_count_based_on_site_status "G2T" $site_index
        get_average_variant_count_based_on_site_status "G2A" $site_index
        get_average_variant_count_based_on_site_status "G2C" $site_index
        get_average_variant_count_based_on_site_status "G_single_del" $site_index
    fi
}

get_average_variant_count_for_each_motif_global_variants() {
    site_index=$1
    motif_seq=$2
    get_average_variant_count_for_global_variants "allvar"
    get_average_variant_count_for_global_variants "snps"
    get_average_variant_count_for_global_variants "indels"
    get_average_variant_count_for_global_variants "single_del"

}

get_average_variant_count_based_on_site_status() {
    var_type=$1
    site_index=$2

    plus_first_base=$(echo $var_type | awk -F"2|_" '{print $1}')
    plus_second_base=$(echo $var_type | sed "s/_/2/1" | awk -F"2" '{print $2}')
    separator=$(echo $var_type | grep -o . | awk 'NR==2{print}')
    minus_var_type=$(get_minus_var_type $plus_first_base)$separator$(get_minus_var_type $plus_second_base)

    mean_variant_count=$mutation_density_dir/site${site_index}.${var_type}_mean_variant_count

    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        base_count_on_this_site_peak=$(cat $motif_site_dir/site${site_index}/${chrom}.peak.* | wc -l)
        chrom_var_count_peak_plus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.peak.plus.bed $var_type)"
        chrom_var_count_peak_minus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.peak.minus.bed $minus_var_type)"
        let chrom_var_count_peak=chrom_var_count_peak_plus+chrom_var_count_peak_minus

        base_count_on_this_site_nopeak=$(cat $motif_site_dir/site${site_index}/${chrom}.nopeak.* | wc -l)
        chrom_var_count_nopeak_plus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.nopeak.plus.bed $var_type)"
        chrom_var_count_nopeak_minus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.nopeak.minus.bed $minus_var_type)"
        let chrom_var_count_nopeak=chrom_var_count_nopeak_plus+chrom_var_count_nopeak_minus

        echo "$chrom_var_count_peak $base_count_on_this_site_peak $chrom_var_count_nopeak $base_count_on_this_site_nopeak" >>$mean_variant_count
    done

}

get_average_variant_count_for_global_variants() {
    var_type=$1
    mean_variant_count=$mutation_density_dir/site${site_index}.${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        base_count_on_this_site_peak=$(cat $motif_site_dir/site${site_index}/${chrom}.peak.* | wc -l)
        chrom_var_count_peak_plus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.peak.plus.bed $var_type)"
        chrom_var_count_peak_minus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.peak.minus.bed $var_type)"
        let chrom_var_count_peak=chrom_var_count_peak_plus+chrom_var_count_peak_minus

        base_count_on_this_site_nopeak=$(cat $motif_site_dir/site${site_index}/${chrom}.nopeak.* | wc -l)
        chrom_var_count_nopeak_plus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.nopeak.plus.bed $var_type)"
        chrom_var_count_nopeak_minus="$(get_average_variant_count_for_chrom_base_count $chrom $motif_site_dir/site${site_index}/${chrom}.nopeak.minus.bed $var_type)"
        let chrom_var_count_nopeak=chrom_var_count_nopeak_plus+chrom_var_count_nopeak_minus

        echo "$chrom_var_count_peak $base_count_on_this_site_peak $chrom_var_count_nopeak $base_count_on_this_site_nopeak" >>$mean_variant_count
    done
}

get_minus_var_type() {
    base=$1
    if [ $base == "A" ]; then
        echo "T"
    elif [ $base == "T" ]; then
        echo "A"
    elif [ $base == "C" ]; then
        echo "G"
    elif [ $base == "G" ]; then
        echo "C"
    elif [ $base == "single_del" ]; then
        echo "single_del"
    fi
}

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
motif_seq=motif_seq_template

global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
variant_dir=$global_output_dir/$popu_symbol/variant_two_group
mutation_density_dir=$global_output_dir"/${popu_symbol}/motif/mutation_dentisy_two_outgroup_consistent/each_site/$motif_seq"

macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
motif_site_global_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
motif_site_dir=$motif_site_global_dir/$motif_seq

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
mkdir -p $mutation_density_dir
for site_index in $(seq 1 6); do
    # get_average_variant_count_based_on_motif_each_site_status ${site_index} ${motif_seq}
    get_average_variant_count_for_each_motif_global_variants ${site_index} ${motif_seq}
done
