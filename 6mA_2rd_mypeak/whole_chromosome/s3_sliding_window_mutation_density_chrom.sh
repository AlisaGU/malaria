#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
# FUNCTIONS TODO:
get_average_variant_count() {
    var_type=$1

    for por in $(echo $relative_por | tr " " "\n"); do
        mutation_density_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$por/mutation_density
        mkdir -p $mutation_density_dir

        for chrom in $(echo $autosomes | tr " " "\n"); do
            mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_${autosomes}
            rm -rf $mean_variant_count
            chrom_var_count=""

            regions="$chrom_bed_dir/${chrom}.nopeak.bed"

            for w in $(seq 1 $window_number); do
                regions+=" $chrom_bed_dir/${sliding_window_type}/relative.${por}/${chrom}/w$w.bed"
            done
            for region_bed in $(echo $regions | tr " " "\n"); do
                if [[ $(basename $region_bed) =~ "nopeak" ]] || [ $sliding_window_type == "sliding_window" ]; then
                    chrom_var_count+=" $(get_average_variant_count_for_chrom $chrom $region_bed)"
                else
                    chrom_var_count+=" $(get_average_variant_count_for_chrom_mean_window $chrom $region_bed)"
                fi
            done
            echo $chrom_var_count >>$mean_variant_count
        done
    done

}

get_average_variant_count_for_chrom() {
    chrom=$1
    region_bed=$2
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
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

get_average_variant_count_for_chrom_mean_window() {
    chrom=$1
    region_bed=$2
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    all_mean=$(echo $region_bed | sed 's/bed/mean.txt/g')
    rm -rf $all_mean

    while read -r window_bed; do
        b=$(get_one_window_mutation_density_chrom "$window_bed")
        echo $b >>$all_mean
    done <"$region_bed"

    mean_window_mutation_density=$($code_dir/ss3_1_average.r $all_mean)
    echo $mean_window_mutation_density
}

get_one_window_mutation_density_chrom() {
    window_bed=$1
    chrom_region_len=$(echo $window_bed | awk '{print $3-$2}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}
# VARIABLE NAMING TODO:
popu_symbol=popu_symbol_template
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd/whole_chromosome"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"

variant_dir=$global_output_dir/$popu_symbol/variant
chrom_bed_dir=$macs2_out_dir/chrom

relative_por="0.05"
window_number="10"
sliding_window_type="sliding_window_noCollapseWindow"
# sliding_window_type="sliding_window"
# autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
autosomes="autosomes_template"
variant_type="var_type_template"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# mkdir -p $mutation_density_dir
get_average_variant_count "$variant_type"
