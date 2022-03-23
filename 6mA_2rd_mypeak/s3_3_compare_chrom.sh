#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
# FUNCTIONS TODO:
get_average_variant_count() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        regions="$chrom_bed_dir/${chrom}.nopeak.bed"
        for por in $(echo $relative_por | tr " " "\n"); do
            regions+=" $chrom_bed_dir/relative.${por}/${chrom}.peak.relative.${por}.bed"
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            regions+=" $chrom_bed_dir/absolute.${len}/${chrom}.peak.absolute.${len}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom $chrom $region_bed)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done
    # $code_dir/s6_plot.R $mean_variant_count
}

get_average_variant_count_for_chrom() {
    chrom=$1
    region_bed=$2
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_mutation_count_bed=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_snp_polysite.bed
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_indel_polysite.bed
    fi

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
    chrom_mutation_count=$($bedtools intersect -a $chrom_mutation_count_bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    chrom_region_mutation_count_mean=$(echo "scale=4;${chrom_mutation_count}/${chrom_region_len}" | bc)
    echo $chrom_region_mutation_count_mean
}

# VARIABLE NAMING TODO:
popu_symbol=$1
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"

variant_dir=$global_output_dir/$popu_symbol/variant
mutation_density_dir=$global_output_dir/$popu_symbol/mutation_density
chrom_bed_dir=$macs2_out_dir/chrom

absolute_len="10 30 50 100 200"
relative_por="0.01 0.05 0.1 0.2 0.5"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
mkdir -p $mutation_density_dir
get_average_variant_count "allvar"
get_average_variant_count "snps"
get_average_variant_count "indels"
