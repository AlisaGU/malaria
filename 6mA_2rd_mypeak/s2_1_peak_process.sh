#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
# FUNCTIONS TODO:
get_chrom_treatment_control_level_bed() {
    cd $modifi_level_dir

    awk -F "\t" 'NR>1{print>$1}' $treatment_level
    for file in Pf3D7_*_v3; do
        mv $file $file".treat.bed"
    done
    mv Pf_M76611 Pf_M76611.treat.bed

    awk -F "\t" 'NR>1{print>$1}' $control_level
    for file in Pf3D7_*_v3; do
        mv $file $file".cont.bed"
    done
    mv Pf_M76611 Pf_M76611.cont.bed
    cd -
}

select_chrom_peak_nopeak_bed() {
    cd $chrom_bed_dir

    awk -F "\t" '{print>$1}' $core_peak_bed
    for file in Pf3D7_*_v3; do
        mv $file $file".peak.bed"
    done

    awk -F "\t" '{print>$1}' $core_nopeak_bed
    for file in Pf3D7_*_v3; do
        mv $file $file".nopeak.bed"
    done

    cd -
}

get_peak_submit_flank_region_pos() {
    # for absolute in $(echo $absolute_len | tr " " "\n"); do
    #     for chrom in $(echo $autosomes | tr " " "\n"); do
    #         get_peak_submit_flank_region_pos_for_chrom $chrom "absolute" $absolute
    #     done
    # done

    for relative in $(echo $relative_por | tr " " "\n"); do
        for chrom in $(echo $autosomes | tr " " "\n"); do
            get_peak_submit_flank_region_pos_for_chrom $chrom "relative" $relative
        done
    done
}

get_peak_submit_flank_region_pos_for_chrom() {
    chrom=$1
    flank_length_type=$2
    flank_length=$3
    output_dir=$chrom_bed_dir/${flank_length_type}.${flank_length}
    mkdir -p $output_dir

    peak_pos=$chrom_bed_dir/${chrom}.peak.bed
    treatment_level=$modifi_level_dir/${chrom}.treat.bed
    control_level=$modifi_level_dir/${chrom}.cont.bed
    output=$output_dir/${chrom}.peak.${flank_length_type}.${flank_length}.bed

    $code_dir/s6_get_peak_submit_flank_pos.R $chrom $peak_pos $flank_length_type $flank_length $treatment_level $control_level $output
}

get_average_variant_count() {
    var_type=$1

    mean_variant_count=$chrom_bed_dir/${var_type}_mean_variant_count
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
        chrom_mutation_count_bed=$variant_dir/WSEA_PASS_${chrom_number}_polysite.bed
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count_bed=$variant_dir/WSEA_PASS_${chrom_number}_snp_polysite.bed
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count_bed=$variant_dir/WSEA_PASS_${chrom_number}_indel_polysite.bed
    fi

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
    chrom_mutation_count=$($bedtools intersect -a $chrom_mutation_count_bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    chrom_region_mutation_count_mean=$(echo "scale=4;${chrom_mutation_count}/${chrom_region_len}" | bc)
    echo $chrom_region_mutation_count_mean
}

plot_mutation_density() {
    $code_dir/s6_plot.R
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
ref_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7"
chrom_bed_dir=$macs2_out_dir/chrom
modifi_level_dir=$macs2_out_dir/modifi_level

core_bed=$ref_dir/anno/3D7.core.bed
peak_bed=$macs2_out_dir/"peak.bed"
core_peak_bed=$macs2_out_dir/"peak.core.bed"
nopeak_bed=$macs2_out_dir/"nopeak.bed"
core_nopeak_bed=$macs2_out_dir/"nopeak.core.bed"
chrom_len=$ref_dir/"genome/chrom.length"
treatment_level=$macs2_out_dir/"3D7_2rd_treat_pileup.bdg"
control_level=$macs2_out_dir/"3D7_2rd_control_lambda.bdg"
absolute_len="10 30 50 100 200"
relative_por="0.01 0.05 0.1 0.2 0.5"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
$bedtools complement -i $peak_bed -g $chrom_len -L >$nopeak_bed # 非peak区位置
$bedtools intersect -a $peak_bed -b $core_bed >$core_peak_bed
$bedtools intersect -a $nopeak_bed -b $core_bed >$core_nopeak_bed
select_chrom_peak_nopeak_bed
get_chrom_treatment_control_level_bed
get_peak_submit_flank_region_pos
