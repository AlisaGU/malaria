#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
split_intergenic_region_by_chrom() {
    mkdir -p $intergenic_dir
    cd $intergenic_dir
    awk '{print >$1}' $intergenic_bed
    cd -
}

get_peak_submit_and_control_in_integenic_region() {
    for autosome in $(echo $autosomes | tr " " "\n"); do
        $bedtools intersect -a $intergenic_dir/$autosome -b $peak_nopeak_bed_dir/$autosome.nopeak.bed >$working_dir/${autosome}.nopeak.intergic.bed
    done

    for absolute in $(echo $absolute_len | tr " " "\n"); do
        specific_peak_flank_dir=$working_dir/"absolute."$absolute
        mkdir -p $specific_peak_flank_dir
        for autosome in $(echo $autosomes | tr " " "\n"); do
            $bedtools intersect -a $intergenic_dir/$autosome -b $peak_nopeak_bed_dir/absolute.$absolute/$autosome.* >$specific_peak_flank_dir/${autosome}.peak.intergic.bed
        done
    done

    for relative in $(echo $relative_por | tr " " "\n"); do
        specific_peak_flank_dir=$working_dir/"relative."$relative
        mkdir -p $specific_peak_flank_dir
        for autosome in $(echo $autosomes | tr " " "\n"); do
            $bedtools intersect -a $intergenic_dir/$autosome -b $peak_nopeak_bed_dir/relative.$relative/$autosome.* >$specific_peak_flank_dir/${autosome}.peak.intergic.bed
        done
    done
}

get_average_variant_count() {
    for popu_symbol in $(echo $popu_symbols); do
        # get_average_variant_count_each_popu $popu_symbol "allvar"
        # get_average_variant_count_each_popu $popu_symbol "snps"
        # get_average_variant_count_each_popu $popu_symbol "indels"
        get_average_variant_count_each_popu $popu_symbol "transition"
        get_average_variant_count_each_popu $popu_symbol "transversion"
        $code_dir/s2_get_ts_tv_ratio_each_popu.R $popu_symbol
        get_average_variant_count_each_popu $popu_symbol "A2G"
        get_average_variant_count_each_popu $popu_symbol "G2A"
        get_average_variant_count_each_popu $popu_symbol "T2C"
        get_average_variant_count_each_popu $popu_symbol "C2T"
        get_average_variant_count_each_popu $popu_symbol "A2C"
        get_average_variant_count_each_popu $popu_symbol "A2T"
        get_average_variant_count_each_popu $popu_symbol "G2C"
        get_average_variant_count_each_popu $popu_symbol "G2T"
        get_average_variant_count_each_popu $popu_symbol "C2A"
        get_average_variant_count_each_popu $popu_symbol "C2G"
        get_average_variant_count_each_popu $popu_symbol "T2A"
        get_average_variant_count_each_popu $popu_symbol "T2G"
        get_average_variant_count_each_popu $popu_symbol "fixation"
        $code_dir/s2_compare_mutation.R $popu_symbol
    done
}

get_average_variant_count_each_popu() {
    popu_symbol=$1
    var_type=$2

    mutation_density_dir=$working_dir/$popu_symbol/mutation_density
    mkdir -p $mutation_density_dir

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        regions="$working_dir/${chrom}.nopeak.intergic.bed"
        for por in $(echo $relative_por | tr " " "\n"); do
            regions+=" $working_dir/relative.${por}/${chrom}.peak.intergic.bed"
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            regions+=" $working_dir/absolute.${len}/${chrom}.peak.intergic.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom $region_bed $popu_symbol)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom() {
    region_bed=$1
    popu_symbol=$2
    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}/variant"

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
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/partition"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_narrow_peak_out=$macs2_out_dir/"3D7_2rd_peaks.narrowPeak"
ref_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7"
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic"
intergenic_dir=$working_dir/intergenic
peak_nopeak_bed_dir=$macs2_out_dir/chrom

chrom_len=$ref_dir/"genome/chrom.length"
ref_anno=$ref_dir/"anno/PlasmoDB-36_Pfalciparum3D7.gff"

peak_bed=$macs2_out_dir/"peak.bed"
nopeak_bed=$macs2_out_dir/"nopeak.bed"
intergenic_bed=$ref_dir/"anno"/intergenic.bed
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
absolute_len="10 30 50 100 200"
relative_por="0.01 0.05 0.1 0.2 0.5"
# popu_symbols="CAF EAF WAF WAF_CAF_EAF WSEA SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE"
popu_symbols="CAF"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
mkdir -p $working_dir
split_intergenic_region_by_chrom
get_peak_submit_and_control_in_integenic_region
get_average_variant_count
