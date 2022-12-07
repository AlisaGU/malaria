#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_mutation_density_peak() {
    ## 计算长度是6bp,9bp,12bp和15bp的合并后的motif片段的突变密度
    var_type=$1
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_peak
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        $bedtools merge -i ${chrom}.motif.peak.bed >${chrom}.motif.peak.merged.bed

        value=""
        for merged_motif_len in $(echo $merged_motif_lens | tr " " "\n"); do
            bed=$(awk -v merged_motif_len=$merged_motif_len '$3-$2==merged_motif_len{print $0}' ${chrom}.motif.peak.merged.bed)
            value+=" $(get_average_variant_count_for_chrom $chrom $var_type "$bed")"
        done
        echo $value >>$mean_variant_count
    done
}

compute_mutation_density_all_len_peak() {
    ## 计算所有可能出现的合并后的motif片段的突变密度
    var_type=$1
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_all_possible_peak
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        all_possible_merged_motif_lens=$(awk '{print $3-$2}' ${chrom}.motif.peak.merged.bed | sort -n | uniq)
        for merged_motif_len in $(echo $all_possible_merged_motif_lens | tr " " "\n"); do
            bed=$(awk -v merged_motif_len=$merged_motif_len '$3-$2==merged_motif_len{print $0}' ${chrom}.motif.peak.merged.bed)
            value="$(get_average_variant_count_for_chrom $chrom $var_type "$bed")"
            echo $chrom" "$merged_motif_len" "$value >>$mean_variant_count
        done
    done
}

compute_mutation_density_all_len_nopeak() {
    ## 计算所有可能出现的合并后的motif片段的突变密度
    var_type=$1
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_all_possible_nopeak
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        $bedtools merge -i ${chrom}.motif.nopeak.bed >${chrom}.motif.nopeak.merged.bed

        all_possible_merged_motif_lens=$(awk '{print $3-$2}' ${chrom}.motif.nopeak.merged.bed | sort -n | uniq)
        for merged_motif_len in $(echo $all_possible_merged_motif_lens | tr " " "\n"); do
            bed=$(awk -v merged_motif_len=$merged_motif_len '$3-$2==merged_motif_len{print $0}' ${chrom}.motif.nopeak.merged.bed)
            value="$(get_average_variant_count_for_chrom $chrom $var_type "$bed")"
            echo $chrom" "$merged_motif_len" "$value >>$mean_variant_count
        done
    done
}

compute_mutation_density_each_motif_peak() {
    ## 计算单个合并后的peak motif片段的突变情况
    var_type=$1
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_each_motif_peak
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        cat ${chrom}.motif.peak.merged.bed | while read line; do
            value="$(get_average_variant_count_for_chrom $chrom $var_type "$line")"
            echo $chrom" "$value >>$mean_variant_count
        done
    done
}

compute_mutation_density_each_motif_nopeak() {
    ## 计算单个合并后的nonpeak motif片段的突变情况
    var_type=$1
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_each_motif_nopeak
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        cat ${chrom}.motif.nopeak.merged.bed | while read line; do
            value="$(get_average_variant_count_for_chrom $chrom $var_type "$line")"
            echo $chrom" "$value >>$mean_variant_count
        done
    done
}

get_average_variant_count_for_chrom() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{print $3-$2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/ESEA_WSEA_OCE_SAM_SAS_PASS_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/ESEA_WSEA_OCE_SAM_SAS_PASS_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        # In this case, the ":" is just a "do-nothing" place holder
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}
# VARIABLE NAMING TODO:
seperate_motif_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/motif_bed"

mutation_density_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/different_length_overlap_motif"
variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/variant_two_group"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
merged_motif_lens="6 9 12 15"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $seperate_motif_bed_dir

## 只考虑6bp,9bp,12bp,15bp
compute_mutation_density_peak "snps"
compute_mutation_density_peak "indels"
## 计算所有可能出现的合并后的peak motif片段的合并的突变密度
compute_mutation_density_all_len_peak "snps"
compute_mutation_density_all_len_peak "indels"
## 计算所有可能出现的合并后的nonpeak motif片段的合并的突变密度
compute_mutation_density_all_len_nopeak "snps"
compute_mutation_density_all_len_nopeak "indels"
## 计算所有可能出现的合并后的单个peak motif片段的突变密度
compute_mutation_density_each_motif_peak "snps"
compute_mutation_density_each_motif_peak "indels"
## 计算所有可能出现的合并后的单个nonpeak motif片段的突变密度
compute_mutation_density_each_motif_nopeak "snps"
compute_mutation_density_each_motif_nopeak "indels"
