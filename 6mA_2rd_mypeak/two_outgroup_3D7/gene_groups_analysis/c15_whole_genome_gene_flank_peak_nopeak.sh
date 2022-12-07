#!/bin/bash
#SBATCH -n 8
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment() {

    var_type=$1

    mutation_density_dir=$mutation_density_global_dir
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_peak_nopeak_flank
    rm -rf $mean_variant_count

    cat $gene_sets_global_dir/gene_flank.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $peak_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA")"

            chrom=$(echo $line | awk '{print $1}')

            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result" "NA")"
        elif [ $peak_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            #peak region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            mean_fold_enrichment=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment -wb | awk '{print $NF}' | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
            value+="$(echo $gene_name "$result")"
            #control region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")

            value+=" $(echo "$result" $mean_fold_enrichment)"
        fi
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{print $3-$2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        # In this case, the ":" is just a "do-nothing" place holder
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}
export -f get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_loc=$m6mA_bam_dir/peak.bed
export peak_with_fold_enrichment=$m6mA_bam_dir/peak_with_fold_enrichment.bed
export mutation_density_global_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/wholeGenome
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group

export gene_sets_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"
export anno="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/PlasmoDB-36_Pfalciparum3D7.gff"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
$parallel -j 2 get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment ::: "snps" "indels"
