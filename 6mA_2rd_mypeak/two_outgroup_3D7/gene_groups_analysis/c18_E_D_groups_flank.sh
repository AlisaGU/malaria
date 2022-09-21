#!/bin/bash
#SBATCH -n 4
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_mutation_density() {
    var_type=$1
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_flank
    rm -rf $mean_variant_count

    cat $whole_genome_bed_dir/gene_flank.bed | while read line; do
        gene_name=$(echo $line | awk '{print $4}')
        chrom=$(echo $line | awk '{print $1}')

        value=""

        bed=$(echo $"$line" | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
        result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
        value+="$(echo $gene_name "$result")"

        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{s+=$3-$2}END{print s}')

    chrom_mutation_count=""
    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}

get_peak_proportion() {
    cd $whole_genome_bed_dir

    peak_control_length_in_gene=peak_control_length_in_flank
    rm -rf $peak_control_length_in_gene

    cat gene_flank.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | wc -l)
        value=""

        gene_name=$(echo $line | awk '{print $4}')
        gene_flank_length=$(grep $gene_name gene_flank.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | awk '{s+=$3-$2}END{print s}')

        gene_flank_length_in_peak=0
        if [ $peak_exist_in_gene_status -gt 0 ]; then
            gene_flank_length_in_peak=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | awk '{s+=$3-$2}END{print s}')
        fi
        let gene_flank_length_in_control=gene_flank_length-gene_flank_length_in_peak
        value+="$(echo $gene_name $gene_flank_length_in_peak $gene_flank_length_in_control)"

        echo $value >>$peak_control_length_in_gene
    done
}
export -f get_mutation_density get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"

export whole_genome_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/E_D_group"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export mutation_density_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/E_D_group
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_with_fold_enrichment=$m6mA_bam_dir/peak_with_fold_enrichment.bed

export gene_flank_bed=$whole_genome_bed_dir/gene_flank.bed
export peak_loc=$m6mA_bam_dir/peak.bed
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"

export genome_length=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length

export parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
$parallel -j 2 get_mutation_density ::: "snps" "indels"
get_peak_proportion
