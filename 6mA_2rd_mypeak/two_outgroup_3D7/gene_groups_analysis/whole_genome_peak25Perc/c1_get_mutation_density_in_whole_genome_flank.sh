#!/bin/bash
#SBATCH -n 24
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_gene_flank() {
    $bedtools slop -i $gene_bed -g $genome_length -b 2000 >$gene_flank2kb_bed
    $bedtools slop -i $gene_bed -g $genome_length -b 1000 >$gene_flank1kb_bed

}

get_mutation_density() {
    var_type=$1
    flank_length=$2

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_flank${flank_length}
    rm -rf $mean_variant_count

    cat $whole_genome_bed_dir/gene_flank${flank_length}.bed | while read line; do
        gene_name=$(echo $line | awk '{print $4}')
        chrom=$(echo $line | awk '{print $1}')

        value=""

        result=$(get_average_variant_count_for_chrom $chrom $var_type "$line")
        value+="$(echo $gene_name "$result")"

        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{s+=$3-$2}END{print s}')

    chrom_mutation_count=""
    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_correct_spanning_deletion_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_correct_spanning_deletion_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}

get_peak_proportion() {
    flank_length=$1

    cd $whole_genome_bed_dir

    peak_control_length_in_gene=peak_control_length_in_flank${flank_length}
    rm -rf $peak_control_length_in_gene

    cat gene_flank${flank_length}.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc -wb | wc -l)
        value=""

        gene_name=$(echo $line | awk '{print $4}')
        gene_flank_length=$(echo $line | awk '{s+=$3-$2}END{print s}')

        gene_flank_length_in_peak=0
        if [ $peak_exist_in_gene_status -gt 0 ]; then
            gene_flank_length_in_peak=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | awk '{s+=$3-$2}END{print s}')
            gene_flank_peak_count=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc -wb | awk '{print $9}' | sort | uniq | wc -l)
        fi
        let gene_flank_length_in_control=gene_flank_length-gene_flank_length_in_peak
        value+="$(echo $gene_name $gene_flank_length_in_peak $gene_flank_length_in_control $gene_flank_peak_count)"

        echo $value >>$peak_control_length_in_gene
    done
}

export -f get_mutation_density get_average_variant_count_for_chrom get_peak_proportion
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"

export whole_genome_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export mutation_density_dir=$global_output_dir/${popu_symbol}/wholeGenome_peak25Perc/3D7
export variant_dir=$global_output_dir/$popu_symbol/variant_3D7 ## 变异信息是以3D7为祖先型的
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"

export gene_bed=$whole_genome_bed_dir/genes.bed
export gene_flank2kb_bed=$whole_genome_bed_dir/gene_flank2kb.bed
export gene_flank1kb_bed=$whole_genome_bed_dir/gene_flank1kb.bed

export peak_loc=$m6mA_bam_dir/chrom_wholeGenome/sliding_window/relative.0.25/sliding_window_25Perc.bed

export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex

$parallel -j 6 get_mutation_density ::: "snps" "indels" ::: "0kb" "1kb" "2kb"
