#!/bin/bash
#SBATCH -n 4
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_gene_flank() {
    awk '$3=="CDS"' $anno | awk '{split($9,a,";");split(a[1],b,"=");print $1"\t"$4-1"\t"$5"\t"b[2]}' | sort -k1,1 -k2,2n | grep -v -E "Pf3D7_API_v3|Pf_M76611" >$CDS_bed
}

get_mutation_density() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_CDS
    rm -rf $mean_variant_count

    cat $whole_genome_bed_dir/genes.bed | while read line; do
        gene_name=$(echo $line | awk '{print $4}')
        chrom=$(echo $line | awk '{print $1}')

        value=""

        bed=$(grep $gene_name $whole_genome_bed_dir/CDS.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
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

    peak_control_length_in_gene=peak_control_length_in_CDS
    rm -rf $peak_control_length_in_gene

    cat genes.bed | while read line; do
        gene_name=$(echo $line | awk '{print $4}')

        cds_bed_for_gene=$(grep $gene_name CDS.bed)
        peak_exist_in_gene_status=$(echo $"$cds_bed_for_gene" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | wc -l)
        value=""

        CDS_length=$(echo $"$cds_bed_for_gene" | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | awk '{s+=$3-$2}END{print s}')

        CDS_length_in_peak=0
        if [ $peak_exist_in_gene_status -gt 0 ]; then
            CDS_length_in_peak=$(echo $"$cds_bed_for_gene" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | awk '{s+=$3-$2}END{print s}')
        fi
        let CDS_length_in_control=CDS_length-CDS_length_in_peak
        value+="$(echo $gene_name $CDS_length_in_peak $CDS_length_in_control)"

        echo $value >>$peak_control_length_in_gene
    done
}
export -f get_mutation_density get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"

export whole_genome_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export mutation_density_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/wholeGenome
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_with_fold_enrichment=$m6mA_bam_dir/peak_with_fold_enrichment.bed
export anno="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/PlasmoDB-36_Pfalciparum3D7.gff"

export gene_bed=$whole_genome_bed_dir/genes.bed
export CDS_bed=$whole_genome_bed_dir/CDS.bed
export peak_loc=$m6mA_bam_dir/peak.bed
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"

genome_length=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length

parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
$parallel -j 2 get_mutation_density ::: "snps" "indels"
# get_peak_proportion
