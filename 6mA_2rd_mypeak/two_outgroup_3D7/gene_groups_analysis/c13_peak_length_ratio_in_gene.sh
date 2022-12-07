#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_peak_control_length_in_gene_for_each_gene_set_include_nopeak_gene_with_fold_enrichment() {
    gene_set=$1
    cd $gene_sets_global_dir/$gene_set/

    mv ${gene_set}_peak_control_length_in_gene ${gene_set}_peak_control_length_in_exon
    peak_control_length_in_gene=${gene_set}_peak_control_length_in_gene
    rm -rf $peak_control_length_in_gene

    cat $gene_sets_global_dir/$gene_set/${gene_set}_genes.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | wc -l)
        value=""

        gene_name=$(echo $line | awk '{print $4}')
        gene_length=$(grep $gene_name ${gene_set}_genes.bed | awk '{s+=$3-$2}END{print s}')

        gene_length_in_peak=0
        if [ $peak_exist_in_gene_status -gt 0 ]; then
            gene_length_in_peak=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore | awk '{s+=$3-$2}END{print s}')
        fi
        let gene_length_in_control=gene_length-gene_length_in_peak
        value+="$(echo $gene_name $gene_length_in_peak $gene_length_in_control)"

        echo $value >>$peak_control_length_in_gene
    done
}

export -f get_peak_control_length_in_gene_for_each_gene_set_include_nopeak_gene_with_fold_enrichment
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_loc=$m6mA_bam_dir/peak.bed
export peak_with_fold_enrichment=$m6mA_bam_dir/peak_with_fold_enrichment.bed
export mutation_density_global_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_5_sets
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group

export gene_sets_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12"
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"
export anno="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/PlasmoDB-36_Pfalciparum3D7.gff"
gene_sets="DR HDR RNA_translation STEVOR VAR RIF" ## 跟c11_gene_5_sets相比，只是加了一个RIF家族
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
$parallel -j 5 get_peak_control_length_in_gene_for_each_gene_set_include_nopeak_gene_with_fold_enrichment ::: $gene_sets
