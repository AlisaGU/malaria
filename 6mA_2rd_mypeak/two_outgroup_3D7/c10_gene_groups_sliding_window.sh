#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_all_chroms_bam_depth() {
    for autosome in $(echo $autosomes | tr " " "\n"); do
        for types in $(echo "ChIP Input" | tr " " "\n"); do
            sed -e "s/chrom_template/${autosome}/g" -e "s/types_template/${types}/g" $code_dir/s10_get_chrom_bam_depth.sh >$code_dir/s10_get_chrom_bam_depth.sh_${autosome}_${types}
            sbatch $code_dir/s10_get_chrom_bam_depth.sh_${autosome}_${types}
        done
    done
}

get_50_nonoverlap_windows_bed_for_gene_groups() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        cd $pseudo
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            cd $gene
            $bedtools makewindows -b ${gene}_genes.bed -n 50 -i srcwinnum >${gene}_genes.window.bed
            cd ../
        done
        cd ../
    done
}

get_window_6mA_fold_change() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            sed "s/pseudo_template/$pseudo/g" $code_dir/s10_gene_groups_sliding_window_6mA.sh | sed "s/gene_template/$gene/g" >$code_dir/s10_gene_groups_sliding_window_6mA.sh_${pseudo}_${gene}
            sbatch $code_dir/s10_gene_groups_sliding_window_6mA.sh_${pseudo}_${gene}
        done
    done
}

get_window_mutation_density() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            sed "s/pseudo_template/$pseudo/g" $code_dir/s10_gene_groups_sliding_window_mutation_density.sh | sed "s/popu_symbol_template/$popu_symbol/g" | sed "s/gene_template/$gene/g" >$code_dir/s10_gene_groups_sliding_window_mutation_density.sh_${popu_symbol}_${pseudo}_${gene}
            sbatch $code_dir/s10_gene_groups_sliding_window_mutation_density.sh_${popu_symbol}_${pseudo}_${gene}
        done
    done
}

get_fixed_length_window_for_intersection_between_gene_and_peak() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        cd $pseudo
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            cd $gene
            $bedtools intersect -a ${gene}_genes.bed -b $peak_loc | $bedtools makewindows -b stdin -w $window_length -s $window_length -i srcwinnum | awk '{print $0"\t"NR}' | $bedtools subtract -a stdin -b $two_outgroup_inconsistent_region_in_core_and_noncore >${gene}_genes.peak.window.${window_length}bp.bed
            cd ../
        done
        cd ../
    done
}

get_peak_window_6mA_fold_change() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            sed "s/pseudo_template/$pseudo/g" $code_dir/s10_gene_groups_sliding_window_6mA.sh | sed "s/gene_template/$gene/g" | sed "s/window_length_template/$window_length/g" >$code_dir/s10_gene_groups_sliding_window_6mA.sh_${pseudo}_${gene}
            sbatch $code_dir/s10_gene_groups_sliding_window_6mA.sh_${pseudo}_${gene}
        done
    done
}

get_peak_window_mutation_density() {
    for pseudo in $(echo "exclude_pseudo include_pseudo" | tr " " "\n"); do
        for gene in $(echo $gene_group_type | tr " " "\n"); do
            sed "s/pseudo_template/$pseudo/g" $code_dir/s10_gene_groups_sliding_window_mutation_density.sh | sed "s/popu_symbol_template/$popu_symbol/g" | sed "s/gene_template/$gene/g" | sed "s/window_length_template/$window_length/g" >$code_dir/s10_gene_groups_sliding_window_mutation_density.sh_${popu_symbol}_${pseudo}_${gene}
            sbatch $code_dir/s10_gene_groups_sliding_window_mutation_density.sh_${popu_symbol}_${pseudo}_${gene}
        done
    done
}
# VARIABLE NAMING TODO:
window_length=50
popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
mutation_density_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_groups
specific_gene_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups"

two_outgroup_consistent_region_in_core_and_noncore_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore
two_outgroup_inconsistent_region_in_core_and_noncore=$two_outgroup_consistent_region_in_core_and_noncore_dir/inconsistent_region_in_core_and_noncore.bed

m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
treat_bam=$m6mA_bam_dir/3D7_2rd_treat_pileup.bdg
input_bam=$m6mA_bam_dir/3D7_2rd_control_lambda.bdg
peak_loc=$m6mA_bam_dir/peak.bed
gene_group_type="BER HDR MMR  NER invasion rhoptry RIF STEVOR SURF VAR PHIST MC-2TM DnaJ DR"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $specific_gene_global_dir

get_all_chroms_bam_depth
# get_50_nonoverlap_windows_bed_for_gene_groups
# get_window_6mA_fold_change
# get_window_mutation_density

get_fixed_length_window_for_intersection_between_gene_and_peak
get_peak_window_6mA_fold_change
get_peak_window_mutation_density

grep "CDS" VAR_anno | awk '{print $1"\t"$4"\t"$5"\t"$3}' | bedtools intersect -a VAR_genes.peak.window.50bp.bed -b stdin -wa -wb -loj | sed 's/\./noncoding/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9}' >VAR_genes.peak.window.50bp.anno.bed
