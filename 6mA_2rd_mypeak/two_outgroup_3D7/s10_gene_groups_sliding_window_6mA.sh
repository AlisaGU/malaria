#!/bin/bash
#SBATCH -n 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
pseudo=pseudo_template
gene=gene_template
window_length=window_length_template
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
specific_gene_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups"

m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
treat_bdg=$m6mA_bam_dir/3D7_2rd_treat_pileup.bdg
input_bdg=$m6mA_bam_dir/3D7_2rd_control_lambda.bdg
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $specific_gene_global_dir/$pseudo/$gene
# rm -rf ${gene}_window50_6mA_FC
# for window in $(seq 1 50); do
#     window_fc=""

#     window_treat_intensity_cmd="grep \"\_${window}$\" $specific_gene_global_dir/$pseudo/$gene/${gene}_genes.window.bed | $bedtools intersect -a $treat_bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$8}'"
#     window_cont_intensity_cmd="grep \"\_${window}$\" $specific_gene_global_dir/$pseudo/$gene/${gene}_genes.window.bed | $bedtools intersect -a $input_bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$8}'"

#     window_fc+=" $($code_dir/s10_get_average_6mA_fc.r "$window_treat_intensity_cmd" "$window_cont_intensity_cmd")"
#     echo $window_fc >>${gene}_window50_6mA_FC
# done

# rm -rf ${gene}_${window_length}bp_6mA_FC
# for window in $(awk '{print $NF}' ${gene}_genes.peak.window.${window_length}bp.bed | sort -n | uniq); do
#     window_fc=""

#     window_treat_intensity_cmd="grep -e \"\\s${window}$\" $specific_gene_global_dir/$pseudo/$gene/${gene}_genes.peak.window.${window_length}bp.bed | $bedtools intersect -a $treat_bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$8}'"
#     window_cont_intensity_cmd="grep -e \"\\s${window}$\" $specific_gene_global_dir/$pseudo/$gene/${gene}_genes.peak.window.${window_length}bp.bed | $bedtools intersect -a $input_bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$8}'"

#     window_fc+=" $($code_dir/s10_get_average_6mA_fc.r "$window_treat_intensity_cmd" "$window_cont_intensity_cmd")"
#     echo $window_fc >>${gene}_${window_length}bp_6mA_FC
# done

rm -rf ${gene}_${window_length}bp_bam_6mA_FC
for window in $(awk '{print $NF}' ${gene}_genes.peak.window.${window_length}bp.bed | sort -n | uniq); do
    window_fc=""
    chrom=$(grep "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | awk 'NR==1{print $1}')
    grep "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | $bedtools intersect -b stdin -a $m6mA_bam_dir/../bam_depth/${chrom}.bam.ChIP.depth.bed >${window}.treat.bed
    grep "\s${window}$" ${gene}_genes.peak.window.${window_length}bp.bed | $bedtools intersect -b stdin -a $m6mA_bam_dir/../bam_depth/${chrom}.bam.Input.depth.bed >${window}.cont.bed

    window_fc+=" $(paste -d"\t" ${window}.treat.bed ${window}.cont.bed | awk '$NF!=0{print $4/$8}' | awk '{s+=$1}END{print s/NR}')"
    echo $window_fc >>${gene}_${window_length}bp_bam_6mA_FC
    rm -rf ${window}.treat.bed ${window}.cont.bed
done
