#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_TSS_FE() {
    $samtools depth -a -b $TSS_bed $bam_dir/3D7-T3_ChIP.bam >$tss_depth_chip
    $samtools depth -a -b $TSS_bed $bam_dir/3D7-T3_Input.bam >$tss_depth_input
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $tss_depth_chip | $bedtools intersect -wa -a $TSS_bed -wb -b stdin >$working_dir/tss.chip.inter
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $tss_depth_input | $bedtools intersect -wa -a $TSS_bed -wb -b stdin >$working_dir/tss.input.inter
}

compute_gene_SMRT() {
    $bedtools intersect -wa -a $TSS_bed -wb -b $SMRT | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"1}' >$working_dir/tss.smrt.bed
}

compute_TSS_SMRT() {
    $bedtools intersect -wa -a $all_genes_bed -wb -b $SMRT | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"1}' >$working_dir/gene.smrt.bed
}
# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal"
bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd"

TSS_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/genome_parts/TSS_2KB.bed"
all_genes_bed=$working_dir/all_genes_bed
tss_depth_chip=$working_dir/tss_depth_chip
tss_depth_input=$working_dir/tss_depth_input
SMRT="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/3D7-6mA-9Cell_FRACgt20_COVgt25.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
