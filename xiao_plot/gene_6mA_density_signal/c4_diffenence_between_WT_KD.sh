#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD"

TSS_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/genome_parts/TSS_2KB.bed"
tss_depth_chip=$bam_dir/tss_depth_chip
tss_depth_input=$bam_dir/tss_depth_input

samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# $samtools depth -a -b $TSS_bed $bam_dir/6mAKD-T3_ChIP.bam >$tss_depth_chip
# $samtools depth -a -b $TSS_bed $bam_dir/6mAKD-T3_Input.bam >$tss_depth_input
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $tss_depth_chip | $bedtools intersect -wa -a $TSS_bed -wb -b stdin >$bam_dir/tss.chip.inter
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $tss_depth_input | $bedtools intersect -wa -a $TSS_bed -wb -b stdin >$bam_dir/tss.input.inter
