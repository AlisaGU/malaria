#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
prefix="hdr"
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/MMS"
genes="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/all_genes_bed"
TSS2kb="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/genome_parts/TSS_2KB.bed"
TTS2kb="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/genome_parts/TTS_2KB.bed"
WT_chip=$working_dir/3D7-0-3_ChIP.bam
WT_input=$working_dir/3D7-0-3_Input.bam
KD_chip=$working_dir/3D7-mmsL-1_ChIP.bam
KD_input=$working_dir/3D7-mmsL-1_Input.bam

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
grep -f ${prefix}_list $genes | $bedtools makewindows -b stdin -n 300 -i srcwinnum >$prefix.gene.window
grep -f ${prefix}_list $TSS2kb | $bedtools makewindows -b stdin -n 100 -i srcwinnum >$prefix.tss.window
grep -f ${prefix}_list $TTS2kb | $bedtools makewindows -b stdin -n 100 -i srcwinnum >$prefix.tts.window

$samtools depth -a -b $prefix.gene.window $WT_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.gene.window -wb -b stdin >$prefix.gene.window.wt.chip
$samtools depth -a -b $prefix.gene.window $WT_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.gene.window -wb -b stdin >$prefix.gene.window.wt.input
$samtools depth -a -b $prefix.gene.window $KD_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.gene.window -wb -b stdin >$prefix.gene.window.kd.chip
$samtools depth -a -b $prefix.gene.window $KD_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.gene.window -wb -b stdin >$prefix.gene.window.kd.input

$samtools depth -a -b $prefix.tss.window $WT_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tss.window -wb -b stdin >$prefix.tss.window.wt.chip
$samtools depth -a -b $prefix.tss.window $WT_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tss.window -wb -b stdin >$prefix.tss.window.wt.input
$samtools depth -a -b $prefix.tss.window $KD_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tss.window -wb -b stdin >$prefix.tss.window.kd.chip
$samtools depth -a -b $prefix.tss.window $KD_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tss.window -wb -b stdin >$prefix.tss.window.kd.input

$samtools depth -a -b $prefix.tts.window $WT_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tts.window -wb -b stdin >$prefix.tts.window.wt.chip
$samtools depth -a -b $prefix.tts.window $WT_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tts.window -wb -b stdin >$prefix.tts.window.wt.input
$samtools depth -a -b $prefix.tts.window $KD_chip | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tts.window -wb -b stdin >$prefix.tts.window.kd.chip
$samtools depth -a -b $prefix.tts.window $KD_input | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | $bedtools intersect -wa -a $prefix.tts.window -wb -b stdin >$prefix.tts.window.kd.input
