#!/bin/bash
#SBATCH -n 10
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
RS_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos_RS/RS_bam"
binsize=5000

bamCompare="/home/gushanshan/anaconda3/envs/vscode_r/bin/bamCompare"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $RS_bam_dir
$samtools index -@ 10 3D7-4th-R_ChIP.bam
$samtools index -@ 10 3D7-4th-R_Input.bam
$samtools index -@ 10 3D7-S3_ChIP.bam
$samtools index -@ 10 3D7-S3_Input.bam
$samtools index -@ 10 6mAKD-4th-R_ChIP.bam
$samtools index -@ 10 6mAKD-4th-R_Input.bam
$samtools index -@ 10 6mAKD-S3_ChIP.bam
$samtools index -@ 10 6mAKD-S3_Input.bam

$bamCompare -b1 3D7-4th-R_ChIP.bam -b2 3D7-4th-R_Input.bam -o 3D7-4th-R_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
$bamCompare -b1 3D7-S3_ChIP.bam -b2 3D7-S3_Input.bam -o 3D7-S3_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
$bamCompare -b1 6mAKD-4th-R_ChIP.bam -b2 6mAKD-4th-R_Input.bam -o 6mAKD-4th-R_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
$bamCompare -b1 6mAKD-S3_ChIP.bam -b2 6mAKD-S3_Input.bam -o 6mAKD-S3_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 3D7-4th-R_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>2){print $1"\t"$2"\t"$3"\t"2}else{print $0}}' >3D7-4th-R_log2ratio_bin${binsize}_gt0_lt2.bed

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 3D7-S3_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>2){print $1"\t"$2"\t"$3"\t"2}else{print $0}}' >3D7-S3_log2ratio_bin${binsize}_gt0_lt2.bed

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 6mAKD-4th-R_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>2){print $1"\t"$2"\t"$3"\t"2}else{print $0}}' >6mAKD-4th-R_log2ratio_bin${binsize}_gt0_lt2.bed

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 6mAKD-S3_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>2){print $1"\t"$2"\t"$3"\t"2}else{print $0}}' >6mAKD-S3_log2ratio_bin${binsize}_gt0_lt2.bed
