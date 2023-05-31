#!/bin/bash
#SBATCH -n 10
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/offtarget/BAM"
binsize=500

bamCompare="/home/gushanshan/anaconda3/envs/vscode_r/bin/bamCompare"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $bam_dir
# $samtools index -@ 10 6mAKD-T4-Flag_ChIP.bam
# $samtools index -@ 10 6mAKD-T4_Input.bam
# $samtools index -@ 10 3D7-T4_Input.bam
# $samtools index -@ 10 3D7-T4-Flag_ChIP.bam

# $bedtools genomecov -bga -ibam 6mAKD-T4-Flag_ChIP.bam | awk '$4<10' >6mAKD-T4-Flag_ChIP.lowCoverage.bed
# $bedtools genomecov -bga -ibam 6mAKD-T4_Input.bam | awk '$4<10' >6mAKD-T4_Input.lowCoverage.bed
# $bedtools genomecov -bga -ibam 3D7-T4_Input.bam | awk '$4<10' >3D7-T4_Input.lowCoverage.bed
# $bedtools genomecov -bga -ibam 3D7-T4-Flag_ChIP.bam | awk '$4<10' >3D7-T4-Flag_ChIP.lowCoverage.bed

# cat 6mAKD-T4-Flag_ChIP.lowCoverage.bed 6mAKD-T4_Input.lowCoverage.bed | sort -k1,1 -k2n | $bedtools merge -i stdin >6mAKD-T4.lowCoverage.bed
# cat 3D7-T4_Input.lowCoverage.bed 3D7-T4-Flag_ChIP.lowCoverage.bed | sort -k1,1 -k2n | $bedtools merge -i stdin >3D7-T4.lowCoverage.bed

$bamCompare -b1 6mAKD-T4-Flag_ChIP.bam -b2 6mAKD-T4_Input.bam --blackListFileName 6mAKD-T4.lowCoverage.bed -o 6mAKD-T4_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
$bamCompare -b1 3D7-T4-Flag_ChIP.bam -b2 3D7-T4_Input.bam --blackListFileName 3D7-T4.lowCoverage.bed -o 3D7-T4_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 6mAKD-T4_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' >6mAKD-T4_log2ratio_bin${binsize}_gt0.bdg

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 6mAKD-T4_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>4.5){print $1"\t"$2"\t"$3"\t"4.5}else{print $0}}' >6mAKD-T4_log2ratio_bin${binsize}_gt0_lt4.5.bdg

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 3D7-T4_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' >3D7-T4_log2ratio_bin${binsize}_gt0.bdg

sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" 3D7-T4_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>4.5){print $1"\t"$2"\t"$3"\t"4.5}else{print $0}}' >3D7-T4_log2ratio_bin${binsize}_gt0_lt4.5.bdg

## wt-kd
mkdir -p KD_WT_bin${binsize}
cd KD_WT_bin${binsize}/
$bedtools intersect -a ../3D7-T4_log2ratio_bin${binsize}.bdg -b ../6mAKD-T4_log2ratio_bin${binsize}.bdg | sort -k1,1 -k2n | awk '{if($4<0){print $1"\t"$2"\t"$3"\t"0}else{print $0}}' >WT_intersect.bdg
$bedtools intersect -b ../3D7-T4_log2ratio_bin${binsize}.bdg -a ../6mAKD-T4_log2ratio_bin${binsize}.bdg | sort -k1,1 -k2n | awk '{if($4<0){print $1"\t"$2"\t"$3"\t"0}else{print $0}}' >KD_intersect.bdg
paste WT_intersect.bdg KD_intersect.bdg | awk '{print $1"\t"$2"\t"$3"\t"$8-$4}' | sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" | awk '$4>1 && $3-$2>200{print}' >KD_WT_1.bed
