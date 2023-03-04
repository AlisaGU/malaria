#!/bin/bash
#SBATCH -n 10
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
enzyme_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/histone_modification/protein"
WT_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd"
KD_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD"
GAAGAA_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/GAAGAA_genome_position"
mutation_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/mutation_density_distribution/wholeGenome/ESEA_WSEA_OCE_SAM_SAS/variant_3D7"

binsize=5000

bamCompare="/home/gushanshan/anaconda3/envs/vscode_r/bin/bamCompare"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $enzyme_dir
# $samtools index CT300-Flag_ChIP.bam
# $samtools index CT300-2_Input.bam
$bamCompare -b1 CT300-Flag_ChIP.bam -b2 CT300-2_Input.bam -o enzyme_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
# awk '{print $1"\t"$2"\t"$3}' enzyme_log2ratio_bin${binsize}.bdg | md5sum
sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" enzyme_log2ratio_bin${binsize}.bdg | awk '$4>0{print $0}' | awk '{if($4>1){print $1"\t"$2"\t"$3"\t"1}else{print $0}}' >enzyme_log2ratio_bin${binsize}_max1_min0.bed

cd $WT_dir
# $samtools index 3D7-T3_ChIP.bam
# $samtools index 3D7-T3_Input.bam
# $bamCompare -b1 3D7-T3_ChIP.bam -b2 3D7-T3_Input.bam -o WT_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10

# awk '{print $1"\t"$2"\t"$3}' WT_log2ratio_bin${binsize}.bdg | md5sum

cd $KD_dir
# $samtools index 6mAKD-T3_ChIP.bam
# $samtools index 6mAKD-T3_Input.bam
# $bamCompare -b1 6mAKD-T3_ChIP.bam -b2 6mAKD-T3_Input.bam -o KD_log2ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10
$bamCompare -b1 6mAKD-T3_ChIP.bam -b2 6mAKD-T3_Input.bam -o KD_ratio_bin${binsize}.bdg --outFileFormat=bedgraph --binSize=$binsize --numberOfProcessors 10 --Operation=="ratio"

# awk '{print $1"\t"$2"\t"$3}' KD_log2ratio_bin${binsize}.bdg | md5sum

# cd $GAAGAA_dir
# awk '{print $1"\t"$2"\t"$3}' $KD_dir/KD_log2ratio_bin${binsize}.bdg | $bedtools intersect -a stdin -b GAAGAA.bed -c >GAAGAA_bin${binsize}_count.bdg
sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" GAAGAA_bin${binsize}_count.bdg | awk '{if($4>100){print $1"\t"$2"\t"$3"\t"100}else{print $0}}' >GAAGAA_bin${binsize}_count_max100.bed

cd $mutation_dir
cat *bed | sort -k1,1 -k2n >all_chroms.bed
awk '{print $1"\t"$2"\t"$3}' $KD_dir/KD_log2ratio_bin${binsize}.bdg | $bedtools intersect -a stdin -b all_chroms.bed -wa -wb >all_chrom_bin${binsize}_tmp.bed
/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/circos/s1.R $mutation_dir/all_chrom_bin${binsize}_tmp.bed
rm -rf all_chrom_bin${binsize}_tmp.bed all_chroms.bed
sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" all_chrom_bin${binsize}_indel.bed | awk '{if($4>600){print $1"\t"$2"\t"$3"\t"600}else{print $0}}' >all_chrom_bin${binsize}_indel_pf_max600.bed
sed -e "s/Pf3D7_01_v3/pf1/g" -e "s/Pf3D7_02_v3/pf2/g" -e "s/Pf3D7_03_v3/pf3/g" -e "s/Pf3D7_04_v3/pf4/g" -e "s/Pf3D7_05_v3/pf5/g" -e "s/Pf3D7_06_v3/pf6/g" -e "s/Pf3D7_07_v3/pf7/g" -e "s/Pf3D7_08_v3/pf8/g" -e "s/Pf3D7_09_v3/pf9/g" -e "s/Pf3D7_10_v3/pf10/g" -e "s/Pf3D7_11_v3/pf11/g" -e "s/Pf3D7_12_v3/pf12/g" -e "s/Pf3D7_13_v3/pf13/g" -e "s/Pf3D7_14_v3/pf14/g" -e "s/Pf_M76611/pf15/g" -e "s/Pf3D7_API_v3/pf16/g" all_chrom_bin${binsize}_snp.bed | awk '{if($4>400){print $1"\t"$2"\t"$3"\t"400}else{print $0}}' >all_chrom_bin${binsize}_snp_pf_max400.bed
