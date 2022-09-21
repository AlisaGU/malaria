#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
anno_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri"

bed_dir=$global_output_dir/genome_parts
gene_bed=$bed_dir/genes.bed
TSS_2KB=$bed_dir/TSS_2KB.bed ##转录起始调控
TTS_2KB=$bed_dir/TTS_2KB.bed ##转录终止调控
exon_bed=$bed_dir/exons.bed
intron_bed=$bed_dir/introns.bed
intergenic_bed=$bed_dir/intergenic.bed
genome_len=$anno_dir/../genome/chrom.length

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
mkdir -p $bed_dir
cd $bed_dir
##获得转录起始位点TSS前2kb区间
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"".""\t"$5}' $gene_bed | $bedtools flank -i stdin -g $genome_len -s -l 2000 -r 0 >$TSS_2KB
##获得转录终止位点TTS后2kb区间
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"".""\t"$5}' $gene_bed | $bedtools flank -i stdin -g $genome_len -s -l 0 -r 2000 >$TTS_2KB
##获得基因间区
awk '{print $1"\t"0"\t"$2}' $genome_len >${genome_len}.medium
$bedtools slop -i $gene_bed -g $genome_len -b 2000 | $bedtools merge -i stdin | $bedtools subtract -a ${genome_len}.medium -b stdin >$intergenic_bed
rm -rf ${genome_len}.medium
