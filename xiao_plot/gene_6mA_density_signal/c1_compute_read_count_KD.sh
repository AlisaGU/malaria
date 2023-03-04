#!/bin/bash
#SBATCH -n 2
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

get_gene_6mA_density_signal() {
    chip_read_total_count=$(grep "paired in sequencing" $ | awk '{print $1}')
    input_read_total_count=$(grep "paired in sequencing" $input_summary | awk '{print $1}')
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $gene_depth_chip | $bedtools intersect -wa -a $all_genes_bed -wb -b stdin >$bam_dir/chip.inter
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' $gene_depth_input | $bedtools intersect -wa -a $all_genes_bed -wb -b stdin >$bam_dir/input.inter

}

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/gene_6mA_density_signal"
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal"
anno_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno"
bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD"
anno=$anno_dir/PlasmoDB-36_Pfalciparum3D7.gff
all_genes_bed=$working_dir/all_genes_bed
chip_summary=$bam_dir/chip_summary
input_summary=$bam_dir/input_summary

gene_depth_chip=$bam_dir/gene_depth_chip
gene_depth_input=$bam_dir/gene_depth_input

samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
$samtools flagstat $bam_dir/6mAKD-T3_ChIP.bam >$chip_summary
$samtools flagstat $bam_dir/6mAKD-T3_Input.bam >$input_summary
$samtools depth -a -b $all_genes_bed $bam_dir/6mAKD-T3_ChIP.bam >$gene_depth_chip
$samtools depth -a -b $all_genes_bed $bam_dir/6mAKD-T3_Input.bam >$gene_depth_input
get_gene_6mA_density_signal
# $code_dir/s1_compute.R
