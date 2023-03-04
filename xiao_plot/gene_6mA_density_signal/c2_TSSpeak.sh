#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal"
TSS_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/genome_parts/TSS_2KB.bed"
peak_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"
TSS_peakpro=$working_dir/TSS_peakpro.txt

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
$bedtools intersect -wao -a $TSS_bed -b $peak_bed | awk 'BEGIN{print "chrom\tTSS_start\tTSS_end\tgene_name\t.\tstrand\toverlap_base"};{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}' >$TSS_peakpro
