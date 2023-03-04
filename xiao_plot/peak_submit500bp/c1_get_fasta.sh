#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/CT300-F-50"
genome_length="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
grep "^Pf3D7" CT300-F-50_peaks.xls | awk '{print $1"\t"$5-1"\t"$5}' | bedtools slop -i stdin -g $genome_length -b 250 >peak_submit500.bed
bedtools getfasta -fi /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta -bed peak_submit500.bed >peak_submit500.fa
