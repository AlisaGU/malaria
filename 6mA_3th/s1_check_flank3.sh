#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern/flank3_test"
m6mA_3rd_bed="../../3D7-6mA-9Cell.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $working_dir
cd $working_dir

flank3_count=$(wc -l ../motif_3_6_9.txt | awk '{print $1}')
for rowth in $(seq 1 $flank3_count); do
    awk -v rowth=$rowth 'NR==rowth{print $1"\t"$21-1"\t"$21}' ../motif_3_6_9.txt | $bedtools intersect -a stdin -b $m6mA_3rd_bed
done

awk '{print $1"\t"$21-1"\t"$21"\t"$2"\t"$3"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}' ../motif_3_6_9.txt | $bedtools intersect -a stdin -b $m6mA_3rd_bed -wb

awk '{print $1"\t"$22-1"\t"$22}' ../motif_3_6_9.txt | $bedtools intersect -a stdin -b $m6mA_3rd_bed
awk '{print $1"\t"$21-1"\t"$21}' ../motif_3_6_9.txt | $bedtools intersect -a stdin -b $m6mA_3rd_bed
