#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7/site_by_site_whole_genome"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
flank_sites=$(seq 1 10)
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -e

# 8. 计算A、T、C、G四种碱基在信号区的基底速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_each_chrom_wholeGenome_base_mutation_density.sh >$code_dir/s7_compute_each_chrom_wholeGenome_base_mutation_density.sh_$popu_symbol
    sbatch $code_dir/s7_compute_each_chrom_wholeGenome_base_mutation_density.sh_$popu_symbol
done
