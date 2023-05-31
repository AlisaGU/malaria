#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/xiao_plot/gene_6mA_density_signal/6_list"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
for prefix in $(echo "merged MMR NER BER MMEJ HDR" | tr " " "\n"); do
    $code_dir/c3_get_line_mmsl.R $prefix
    $code_dir/c3_get_line_t3.R $prefix
done
