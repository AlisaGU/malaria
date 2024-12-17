#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_ZbynekBozdech_Rstage"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/gene_6mA_density_signal/drug_resistance"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir/

for prefix in $(ls *list | awk -F "." '{print $1}' | sort -u); do
    $code_dir/s2.R $working_dir $prefix
done
