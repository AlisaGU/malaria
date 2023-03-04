#!/bin/bash
#SBATCH -n 4
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/circos"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
time $code_dir/s2.R
