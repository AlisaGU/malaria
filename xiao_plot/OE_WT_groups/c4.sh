#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/xiao_plot/OE_WT_groups"
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/xiao_plot/OE_WT_groups"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
$code_dir/s4_get_line.R "class1"
$code_dir/s4_get_line.R "class5"

$code_dir/s4_get_line.R "RNA_class1"
$code_dir/s4_get_line.R "RNA_class5"
