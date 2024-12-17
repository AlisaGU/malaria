#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/drugResistance/densityPlot"
groups="s16_S1A_6AR_vs_6AS s16_S1B_11CR_vs_11CS s17_S2A_6AR_vs_6AS_ART s17_S2B_11CR_vs_11CS_ART"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
for group in $(echo $groups | tr " " "\n"); do
    $code_dir/s1.Rstage.R $group
done
