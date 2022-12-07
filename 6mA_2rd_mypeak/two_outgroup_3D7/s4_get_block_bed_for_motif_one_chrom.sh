#!/bin/bash
#SBATCH -n 2
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
bed=bed_template
region_bed_dir=region_bed_dir_template
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
$code_dir/s4_get_block_bed_for_motif.R $bed $region_bed_dir
