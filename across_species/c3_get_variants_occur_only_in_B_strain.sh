#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/across_species"
chroms="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
cd $code_dir

# species3
for chrom in $chroms; do
    #     sed -e "s/chrom_template/${chrom}/g" $code_dir/s3_get_mutation_load_for_one_chrom_species3.sh >s3_get_mutation_load_for_one_chrom_species3.sh_${chrom}
    #     sbatch s3_get_mutation_load_for_one_chrom_species3.sh_${chrom}
    sed -e "s/chrom_template/${chrom}/g" $code_dir/s3_get_mutation_load_for_one_chrom_species3_include_indel.sh >s3_get_mutation_load_for_one_chrom_species3_include_indel.sh_${chrom}
    sbatch s3_get_mutation_load_for_one_chrom_species3_include_indel.sh_${chrom}
done

# species2
# for species2 in $(echo "P_falciparum.P_reichenowi P_falciparum.P_praefalciparum"); do
# for species2 in $(echo "P_falciparum.P_reichenowi"); do
#     for chrom in $chroms; do
#         sed -e "s/chrom_template/${chrom}/g" -e "s/species2_template/$species2/g" $code_dir/s3_get_mutation_load_for_one_chrom_species2.sh >s3_get_mutation_load_for_one_chrom_species2.sh_${chrom}
#         sbatch s3_get_mutation_load_for_one_chrom_species2.sh_${chrom}
#     done
# done
