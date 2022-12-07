#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
two_outgroup_consistent_peak_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/allchroms.peak.bed"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
whole_region_length=$(awk '{print $3-$2}' $two_outgroup_consistent_peak_bed | awk '{s+=$1} END {print s}')
A_count=$(seqkit subseq --bed $two_outgroup_consistent_peak_bed $genome | seqkit fx2tab -n -i -H -C A | awk 'NR>1{s+=$2} END {print s}')
T_count=$(seqkit subseq --bed $two_outgroup_consistent_peak_bed $genome | seqkit fx2tab -n -i -H -C T | awk '{s+=$2} END {print s}')
C_count=$(seqkit subseq --bed $two_outgroup_consistent_peak_bed $genome | seqkit fx2tab -n -i -H -C C | awk '{s+=$2} END {print s}')
G_count=$(seqkit subseq --bed $two_outgroup_consistent_peak_bed $genome | seqkit fx2tab -n -i -H -C G | awk '{s+=$2} END {print s}')

$code_dir/s7_create_random_base_combination.r $whole_region_length $A_count $T_count $C_count $G_count
