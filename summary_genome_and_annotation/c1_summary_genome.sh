#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome
seqkit fx2tab -l  -n -i -H PlasmoDB-36_Pfalciparum3D7_Genome.fasta