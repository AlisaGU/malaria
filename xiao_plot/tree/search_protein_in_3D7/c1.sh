#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/tree/search_protein_in_3D7"
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/search_protein_in_3D7"

proteins="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_AnnotatedProteins.fasta"

seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
cat $proteins | $seqkit locate -P -i -r -f proteinPattern.pp --bed >proteinPattern.pp.bed
cat $proteins | $seqkit locate -P -i -d -f proteinPattern.gxg --bed >proteinPattern.gxg.bed
$code_dir/s1.R
