#!/bin/bash
#SBATCH -n 20
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/domain_tree_onlyDNAmethylase"

seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
mafft="/home/gushanshan/anaconda3/envs/vscode_r/bin/mafft"
raxmlHPC="/home/gushanshan/anaconda3/envs/ruby/bin/raxmlHPC"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
# $mafft --auto domain.seq.remove1303100.fa >domain.seq.remove1303100.mafft.fa
$raxmlHPC -f a -x 12345 -p 12345 -# 100 -m PROTGAMMALGX -s domain.seq.remove1303100.mafft.fa -n domain.seq.remove1303100.mafft.raxml.tree -T 20
