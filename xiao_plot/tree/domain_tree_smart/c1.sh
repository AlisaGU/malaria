#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree_smart"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
mafft="/home/gushanshan/anaconda3/envs/vscode_r/bin/mafft"
raxmlHPC="/home/gushanshan/anaconda3/envs/ruby/bin/raxmlHPC"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
# awk 'NR>1{print $1"\t"$4-1"\t"$5}' domain_index_in_protein.txt >domain_index_in_protein.bed
# $seqkit subseq --bed domain_index_in_protein.bed ../domain_tree/tree.OutGroup.remove2.add_PVC50508.1.fasta >domain.seq.fasta
# cp domain.seq.fasta domain.seq.removeOutgroup.fasta
# ## 之后手动删掉
# $mafft --auto domain.seq.removeOutgroup.fasta >domain.seq.removeOutgroup.mafft.fasta

# ## sort序列
# grep "^>" ../domain_tree/remove_outgroup/domain.seq.1.1.removeOutgroup.mafft.mega.fas
