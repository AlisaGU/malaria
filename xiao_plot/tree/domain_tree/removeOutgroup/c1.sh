#!/bin/bash
#SBATCH -n 20
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree/remove_outgroup"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/tree/domain_tree"

mafft="/home/gushanshan/anaconda3/envs/vscode_r/bin/mafft"
raxmlHPC="/home/gushanshan/anaconda3/envs/ruby/bin/raxmlHPC"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
$raxmlHPC -f a -x 12345 -p 12345 -# 100 -o A.thaliana_BAA97534.1_30-362 -m PROTGAMMALGX -s domain.seq.1.1.mafft.fasta -n domain.seq.1.1.mafft.raxml.tree -T 20
$code_dir/consistent.R $working_dir/domain.seq.1.1.removeOutgroup.mafft.tree $working_dir/domain.seq.1.1.removeOutgroup.mafft.fasta
