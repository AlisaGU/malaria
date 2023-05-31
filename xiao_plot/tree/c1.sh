#!/bin/bash
#SBATCH -n 20
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1"
mafft="/home/gushanshan/anaconda3/envs/vscode_r/bin/mafft"
raxmlHPC="/home/gushanshan/anaconda3/envs/ruby/bin/raxmlHPC"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
# $mafft --auto tree.OutGroup.remove2.fasta >tree.OutGroup.remove2.mafft.fasta
# sed "s/ /_/g" tree.OutGroup.remove2.mafft.fasta >tree.OutGroup.remove2.mafft1.fasta
# $mafft --auto tree.OutGroup.fasta >tree.OutGroup.mafft.fasta
# sed "s/ /_/g" tree.OutGroup.mafft.fasta >tree.OutGroup.mafft1.fasta

# $mafft --auto tree.OutGroup.remove2.add_PVC50508.1.fasta >tree.OutGroup.remove2.add_PVC50508.1.mafft.fasta
sed -e 's/ /_/g' -e "s/\[//g" -e "s/]//g" -e "s/(//g" -e "s/)//g" -e "s/,//g" tree.OutGroup.remove2.add_PVC50508.1.mafft.fasta >tree.OutGroup.remove2.add_PVC50508.1.mafft1.fasta
# -f a
# 此参数用于选择 RAxML 运算的算法。可以设定的值非常之多。 a 表示执行快速 Bootstrap 分析并搜索最佳得分的 ML 树。
# -x 12345
# 指定一个 int 数作为随机种子，以启用快速 Bootstrap 算法。
# -p 12345
# 指定一个随机数作为 parsimony inferences 的种子。
# -# 100
# 指定 bootstrap 的次数。
# -m PROTGAMMALGX
# 指定核苷酸或氨基酸替代模型。PROTGAMMALGX 的解释： "PROT" 表示氨基酸替代模型； GAMMA 表示使用 GAMMA 模型； X 表示使用最大似然法估计碱基频率。

$raxmlHPC -f a -x 12345 -p 12345 -# 100 -o A.thaliana_BAA97534.1 -m PROTGAMMALGX -s tree.OutGroup.remove2.add_PVC50508.1.mafft1.fasta -n tree.OutGroup.remove2.add_PVC50508.1.mafft1.raxml.tree -T 20
