#!/bin/bash
#SBATCH -n 20
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/tree/domain_tree"
hmmsearch="/home/gushanshan/anaconda3/bin/hmmsearch"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
mafft="/home/gushanshan/anaconda3/envs/vscode_r/bin/mafft"
raxmlHPC="/home/gushanshan/anaconda3/envs/ruby/bin/raxmlHPC"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
# $hmmsearch --cpu 4 -o /dev/null --domtblout domain_index.dom --noali --cut_ga /picb/evolgen/users/gushanshan/software/Hmmer/Pfam-A.hmm tree.OutGroup.remove2.add_PVC50508.1.fasta
# $code_dir/s1.R

# awk 'NR>1{print $1"\t"$4-1"\t"$5}' domain_index_in_protein_select_domain1.txt >domain_index_in_protein_select_domain1.bed
# awk 'NR>1{print $1"\t"$4-1"\t"$5}' domain_index_in_protein_select_domain2.txt >domain_index_in_protein_select_domain2.bed

# $seqkit subseq --bed domain_index_in_protein_select_domain1.bed tree.OutGroup.remove2.add_PVC50508.1.fasta >domain.seq.1.fasta
# $seqkit subseq --bed domain_index_in_protein_select_domain2.bed tree.OutGroup.remove2.add_PVC50508.1.fasta >domain.seq.2.fasta

# ## 只是改一下序列name，方便raxml建树
# sed -e 's/ /_/g' -e "s/\[//g" -e "s/]//g" -e "s/(//g" -e "s/)//g" -e "s/,//g" domain.seq.1.fasta >domain.seq.1.1.fasta
# sed -e 's/ /_/g' -e "s/\[//g" -e "s/]//g" -e "s/(//g" -e "s/)//g" -e "s/,//g" domain.seq.2.fasta >domain.seq.2.1.fasta
# 之后也手动修改了序列名
# ## 比对
# $mafft --auto domain.seq.1.1.fasta >domain.seq.1.1.mafft.fasta
# $mafft --auto domain.seq.2.1.fasta >domain.seq.2.1.mafft.fasta

## 建树
$raxmlHPC -f a -x 12345 -p 12345 -# 100 -o A.thaliana_BAA97534.1_30-362 -m PROTGAMMALGX -s domain.seq.1.1.mafft.fasta -n domain.seq.1.1.mafft.raxml.tree -T 20
# $raxmlHPC -f a -x 12345 -p 12345 -# 100 -o A.thaliana_BAA97534.1_30-362 -m PROTGAMMALGX -s domain.seq.2.1.mafft.fasta -n domain.seq.2.1.mafft.raxml.tree -T 20

## reorder fasta
$code_dir/consistent.R $working_dir/RAxML_bipartitions.domain.seq.1.1.mafft.raxml.tree $working_dir/domain.seq.1.1.mafft.fasta
$code_dir/consistent.R $working_dir/RAxML_bipartitions.domain.seq.2.1.mafft.raxml.tree $working_dir/domain.seq.2.1.mafft.fasta
