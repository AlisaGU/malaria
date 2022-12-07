# 1. modules and external scripts ----------------------------------------- TODO:
import os
import numpy as np
# 2. functions ------------------------------------------------------------ TODO:

# 3. input ---------------------------------------------------------------- TODO:

# 4. process -------------------------------------------------------------- TODO:
os.chdir("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot")

data = np.load('CT300-2-500-all_3D7-T3_ChIP.npz')
data.files ## 查看npz文件里都有什么文件
np.savetxt("CT300-2-500-all_3D7-T3_ChIP_all_downstream_gene_2kb.txt", data['CT300-2-500-all_downstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_3D7-T3_ChIP_all_upstream_gene_2kb.txt", data['CT300-2-500-all_upstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_3D7-T3_ChIP_all_full_genes.txt", data['CT300-2-500-all_full_genes'], fmt = '%f', delimiter = '\t')


data = np.load('CT300-2-500-all_3D7-T3_Input.npz')
data.files ## 查看npz文件里都有什么文件
np.savetxt("CT300-2-500-all_3D7-T3_Input_all_downstream_gene_2kb.txt", data['CT300-2-500-all_downstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_3D7-T3_Input_all_upstream_gene_2kb.txt", data['CT300-2-500-all_upstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_3D7-T3_Input_all_full_genes.txt", data['CT300-2-500-all_full_genes'], fmt = '%f', delimiter = '\t')

data = np.load('CT300-2-500-all_6mAKD-T3_ChIP.npz')
data.files ## 查看npz文件里都有什么文件
np.savetxt("CT300-2-500-all_6mAKD-T3_ChIP_all_downstream_gene_2kb.txt", data['CT300-2-500-all_downstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_6mAKD-T3_ChIP_all_upstream_gene_2kb.txt", data['CT300-2-500-all_upstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_6mAKD-T3_ChIP_all_full_genes.txt", data['CT300-2-500-all_full_genes'], fmt = '%f', delimiter = '\t')


data = np.load('CT300-2-500-all_6mAKD-T3_Input.npz')
data.files ## 查看npz文件里都有什么文件
np.savetxt("CT300-2-500-all_6mAKD-T3_Input_all_downstream_gene_2kb.txt", data['CT300-2-500-all_downstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_6mAKD-T3_Input_all_upstream_gene_2kb.txt", data['CT300-2-500-all_upstream_gene_2kb'], fmt = '%f', delimiter = '\t')
np.savetxt("CT300-2-500-all_6mAKD-T3_Input_all_full_genes.txt", data['CT300-2-500-all_full_genes'], fmt = '%f', delimiter = '\t')