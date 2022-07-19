#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
gene_group_overlap_peak_related_count <- as.integer(Args[6]) # 取的球里有几个白的
peak_related_gene_count <- as.integer(Args[7]) # 白球
whole_genome_gene_count_except_peak_related <- as.integer(Args[8]) # 黑球
gene_group_count <- as.integer(Args[9]) # 取了几个球

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
p <- 1 - phyper(gene_group_overlap_peak_related_count - 1, peak_related_gene_count, whole_genome_gene_count_except_peak_related, gene_group_count)
cat(p)
