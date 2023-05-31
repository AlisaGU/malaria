#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree/")
doms <- fread(cmd = paste0("zcat /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree/domain_index.dom.gz|grep -v \"^#\""), header = F, stringsAsFactors = F)
colnames(doms) <- c(
    "targetName", "targetAccession", "targetLen", "queryName", "queryAccession", "queryLen",
    "fullSequence.Evalue", "fullSequence.score", "fullSequence.bias",
    "thisDomain.rank", "thisDomain.totalNumber", "thisDomain.c-Evalue", "thisDomain.i-Evalue", "thisDomain.score", "thisDomain.bias",
    "hmm.from", "hmm.to", "ali.from", "ali.to", "env.from", "env.to", "expectedAccuracy", "targetDescription"
)

result <- doms[, c("targetName", "targetLen", "queryName", "env.from", "env.to")]
result <- result[order(result$targetName), ]
colnames(result) <- c("protein", "proteinLength", "domainName", "start_of_domain_in_protein", "end_of_domain_in_protein")
write.table(result, "domain_index_in_protein.txt", row.names = F, col.names = T, sep = "\t")
