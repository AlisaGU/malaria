#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_exclude_pseudo")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
gene_sets <- c("DR", "HDR", "RNA_translation", "STEVOR", "VAR", "RIF")
data <- lapply(gene_sets, function(gene_set) {
    a <- read.table(paste0(gene_set, "/", gene_set, "_motif_info"), header = F, as.is = T)
    result <- cbind(a, gene_set)
    return(result)
})
data <- do.call(rbind, data)
colnames(data)[1:4] <- c("gene", "gene_len", "motif_count", "merged_motif_base_count")
data$average_motif_count <- data$motif_count / data$gene_len
data$average_merged_motif_base_count <- data$merged_motif_base_count / data$gene_len
ggplot(data) +
    geom_boxplot(aes(x = gene_set, y = average_motif_count, color = gene_set))

ggplot(data) +
    geom_boxplot(aes(x = gene_set, y = merged_motif_base_count, color = gene_set))
