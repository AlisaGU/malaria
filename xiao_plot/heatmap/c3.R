#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/heatmap")

# 5. process -------------------------------------------------------------- TODO:
data <- read.table("Genelist202545rhoptry.txt", header = T, as.is = T)
data$order <- 1
data$gene_id <- factor(data$gene_id, levels = data$gene_id)


ggplot(data = data, aes(x = gene_id, y = order)) +
    geom_point(aes(color = log2FoldChange, size = -log10(pvalue))) +
    scale_color_continuous(low = "#864f30", high = "#c7a963") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Genelist202545rhoptry.pdf", width = 13, height = 5)
