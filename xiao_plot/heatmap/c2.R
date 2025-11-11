#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(pheatmap)

# 2. functions ------------------------------------------------------------ TODO:


# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/heatmap")
data <- read.csv("heatmap_zscore_data.csv", header = T, as.is = T, row.names = 1)

# 5. process -------------------------------------------------------------- TODO:
data <- t(data)

bk <- c(
    min(data), seq(-1.45, -0.01, length.out = 49),
    0, seq(0.01, 1.45, length.out = 49), max(data)
)
color <- c(
    "#7A7ABC", colorRampPalette(c("#7A7ABC", "white"))(48), "white",
    colorRampPalette(c("white", "firebrick3"))(49), "firebrick3"
)

pdf("heatmap_zscore_data.pdf", height = 4, width = 8)
pheatmap(data,
    cluster_rows = F, cluster_cols = F, color = color,
    breaks = bk, border_color = "white",
    cellwidth = 10, cellheight = 15
)
dev.off()
