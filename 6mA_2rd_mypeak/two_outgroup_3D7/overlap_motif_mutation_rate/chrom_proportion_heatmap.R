#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(pheatmap)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/different_length_overlap_motif/snps_mean_variant_count_each_motif_peak", select = c(1, 3), stringsAsFactors = F)
summary_table <- as.matrix(table(data$V1, data$V3))
percent_summary_table <- t(apply(summary_table, 1, function(x) {
    x / sum(x)
}))
percent_summary_table1 <- t(apply(percent_summary_table, 2, function(x) {
    x / x[1]
}))
bk <- c(
    seq(0, 0.1, length.out = 50000),
    seq(0.100001, 0.2, length.out = 50), seq(0.17, 0.9, length.out = 2)
)
color <- c(
    colorRampPalette(c("white", "blue"))(50000),
    colorRampPalette(c("blue", "firebrick3"))(50), "firebrick3", "#7f0a0a"
)
pdf("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/different_length_overlap_motif/motif_proportion_in_chromosome.pdf", width = 14, height = 8)
# percent_summary_table[percent_summary_table == 0] <- NA
# pheatmap(percent_summary_table, color = colorRampPalette(c("#e5f5f9", "#99d8c9", "#2ca25f"))(50), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.1e",legend=F)

percent_summary_table[percent_summary_table == 0] <- NA
pheatmap(percent_summary_table[, -c(1, 2)],
    # color = colorRampPalette(c("#e5f5f9", "#99d8c9", "#2ca25f"))(50),
    color = color, breaks = bk,
    display_numbers = TRUE, number_format = "%.1e", cluster_cols = FALSE, legend = F, cluster_rows = F
)

dev.off()
