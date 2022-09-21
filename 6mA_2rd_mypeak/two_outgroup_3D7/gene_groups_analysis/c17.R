#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(pheatmap)

# 2. functions ------------------------------------------------------------ TODO:
read_data <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(mutation_density_global_dir, "/", gene_set, "/", "snps_mean_variant_count_window"), header = F, as.is = T, sep = ";")
        result <- lapply(1:nrow(data), function(ith) {
            x <- unlist(strsplit(data[ith, ], " "))
            space_index <- which(x == "space")

            if (space_index != 2) {
                # RIF 家族的基因PF3D7_0222600与信号区只有2bp的交集，没办法滑窗，只能去掉
                peak_data <- x[2:(space_index - 1)]
                control_data <- x[(space_index + 1):length(x)]

                peak_data_converted <- transvert_to_10(value = as.numeric(peak_data))
                control_data_converted <- transvert_to_10(value = as.numeric(control_data))
                result <- c(x[1], peak_data_converted, control_data_converted, gene_set)
                return(result)
            }
        })
        result <- do.call(rbind, result)
        return(result)
    })
    result <- do.call(rbind, result)
}
transvert_to_10 <- function(value = NULL) {
    fold <- length(value) / 20
    if (fold == 1) {
        return(value)
    } else {
        result <- as.vector(sapply(1:10, function(y) {
            data1 <- value[(1 + 2 * fold * (y - 1)):(2 * fold * y)]
            return(c(
                sum(data1[seq(1, length(data1), by = 2)]), sum(data1[seq(2, length(data1), by = 2)])
            ))
        }))
        return(result)
    }
}
# 3. input ---------------------------------------------------------------- TODO:
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets"

gene_sets <- c("HDR", "RNA_translation", "VAR", "RIF", "STEVOR")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read_data()
ann_colors <- list(
    Type = c("6mA" = "#ffa500", Control = "#BEBEBE"),
    gene_Class = c(HDR = "#984ea3", RNA_translation = "#4DBBD5", RIF = "#E64B35", STEVOR = "#FDBF6F", VAR = "#4DAF4A")
)
## 全部数据
data_for_plot <- data[, 2:41]
data_for_plot <- apply(data_for_plot, 2, function(x) {
    return(as.integer(x))
})
data_for_plot <- data_for_plot[, seq(1, ncol(data_for_plot), by = 2)] / data_for_plot[, seq(2, ncol(data_for_plot), by = 2)]
colnames(data_for_plot) <- c(paste("peak_window", 1:10, sep = ""), paste("control_window", 1:10, sep = ""))
rownames(data_for_plot) <- data[, 1]

annotation_col <- data.frame("Type" = rep(c("6mA", "Control"), each = 10))
rownames(annotation_col) <- colnames(data_for_plot)
annotation_row <- data.frame("gene_Class" = data[, 42])
rownames(annotation_row) <- rownames(data_for_plot)
rownames(annotation_col) <- colnames(data_for_plot)

pheatmap(data_for_plot,
    cluster_rows = F, cluster_cols = F,
    annotation_col = annotation_col, annotation_row = annotation_row,
    show_rownames = F, show_colnames = F, breaks = seq(0, 0.1, length.out = 100),
    color = colorRampPalette(c("#cde3f6", "#083572"))(100),
    na_col = "white"
)

# pheatmap(data_for_plot,
#     cluster_rows = F, cluster_cols = F,
#     annotation_col = annotation_col, annotation_row = annotation_row,
#     show_rownames = F, show_colnames = F, breaks = seq(0, 1.5, length.out = 100),
#     color = colorRampPalette(c("#cde3f6", "#083572"))(100),
#     na_col = "white", border_color = "grey",annotation_colors = ann_colors,filename="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets/heatmap.pdf",width=8,height=10
# )

## 有peak的数据
data <- data[!is.na(data[, 2]), ]
data_for_plot <- data[, 2:41]
data_for_plot <- apply(data_for_plot, 2, function(x) {
    return(as.integer(x))
})
data_for_plot <- data_for_plot[, seq(1, ncol(data_for_plot), by = 2)] / data_for_plot[, seq(2, ncol(data_for_plot), by = 2)]
colnames(data_for_plot) <- c(paste("peak_window", 1:10, sep = ""), paste("control_window", 1:10, sep = ""))
rownames(data_for_plot) <- data[, 1]

annotation_col <- data.frame("Type" = rep(c("6mA", "Control"), each = 10))
rownames(annotation_col) <- colnames(data_for_plot)
annotation_row <- data.frame("gene_Class" = data[, 42])
rownames(annotation_row) <- rownames(data_for_plot)
rownames(annotation_col) <- colnames(data_for_plot)

pheatmap(data_for_plot,
    cluster_rows = F, cluster_cols = F,
    annotation_col = annotation_col, annotation_row = annotation_row,
    show_rownames = F, show_colnames = F, breaks = seq(0, 1.5, length.out = 100),
    color = colorRampPalette(c("white", "#083572"))(100),
    na_col = "white", border_color = "grey", annotation_colors = ann_colors, filename = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets/heatmap_only_with_peak.pdf", width = 8, height = 10
)

pheatmap(data_for_plot,
    cluster_rows = T, cluster_cols = F,
    annotation_col = annotation_col, annotation_row = annotation_row,
    show_rownames = F, show_colnames = F, breaks = seq(0, 1.5, length.out = 100),
    color = colorRampPalette(c("white", "#083572"))(100),
    na_col = "white", border_color = "grey", annotation_colors = ann_colors, filename = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets/heatmap_only_with_peak_cluster.pdf", width = 8, height = 10
)
