#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(data.table)
library(purrr)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function(gene_groups_pattern = NULL) {
    file_names <- paste(gene_groups_pattern, c("allvar", "snps", "indels"), c("mean_variant_count"), sep = "_")
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        no_zero_chrom_index <- apply(a, 1, function(x) {
            all(x[seq(2, length(x), by = 2)] != 0)
        })
        a <- a[no_zero_chrom_index, ]

        data <- lapply(seq(2, ncol(a), by = 2), function(y) {
            a[, y - 1] / a[, y]
        })
        data <- do.call(cbind, data)
        colnames(data) <- unlist(strsplit(gene_groups_pattern, "_"))
        rownames(data) <- rownames(a)

        data_for_plot <- data.frame(
            value = as.vector(data),
            chrom = rep(rownames(data), times = ncol(data)),
            source = rep(colnames(data), each = nrow(data))
        )
        a <- gsub("chrom", "", data_for_plot$chrom)
        a[a != "1"] <- ""
        data_for_plot$label <- a
        data_for_plot$variant <- gsub("_mean_variant_count", "", gsub(paste0(gene_groups_pattern, "_"), "", x))
        data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
        data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
            "All variants", "SNPs", "Indels",
            "A to G", "A to T", "A to C", "A to sDel",
            "G to A", "G to T", "G to C", "G to sDel",
            "C to T", "C to A", "C to G", "C to sDel",
            "T to C", "T to A", "T to G", "T to sDel"
        ))
        return(data_for_plot)
    })
    result <- do.call(rbind, data)
    return(result)
}

myplot <- function(data_for_plot = NULL, gene_groups_pattern = NULL) {
    # color_value_N2N <- c("control_motif" = "black", "peak_motif" = "#b2182b")
    gene_groups <- unlist(strsplit(gene_groups_pattern, "_"))
    gene_groups_combi <- combn(gene_groups, 2)
    my_comparisons <- lapply(1:ncol(gene_groups_combi), function(x) {
        gene_groups_combi[, x]
    })
    p <- ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(size = 8) +
        # scale_color_manual(values = color_value_N2N, labels = c("Control motif", "Peak motif")) +
        stat_compare_means(
            comparisons = my_comparisons,
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 15
        ) +
        facet_wrap(~variant_f, nrow = 1, ncol = length(gene_groups), scale = "free_x") +
        theme_bw() +
        labs(y = "Mutation density") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_groups"
setwd(path)

Args <- commandArgs()
gene_groups_pattern <- Args[6]
# 4. variable setting of test module--------------------------------------- TODO:
# gene_groups_pattern <- "VAR_HDR_invasion"


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read_data(gene_groups_pattern = gene_groups_pattern)
myplot(data_for_plot = data_for_plot, gene_groups_pattern = gene_groups_pattern)
ggsave(paste0(gene_groups_pattern, "_global_variants.pdf"), width = 24, height = 8)
