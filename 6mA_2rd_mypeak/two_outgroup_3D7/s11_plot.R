#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_motif_data <- function(filename = NULL, variant_type = NULL) {
    # files <- "snps_mean_variant_count"
    files <- filename
    data <- read.table(files, as.is = T)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    colnames(data1) <- c("peak_motif", "control_motif", "whole_peak_region", "whole_control_region")
    rownames(data1) <- paste("chrom", 1:14, sep = "")

    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    a <- gsub("chrom", "", data_for_plot$chrom)
    a[a != "1"] <- ""
    data_for_plot$label <- a
    # data_for_plot$variant <- "SNPs"
    data_for_plot$variant <- variant_type

    data_for_plot$source <- factor(data_for_plot$source, levels = c("peak_motif", "whole_peak_region", "control_motif", "whole_control_region"))
    return(data_for_plot)
}

plot_motif <- function(data_for_plot = NULL, title = NULL) {
    color_value_N2N <- c("control_motif" = "black", "peak_motif" = "#b2182b", "whole_peak_region" = "#d75665", "whole_control_region" = "grey")
    p <- ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 8) +
        scale_color_manual(values = color_value_N2N, labels = c("Peak motif", "Whole peak region", "Control motif", "Whole control region")) +
        scale_x_discrete(labels = c("Peak motif", "Whole peak region", "Control motif", "Whole control region")) +
        stat_compare_means(
            comparisons = list(c("peak_motif", "whole_peak_region"), c("peak_motif", "control_motif"), c("control_motif", "whole_control_region"), c("peak_motif", "whole_control_region")),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 14
        ) +
        labs(y = "Mutation density", title = title) +
        theme_bw() +
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
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}


# 3. input ---------------------------------------------------------------- TODO:

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_whole_region_compare")
data_for_plot <- read_motif_data(filename = "snps_mean_variant_count", variant_type = "SNPs")
snps <- plot_motif(data_for_plot = data_for_plot, title = "SNPs")

data_for_plot <- read_motif_data(filename = "single_del_mean_variant_count", variant_type = "Single_deletion")
single_del <- plot_motif(data_for_plot = data_for_plot, title = "Single deletion")
snps / single_del
ggsave("motif_whole_region_compare.pdf", width = 16, height = 14)
