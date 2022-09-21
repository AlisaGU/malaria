#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(data.table)
# 2. functions ------------------------------------------------------------
get_data_for_plot <- function(dataname = NULL) {
    data <- read.table(dataname, stringsAsFactors = F, sep = " ")
    colnames(data) <- c(
        "nopeak", paste("w", 1:10, sep = "")
    )
    rownames(data) <- paste("chrom", 1:14, sep = "")

    data_for_plot <- data.frame(
        value = unlist(data),
        chrom = rep(rownames(data), times = ncol(data)),
        source = rep(colnames(data), each = nrow(data))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c("nopeak", paste("w", 1:10, sep = "")))

    a <- gsub("chrom", "", data_for_plot$chrom)
    a[a != "1"] <- ""
    data_for_plot$label <- a
    return(data_for_plot)
}


plot_sliding_window_mutation_density <- function(data_for_plot = NULL) {
    cols <- brewer.pal(3, "YlOrRd")
    pal <- colorRampPalette(cols)
    mycolors <- rev(pal(10))
    color_value <- c("#000000", "#000000", mycolors)
    names(color_value) <- c("nopeak", "nopeak500", paste("w", 1:10, sep = ""))

    p <- ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
        geom_boxplot(
            position = position_dodge(), outlier.shape = NA, width = 0.5, size = 1
        ) +
        geom_point(
            position = position_jitterdodge(), size = 6
        ) +
        facet_wrap(~type, nrow = 1, ncol = 2, scale = "free") +
        scale_color_manual(values = color_value) +
        stat_compare_means(
            method = "t.test", tip.length = 0, bracket.size = 0,
            comparisons = list(c("nopeak", "w1")),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p = {p.format}", paired = T, size = 15
        ) +
        scale_x_discrete(labels = c("Control", "Submit~5%", "5~10%", "10~15%", "15~20%", "20~25%", "25~30%", "30~35%", "35~40%", "40~45%", "45~50%")) +
        labs(y = "Mutation density") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black",
                angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

# 3. input ----------------------------------------------------------------
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/sliding_window_mutation_density")
# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
data_snp <- get_data_for_plot(dataname = "snps_mean_variant_count")
data_indel <- get_data_for_plot(dataname = "indels_mean_variant_count")
data <- rbind(cbind(data_snp, type = "SNPs"), cbind(data_indel, type = "Indels"))
data$type <- factor(data$type, levels = c("SNPs", "Indels"))
plot_sliding_window_mutation_density(data_for_plot = data)
ggsave("motif_mutation_density_sliding_window.pdf", height = 10, width = 20)
