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

# 2. functions ------------------------------------------------------------ TODO:
read_sliding_window_data <- function(dataname = NULL) {
    data <- read.table(dataname, as.is = T, header = F)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    rownames(data1) <- paste("chrom", 1:14, sep = "")
    colnames(data1) <- colnames(data) <- c("nopeak", "nopeak500", paste("w", 1:10, sep = ""))
    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c("nopeak", "nopeak500", paste("w", 1:10, sep = "")))
    return(data_for_plot)
}


plot_and_cor_sliding_window <- function(data_for_plot = NULL, title = NULL) {
    cols <- brewer.pal(3, "YlOrRd")
    pal <- colorRampPalette(cols)
    mycolors <- rev(pal(10))
    color_value <- c("#000000", "#000000", mycolors)
    names(color_value) <- c("nopeak", "nopeak500", paste("w", 1:10, sep = ""))

    data_for_plot$distance <- as.integer(data_for_plot$source) - 2
    data_for_plot$distance[data_for_plot$distance == 0] <- 11
    data_for_plot$distance[data_for_plot$distance == -1] <- 12

    p <- ggplot(data_for_plot) +
        geom_boxplot(aes(x = distance, y = value, color = source),
            position = position_dodge(), outlier.shape = NA
        ) +
        geom_point(aes(x = distance, y = value, color = source),
            position = position_jitterdodge(), size = 4
        ) +
        # geom_smooth(aes(x = distance, y = value),
        #     method = "lm", se = TRUE, group = 1,
        #     fullrange = FALSE, level = 0.95
        # ) +
        # stat_cor(aes(x = distance, y = value), ,
        #     size = 20, label.x = 3.5, label.y = 0.0045
        # ) +
        scale_color_manual(values = color_value) +
        scale_x_continuous(breaks = c(1:12), labels = c("Submit~5%", "5~10%", "10~15%", "15~20%", "20~25%", "25~30%", "30~35%", "35~40%", "40~45%", "45~50%", "Control", "Control(500)")) +
        labs(y = "6mA sites density") +
        ggtitle(title) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 30, color = "black",
                angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 30, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

read_motif_flank_data <- function(dataname = NULL) {
    data <- read.table(dataname, as.is = T, header = F)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    rownames(data1) <- paste("chrom", 1:14, sep = "")
    colnames(data1) <- c(paste("site", 1:6, sep = ""), paste("flank", 1:10, sep = ""))
    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c(paste("site", 1:6, sep = ""), paste("flank", 1:10, sep = "")))
    data_for_plot$class <- NA
    data_for_plot$class[grep("site", data_for_plot$source)] <- "motif"
    data_for_plot$class[grep("flank", data_for_plot$source)] <- "motif_flank"
    return(data_for_plot)
}

plot_motif_flank <- function(data_for_plot = NULL, title = NULL) {
    color_value <- c("motif" = "#b2182b", "motif_flank" = "#2166ac")
    p <- ggplot(data_for_plot) +
        geom_boxplot(aes(x = source, y = value, color = class),
            position = position_dodge(), outlier.shape = NA
        ) +
        geom_point(aes(x = source, y = value, color = class),
            position = position_jitterdodge(), size = 4
        ) +
        stat_summary(aes(x = source, y = value, color = class),
            fun = median, colour = "#fabf12",
            geom = "point", position = position_dodge(width = 0.75)
        ) +
        stat_summary(aes(x = source, y = value, color = class, group = 1),
            fun = median, colour = "#fabf12",
            geom = "line", lwd = 1, lty = 1
        ) +
        labs(y = "6mA sites density") +
        ggtitle(title) +
        scale_color_manual(values = color_value) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 30, color = "black",
                angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 30, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}


read_motif_and_averaged_flank_data <- function(dataname = NULL) {
    data <- read.table(dataname, as.is = T, header = F)
    data_motif <- data[, 1:12]
    data_flank <- data[, -c(1:12)]
    data1 <- data_motif[, seq(1, ncol(data_motif), by = 2)] / data_motif[, seq(2, ncol(data_motif), by = 2)]
    data2 <- apply(data_flank[, seq(1, ncol(data_flank), by = 2)], 1, sum) / apply(data_flank[, seq(2, ncol(data_flank), by = 2)], 1, sum)

    data1 <- cbind(data1, matrix(data2, ncol = 1))
    rownames(data1) <- paste("chrom", 1:14, sep = "")
    colnames(data1) <- c(paste("site", 1:6, sep = ""), "flank1-10")
    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c(paste("site", 1:6, sep = ""), "flank1-10"))
    data_for_plot$class <- NA
    data_for_plot$class[grep("site", data_for_plot$source)] <- "motif"
    data_for_plot$class[grep("flank", data_for_plot$source)] <- "motif_flank"
    return(data_for_plot)
}

plot_motif_averaged_flank <- function(data_for_plot = NULL, title = NULL) {
    color_value <- c("motif" = "#b2182b", "motif_flank" = "#2166ac")
    p <- ggplot(data_for_plot) +
        geom_boxplot(aes(x = source, y = value, color = class),
            position = position_dodge(), outlier.shape = NA
        ) +
        geom_point(aes(x = source, y = value, color = class),
            position = position_jitterdodge(), size = 4
        ) +
        stat_summary(aes(x = source, y = value, color = class),
            fun = median, colour = "#fabf12",
            geom = "point", position = position_dodge(width = 0.75)
        ) +
        stat_summary(aes(x = source, y = value, color = class, group = 1),
            fun = median, colour = "#fabf12",
            geom = "line", lwd = 1, lty = 1
        ) +
        labs(y = "6mA sites density") +
        ggtitle(title) +
        scale_color_manual(values = color_value) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 30, color = "black"
            ),
            axis.text.y = element_text(size = 30, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}


# 3. input ---------------------------------------------------------------- TODO:
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/sliding_window")
pdf("6mA_distri_sliding_window.pdf", width = 14, height = 7)
for (prefix in c("", "_COVgt25", "_COVgt50", "_QVgt20_COVgt10", "_FRACgt20_COVgt25")) {
    data_for_plot <- read_sliding_window_data(dataname = paste0("sliding_window_6mA_distri", prefix))
    p <- plot_and_cor_sliding_window(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/sliding_window")
pdf("6mA_distri_sliding_window_ATcountScaled.pdf", width = 14, height = 7)
for (prefix in c("", "_COVgt25", "_COVgt50", "_QVgt20_COVgt10", "_FRACgt20_COVgt25")) {
    data_for_plot <- read_sliding_window_data(dataname = paste0("sliding_window_6mA_distri_ATcountScaled", prefix))
    p <- plot_and_cor_sliding_window(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern")
pdf("6mA_distri_motif_flank.pdf", width = 14, height = 7)
# for (prefix in c("", "_COVgt25", "_COVgt50", "_QVgt20_COVgt10", "_FRACgt20_COVgt25")) {
for (prefix in c("_FRACgt20_COVgt25")) {
    data_for_plot <- read_motif_flank_data(dataname = paste0("motif_flank_pattern_6mA_distri", prefix))
    p <- plot_motif_flank(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern")
pdf("6mA_distri_motif_flank_ATcountScaled.pdf", width = 14, height = 7)
for (prefix in c("", "_COVgt25", "_COVgt50", "_QVgt20_COVgt10", "_FRACgt20_COVgt25")) {
    data_for_plot <- read_motif_flank_data(dataname = paste0("motif_flank_pattern_6mA_distri_ATcountScaled", prefix))
    p <- plot_motif_flank(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern")
pdf("6mA_distri_motif_averaged_flank.pdf", width = 14, height = 7)
for (prefix in c("_FRACgt20_COVgt25")) {
    data_for_plot <- read_motif_and_averaged_flank_data(dataname = paste0("motif_flank_pattern_6mA_distri", prefix))
    p <- plot_motif_averaged_flank(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()








setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/private/sliding_window")
pdf("6mA_distri_sliding_window.pdf", width = 14, height = 7)
for (prefix in c("")) {
    data_for_plot <- read_sliding_window_data(dataname = paste0("sliding_window_6mA_distri", prefix))
    p <- plot_and_cor_sliding_window(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/private/sliding_window")
pdf("6mA_distri_sliding_window_ATcountScaled.pdf", width = 14, height = 7)
for (prefix in c("")) {
    data_for_plot <- read_sliding_window_data(dataname = paste0("sliding_window_6mA_distri_ATcountScaled", prefix))
    p <- plot_and_cor_sliding_window(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/private/motif_flank_pattern")
pdf("6mA_distri_motif_flank.pdf", width = 14, height = 7)
for (prefix in c("")) {
    data_for_plot <- read_motif_flank_data(dataname = paste0("motif_flank_pattern_6mA_distri", prefix))
    p <- plot_motif_flank(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/private/motif_flank_pattern")
pdf("6mA_distri_motif_flank_ATcountScaled.pdf", width = 14, height = 7)
for (prefix in c("")) {
    data_for_plot <- read_motif_flank_data(dataname = paste0("motif_flank_pattern_6mA_distri_ATcountScaled", prefix))
    p <- plot_motif_flank(data_for_plot = data_for_plot, title = prefix)
    print(p)
}
dev.off()
