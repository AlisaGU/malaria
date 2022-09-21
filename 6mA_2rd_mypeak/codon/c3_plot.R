#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)

# 2. functions ------------------------------------------------------------ TODO:
read.data.boxplot <- function(filename = NULL) {
    data <- read.table(filename, header = F, as.is = T)
    rownames(data) <- paste("chrom", 1:14, sep = "")
    colnames(data) <- c("1st Base", "2nd Base", "3rd Base")

    data <- data / apply(data, 1, sum)
    result <- data.frame(value = unlist(data), chrom = rep(rownames(data), ncol(data)), source = rep(colnames(data), each = nrow(data)))
    result$chrom <- factor(result$chrom, levels = paste("chrom", 1:14, sep = ""))
    return(result)
}
read.data.donut <- function(filename = NULL) {
    data <- read.table(filename, header = F, as.is = T)
    rownames(data) <- paste("chrom", 1:14, sep = "")
    colnames(data) <- c("1st Base", "2nd Base", "3rd Base")
    data <- apply(data, 2, sum) / sum(apply(data, 2, sum))
    result <- data.frame(category = names(data), fraction = data)
    return(result)
}

read.data.bar <- function(filename = NULL) {
    data <- read.table(filename, header = F, as.is = T)
    rownames(data) <- paste("chrom", 1:14, sep = "")
    colnames(data) <- c("1st Base", "2nd Base", "3rd Base")
    result <- t(apply(data, 2, function(x) {
        return(c(mean(x), sd(x)))
    }))
    result <- data.frame(class = rownames(result), MEAN = result[, 1], SD = result[, 2])
    return(result)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/motif_start_site")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
## 箱式图
data_for_plot <- read.data.boxplot(filename = "motif_each_site")
ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
    geom_boxplot(aes(color = source), width = 0.5, size = 2) +
    geom_point(aes(color = source), position = position_jitterdodge(), size = 6) +
    stat_compare_means(
        comparisons = list(c("1st Base", "2nd Base"), c("1st Base", "3rd Base"), c("2nd Base", "3rd Base")),
        label = "p.signif", paired = T, size = 10, tip.length = 0.01, bracket.size = 0.1
    ) +
    scale_color_manual(values = c(
        "1st Base" = rgb(241, 213, 124, maxColorValue = 255),
        "2nd Base" = rgb(234, 130, 19, maxColorValue = 255),
        "3rd Base" = rgb(47, 0, 53, maxColorValue = 255)
    )) +
    scale_x_discrete(labels = c("1st Base", "2nd Base", "3rd Base")) +
    labs(x = "", y = "Motif first base occurrence") +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        panel.border = element_blank(),
        # panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
        # axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_text(size = 40, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        legend.position = "bottom", legend.direction = "horizontal",
        # legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )

## 甜甜圈图
data_for_plot <- read.data.donut(filename = filename)
data_for_plot$ymax <- cumsum(data_for_plot$fraction)
data_for_plot$ymin <- c(0, head(data_for_plot$ymax, n = -1))
data_for_plot$labelPosition <- (data_for_plot$ymax + data_for_plot$ymin) / 2
data_for_plot$fraction1 <- scales::percent(data_for_plot$fraction)

ggplot(data_for_plot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = category)) +
    geom_rect() +
    geom_text(x = 3.5, aes(y = labelPosition, label = fraction1), size = 10) +
    scale_fill_manual(values = c(
        "1st Base" = rgb(54, 69, 155, maxColorValue = 255),
        "2nd Base" = rgb(177, 223, 229, maxColorValue = 255),
        "3rd Base" = rgb(250, 163, 27, maxColorValue = 255)
    )) +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = c(0.5, 0.5), legend.title = element_blank(), legend.text = element_text(size = 20))
ggsave("nopeak_motif_start_site.pdf", height = 7, width = 7)
## 柱状图
data_for_plot <- read.data.bar(filename = filename)
ggplot(data_for_plot, aes(x = class, y = MEAN, fill = class)) +
    geom_bar(
        stat = "identity", color = "black",
        position = position_dodge()
    ) +
    geom_errorbar(aes(ymin = MEAN - SD, ymax = MEAN + SD),
        width = .2,
        position = position_dodge(.9)
    ) +
    scale_fill_manual(values = c(
        "1st Base" = rgb(54, 69, 155, maxColorValue = 255),
        "2nd Base" = rgb(177, 223, 229, maxColorValue = 255),
        "3rd Base" = rgb(250, 163, 27, maxColorValue = 255)
    )) +
    labs(x = "", y = "Motif first base occurrence") +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        panel.border = element_blank(),
        # panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
        # axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_text(size = 40, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        legend.position = "bottom", legend.direction = "horizontal",
        # legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
