#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:

read.data.bar <- function(filename = NULL, se = NULL) {
    data <- read.table(filename, header = F, as.is = T)
    rownames(data) <- paste("chrom", 1:14, sep = "")
    colnames(data) <- c(
        "motif_3_6_in_codon1", "codon1_count",
        "motif_3_6_in_codon2", "codon2_count",
        "motif_3_6_in_codon3", "codon3_count"
    )

    data <- data[-14, ]
    data <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    if (se == TRUE) {
        result <- t(apply(data, 2, function(x) {
            return(c(mean(x), sd(x)))
        }))
        result <- data.frame(class = rownames(result), MEAN = result[, 1], SD = result[, 2])
        return(result)
    } else {
        return(data)
    }
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/motif_start_site")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

## 柱状图

peak_data <- read.data.bar(filename = "peak_motif_3_6_in_codon123_var", se = FALSE)
nopeak_data <- read.data.bar(filename = "nopeak_motif_3_6_in_codon123_var", se = FALSE)
sapply(1:3, function(i) {
    t.test(peak_data[, i], nopeak_data[, i], paired = TRUE)$p.value
})


ratio <- peak_data / nopeak_data
apply(ratio, 2, function(x) {
    mean(x, na.rm = TRUE)
})
ratio <- data.frame(
    value = unlist(ratio),
    chrom = rep(paste("chrom", 1:13, sep = ""), 3), ## 第14号染色体上没有var基因（非假基因）
    source = rep(c("Codon1", "Codon2", "Codon3"), each = 13)
)
ratio_p <- ggplot(ratio, aes(x = source, y = value, color = source)) +
    geom_boxplot(width = 0.5, size = 1, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(), size = 6) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    scale_color_manual(values = c(
        rgb(92, 167, 147, maxColorValue = 255),
        rgb(92, 167, 147, maxColorValue = 255),
        rgb(92, 167, 147, maxColorValue = 255)
    )) +
    scale_x_discrete(labels = c("First base", "Second base", "Third base")) +
    labs(x = "The position of base in codon", y = "Site proportion ratio\n(6mA vs control)") +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        # panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 36, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        # axis.text.x = element_text(size = 36, color = "black"),
        axis.text.y = element_text(size = 36, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        # legend.position = "bottom", legend.direction = "horizontal",
        legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )




data_for_plot_peak <- read.data.bar(filename = "peak_motif_3_6_in_codon123_var", se = TRUE)
data_for_plot_nopeak <- read.data.bar(filename = "nopeak_motif_3_6_in_codon123_var", se = TRUE)
data_for_plot <- rbind(cbind(data_for_plot_peak, source = "6mA"), cbind(data_for_plot_nopeak, source = "Control"))
data_for_plot$source <- factor(data_for_plot$source, levels = c("6mA", "Control"))

bar_p <- ggplot(data = data_for_plot, aes(x = class, y = MEAN, fill = source)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = MEAN - SD, ymax = MEAN + SD),
        width = .2,
        position = position_dodge(.9)
    ) +
    geom_signif(
        annotations = c("0.31", "0.017", "0.12"), y_position = c(0.007, 0.01, 0.035), xmin = c(0.75, 1.75, 2.75), xmax = c(1.25, 2.25, 3.25),
        tip_length = c(c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01)), vjust = -0.5, size = 1, textsize = 14
    ) +
    scale_y_continuous(limits = c(0, 0.04)) +
    scale_x_discrete(labels = c("First base", "Second base", "Third base")) +
    scale_fill_manual(values = c(
        "6mA" = rgb(236, 174, 34, maxColorValue = 255),
        "Control" = rgb(137, 157, 159, maxColorValue = 255)
    )) +
    labs(x = "", y = "Proportion of motif\n3rd and 6th site") +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        # panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 36, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        # axis.text.x = element_text(size = 36, color = "black"),
        axis.text.y = element_text(size = 36, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_text(size = 40, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, color = "black"),
        # legend.position = "bottom", legend.direction = "horizontal",
        legend.position = c(0.15, 0.9),
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )

pdf("peak_nopeak_3_6_codon123_var.pdf", width = 18, height = 10)
bar_p | ratio_p
grid::grid.draw(grid::textGrob("The position of base in codon", x = 0.6, y = 0.04, gp = grid::gpar(fontsize = 40, col = "black")))
dev.off()
