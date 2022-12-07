#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(pheatmap)
library(ggplot2)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/", gene_set, "_peak_control_length_in_gene"), header = F, as.is = T)

        data$SUM <- data[, 2] + data[, 3]
        data <- data.frame(data[, 1], data[, 2] / data$SUM, data[, 3] / data$SUM, gene_set)
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "peak_proportion", "control_proportion", "gene_set_type")
    result <- as.data.frame(result)
    return(result)
}

read_data_for_ggplot <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/", gene_set, "_peak_control_length_in_gene"), header = F, as.is = T)

        data$SUM <- data[, 2] + data[, 3]
        data <- data.frame(rep(data[, 1], 2), rbind(
            data.frame(values = data[, 2] / data$SUM, peak_type = "peak_proportion"),
            data.frame(values = data[, 3] / data$SUM, peak_type = "control_proportion")
        ), gene_set)
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "proportion", "peak_type", "gene_set_type")
    result <- as.data.frame(result)
    return(result)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_exclude_pseudo")
gene_sets <- c("DR", "HDR", "RNA_translation", "STEVOR", "VAR", "RIF")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read_data()
sapply(split(data, f = data$gene_set_type), function(x) {
    length(which(x[, 2] >= 0.2)) / nrow(x)
})
data %>%
    group_by(gene_set_type) %>%
    summarise()
data_for_plot <- as.matrix(data[, 2:3])
rownames(data_for_plot) <- data[, 1]
annotation_row <- data.frame(Gene_set = c(
    rep("DR", 7),
    rep("HDR", 22),
    rep("RNA_translation", 44),
    rep("STEVOR", 32),
    rep("VAR", 65),
    rep("RIF", 157)
))
rownames(annotation_row) <- data[, 1]

bk <- seq(0, 1, length.out = 100)
# color <- c(colorRampPalette(c("blue", "white"))(50), colorRampPalette(c("white", "red"))(50))
pdf("peak_control_ratio.pdf", width = 15, height = 5)
pheatmap(t(data_for_plot), annotation_col = annotation_row, show_colnames = F, breaks = bk, color = colorRampPalette(c("white", "red"))(100), cluster_rows = F)
dev.off()

ggplot_data <- read_data_for_ggplot()
ggplot_data$gene_set_type <- factor(ggplot_data$gene_set_type,
    levels = c("RIF", "STEVOR", "VAR", "HDR", "RNA_translation")
)
ggplot(ggplot_data %>%
    filter(peak_type == "peak_proportion") %>%
    filter(gene_set_type != "DR"), aes(x = gene_set_type, y = proportion, color = gene_set_type)) +
    geom_violin() +
    scale_color_manual(
        values = c(
            "RIF" = rgb(148, 84, 151, maxColorValue = 255),
            "STEVOR" = rgb(113, 184, 196, maxColorValue = 255),
            "VAR" = rgb(114, 172, 69, maxColorValue = 255),
            "HDR" = rgb(192, 4, 54, maxColorValue = 255),
            "RNA_translation" = rgb(189, 151, 62, maxColorValue = 255)
        )
    ) +
    scale_x_discrete(labels = c("RIF", "STEVOR", "VAR", "HDR", "RNA\ntranslation")) +
    labs(y = "Proportion of 6mA peak") +
    geom_point(position = position_jitterdodge(), size = 6) +
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
            size = 36, color = "black"
            # angle = 45, hjust = 0.5, vjust = 0.5
        ),
        axis.text.y = element_text(size = 36, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave("proportion_of_peak_in_each_gene.pdf", width = 10, height = 8)
