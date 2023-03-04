#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)
# 2. functions ------------------------------------------------------------ TODO:
read.data <- function(part = NULL, type = NULL) {
    chip_filename <- paste0(prefix, ".", part, ".window.", type, ".chip")
    input_filename <- paste0(prefix, ".", part, ".window.", type, ".input")

    chip_data <- fread(chip_filename, header = FALSE, stringsAsFactors = F)
    input_data <- fread(input_filename, header = FALSE, stringsAsFactors = F)
    colnames(chip_data) <- colnames(input_data) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

    chip_read_count <- chip_data %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))

    input_read_count <- input_data %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))
    intersect_gene <- intersect(chip_read_count$gene_name, input_read_count$gene_name)
    chip_read_count_intersect <- chip_read_count[match(intersect_gene, chip_read_count$gene_name), ]
    input_read_count_intersect <- input_read_count[match(intersect_gene, input_read_count$gene_name), ]


    chip_depth <- ""
    input_depth <- ""

    if (type == "wt") {
        chip_depth <- 36971023
        input_depth <- 21500156
    } else if (type == "kd") {
        chip_depth <- 39189761
        input_depth <- 27773718
    }

    if (all(chip_read_count_intersect$gene_name == input_read_count_intersect$gene_name)) {
        density_signal <- (chip_read_count_intersect$depth / chip_depth) / (input_read_count_intersect$depth / input_depth)
        data <- data.frame(gene = chip_read_count_intersect$gene_name, density_signal = density_signal)
        return(data)
    }
}


heatmap <- function(tss_data = NULL, gene_data = NULL, tts_data = NULL) {
    tss_data1 <- convert_data_format(data = tss_data)
    gene_data1 <- convert_data_format(data = gene_data)
    tts_data1 <- convert_data_format(data = tts_data)
    if (all(rownames(tss_data1) == rownames(gene_data1)) & all(rownames(tss_data1) == rownames(tts_data1)) & all(rownames(gene_data1) == rownames(tts_data1))) {
        data <- cbind(tss_data1, gene_data1, tts_data1)
        colnames(data) <- c(paste("tss", 1:100, sep = "_"), paste("gene", 1:100, sep = "_"), paste("tts", 1:100, sep = "_"))
        annotation_col <- data.frame(Part = c(rep("TSS", 100), rep("Gene", 100), rep("TTS", 100)))
        rownames(annotation_col) <- colnames(data)
        bk <- seq(0, 35, length.out = 100)
        color <- colorRampPalette(c(rgb(253, 244, 236, maxColorValue = 255), rgb(128, 40, 4, maxColorValue = 255)))(100)

        p <- pheatmap(data, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation_col, color = color, breaks = bk, fontsize = 20, show_colnames = FALSE)
        return(p)
    }
}

line_data <- function(tss_data = NULL, gene_data = NULL, tts_data = NULL) {
    tss_data1 <- convert_data_format(data = tss_data)
    gene_data1 <- convert_data_format(data = gene_data)
    tts_data1 <- convert_data_format(data = tts_data)
    if (all(rownames(tss_data1) == rownames(gene_data1)) & all(rownames(tss_data1) == rownames(tts_data1)) & all(rownames(gene_data1) == rownames(tts_data1))) {
        data <- cbind(tss_data1, gene_data1, tts_data1)
        colnames(data) <- c(paste("tss", 1:100, sep = "_"), paste("gene", 1:100, sep = "_"), paste("tts", 1:100, sep = "_"))
        value <- apply(data, 2, median)
        return(value)
    }
}


convert_data_format <- function(data = NULL) {
    genename_order <- t(sapply(data$gene, function(x) {
        elements <- unlist(strsplit(x, "_"))
        genename <- paste0(elements[1], "_", elements[2])
        order <- elements[3]
        return(c(genename, order))
    }))
    data$genename <- genename_order[, 1]
    data$order <- genename_order[, 2]
    data$genename <- as.factor(data$genename)
    data <- as.data.table(data)
    data_split <- split(data, f = data$genename, drop = T)
    data_value <- lapply(data_split, function(x) {
        value <- x$density_signal[match(1:100, x$order)]
        return(value)
    })
    data_value <- do.call(rbind, data_value)
    return(data_value)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
prefix <- Args[6]

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/RNA_T_DNA_and_rhoptry")
gene_wt <- read.data(part = "gene", type = "wt")
gene_kd <- read.data(part = "gene", type = "kd")
tss_wt <- read.data(part = "tss", type = "wt")
tss_kd <- read.data(part = "tss", type = "kd")
tts_wt <- read.data(part = "tts", type = "wt")
tts_kd <- read.data(part = "tts", type = "kd")
pdf(file = paste0(prefix, "_wt.pdf"), height = 20, width = 20)
heatmap(tss_data = tss_wt, gene_data = gene_wt, tts_data = tts_wt)
dev.off()

pdf(file = paste0(prefix, "_kd.pdf"), height = 20, width = 20)
heatmap(tss_data = tss_kd, gene_data = gene_kd, tts_data = tts_kd)
dev.off()


wt_line_data <- line_data(tss_data = tss_wt, gene_data = gene_wt, tts_data = tts_wt)
kd_line_data <- line_data(tss_data = tss_kd, gene_data = gene_kd, tts_data = tts_kd)
data_for_plot <- data.frame(group = c(rep("WT", 300), rep("KD", 300)), order = c(1:300, 1:300), value = c(wt_line_data, kd_line_data))
# data_for_plot$order<-factor(data_for_plot$order,1:300)
ggplot(data_for_plot, aes(x = order, y = value, group = group, color = group)) +
    geom_line()
#     geom_point() +
geom_smooth(method = "loess", span = 0.1, se = FALSE) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(
            size = 30, color = "black",
        ),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 36, color = "black"),
        axis.title.x = element_text(size = 36, color = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
