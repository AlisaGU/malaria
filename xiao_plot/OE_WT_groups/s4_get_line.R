#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))

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
        chip_depth <- 35068263
        input_depth <- 27989115
    } else if (type == "kd") {
        chip_depth <- 32038632
        input_depth <- 8218045
    }

    if (all(chip_read_count_intersect$gene_name == input_read_count_intersect$gene_name)) {
        density_signal <- (chip_read_count_intersect$depth / chip_depth) / (input_read_count_intersect$depth / input_depth)
        data <- data.frame(gene = chip_read_count_intersect$gene_name, density_signal = density_signal)
        return(data)
    }
}

line_data <- function(tss_data = NULL, gene_data = NULL, tts_data = NULL) {
    tss_data1 <- convert_data_format(data = tss_data)
    gene_data1 <- convert_data_format(data = gene_data)
    tts_data1 <- convert_data_format(data = tts_data)

    common_gene <- intersect(rownames(tts_data1), intersect(rownames(tss_data1), rownames(gene_data1)))
    tss_data1 <- tss_data1[match(common_gene, rownames(tss_data1)), ]
    gene_data1 <- gene_data1[match(common_gene, rownames(gene_data1)), ]
    tts_data1 <- tts_data1[match(common_gene, rownames(tts_data1)), ]


    if (all(rownames(tss_data1) == rownames(gene_data1)) & all(rownames(tss_data1) == rownames(tts_data1)) & all(rownames(gene_data1) == rownames(tts_data1))) {
        data <- cbind(tss_data1, gene_data1, tts_data1)
        colnames(data) <- c(paste("tss", 1:100, sep = "_"), paste("gene", 1:300, sep = "_"), paste("tts", 1:100, sep = "_"))
        value <- apply(data, 2, function(x) {
            mean(x[which(!is.nan(x) & !is.infinite(x))])
        })
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
        value <- x$density_signal[match(1:max(as.numeric(unlist(x$order))), x$order)]
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
setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/", prefix))
gene_wt <- read.data(part = "gene", type = "wt")
gene_kd <- read.data(part = "gene", type = "kd")
tss_wt <- read.data(part = "tss", type = "wt")
tss_kd <- read.data(part = "tss", type = "kd")
tts_wt <- read.data(part = "tts", type = "wt")
tts_kd <- read.data(part = "tts", type = "kd")


wt_line_data <- line_data(tss_data = tss_wt, gene_data = gene_wt, tts_data = tts_wt)
kd_line_data <- line_data(tss_data = tss_kd, gene_data = gene_kd, tts_data = tts_kd)
data_for_plot <- data.frame(group = c(rep("WT", 500), rep("OE", 500)), order = c(1:500, 1:500), value = c(wt_line_data, kd_line_data), prefix = prefix)
save(data_for_plot, file = paste0(prefix, ".Rdata"))
