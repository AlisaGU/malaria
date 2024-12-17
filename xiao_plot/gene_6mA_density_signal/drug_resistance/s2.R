#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))

# 2. functions ------------------------------------------------------------ TODO:
read.data <- function(direction = NULL, part = NULL) {
    chip_filename <- paste0(prefix, ".", direction, ".", part, ".window.wt.chip")
    input_filename <- paste0(prefix, ".", direction, ".", part, ".window.wt.input")

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


    chip_depth <- 36971023
    input_depth <- 21500156

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

    common_gene <- intersect(intersect(rownames(tss_data1), rownames(gene_data1)), rownames(tts_data1))
    tss_data1 <- tss_data1[match(common_gene, rownames(tss_data1)), ]
    gene_data1 <- gene_data1[match(common_gene, rownames(gene_data1)), ]
    tts_data1 <- tts_data1[match(common_gene, rownames(tts_data1)), ]

    data <- cbind(tss_data1, gene_data1, tts_data1)
    colnames(data) <- c(paste("tss", 1:100, sep = "_"), paste("gene", 1:300, sep = "_"), paste("tts", 1:100, sep = "_"))
    value <- apply(data, 2, function(x) {
        mean(x[which(!is.nan(x) & !is.infinite(x))])
    })
    return(value)
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
working_dir <- Args[6]
prefix <- Args[7]

# 4. variable setting of test module--------------------------------------- TODO:
# working_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_ZbynekBozdech_Rstage"
# prefix <- "s16_S1A_6AR_vs_6AS"

# 5. process -------------------------------------------------------------- TODO:
setwd(working_dir)
gene_up <- read.data(part = "gene", direction = "up")
gene_down <- read.data(part = "gene", direction = "down")

tss_up <- read.data(part = "tss", direction = "up")
tss_down <- read.data(part = "tss", direction = "down")

tts_up <- read.data(part = "tts", direction = "up")
tts_down <- read.data(part = "tts", direction = "down")

up_line_data <- line_data(tss_data = tss_up, gene_data = gene_up, tts_data = tts_up)
down_line_data <- line_data(tss_data = tss_down, gene_data = gene_down, tts_data = tts_down)

data_for_plot <- data.frame(group = c(rep("Up", 500), rep("Down", 500)), order = c(1:500, 1:500), value = c(up_line_data, down_line_data))
# data_for_plot$order<-factor(data_for_plot$order,1:300)
ggplot(data_for_plot, aes(x = order, y = value, group = group, color = group)) +
    # geom_line()
    # geom_point() +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, size = 3) +
    scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
    scale_color_manual(values = c("Up" = "red", "Down" = "#0000fc"), ) +
    labs(y = "6mA signal density") +
    theme_classic() +
    theme(
        axis.ticks.length = unit(.25, "cm"),
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
        legend.position = c(0.15, 0.9),
        legend.text = element_text(size = 36, color = "black"),
        legend.key.width = unit(2, "cm"),
        plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave(paste0(prefix, ".pdf"), width = 10, height = 7)
