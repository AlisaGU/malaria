#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(all_genes_info = NULL, chip_info = NULL, input_info = NULL, chip_depth = NULL, input_depth = NULL) {
    chip_read_count <- chip_info %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))

    input_read_count <- input_info %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))

    intersect_gene <- intersect(chip_read_count$gene_name, input_read_count$gene_name)
    chip_read_count_intersect <- chip_read_count[match(intersect_gene, chip_read_count$gene_name), ]
    input_read_count_intersect <- input_read_count[match(intersect_gene, input_read_count$gene_name), ]

    if (all(chip_read_count_intersect$gene_name == input_read_count_intersect$gene_name)) {
        density_signal <- (chip_read_count_intersect$depth / chip_depth) / (input_read_count_intersect$depth / input_depth)
        data <- data.frame(gene = chip_read_count_intersect$gene_name, density_signal = density_signal)
        return(data)
    }
}

group_gene <- function(data = NULL) {
    data_sort <- data[order(data$AV, decreasing = TRUE), ]
    data_sort$label <- ceiling((1:nrow(data_sort)) / (nrow(data_sort) / 10))
    data_sort$label <- factor(data_sort$label, levels = 10:1)
    return(data_sort)
}

plot_boxplot <- function(data_label = NULL) {
    p <- ggplot(data_label) +
        geom_violin(aes(x = label, y = density_signal, color = label), position = position_dodge(), size = 2) +
        geom_boxplot(aes(x = label, y = density_signal, group = label), position = position_dodge(), outlier.shape = NA, alpha = 0.1) +
        geom_point(aes(x = label, y = density_signal, color = label), position = position_jitterdodge(), size = 1) +
        scale_y_continuous(trans = "log10") +
        scale_color_manual(values = colors) +
        scale_x_discrete(labels = rev(c(paste(paste(seq(0, 90, length.out = 10), c(seq(10, 100, length.out = 10)), sep = "-"), "%", sep = "")))) +
        labs(y = "6mA density signal", x = "Gene groups with different expression levels") +
        theme_bw() +
        coord_flip() +
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
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_ridge <- function(data_label = NULL) {
    p <- ggplot(data_label[!is.na(data_label$density_signal) & !is.infinite(data_label$density_signal) & data_label$density_signal != 0, ], aes(y = label, x = density_signal)) +
        geom_density_ridges2(aes(fill = label), color = "black") +
        scale_x_continuous(trans = "log10") +
        scale_y_discrete(labels = rev(c(paste(paste(seq(0, 90, length.out = 10), c(seq(10, 100, length.out = 10)), sep = "-"), "%", sep = "")))) +
        scale_fill_manual(values = colors) +
        labs(x = "6mA density signal", y = "Gene groups with different expression levels") +
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
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/all_genes_bed"
chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"
chip_depth <- 36971023
input_depth <- 21500156
RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_T_RNA expression.txt"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
chip_info <- fread(chip_info_filename, stringsAsFactors = F, header = F)
input_info <- fread(input_info_filename, stringsAsFactors = F, header = F)
RNA <- fread(RNA_filename, stringsAsFactors = F, header = T)

colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(RNA)[1] <- "gene"

density_signal_data <- get_6mA(
    all_genes_info = all_genes_info,
    chip_info = chip_info, input_info = input_info,
    chip_depth = chip_depth, input_depth = input_depth
)

data <- merge(density_signal_data, RNA, by = "gene")
data_label <- group_gene(data = data)

colors <- rev(c("#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb", "#d0d1e6", "#ece2f0", "#f9ebf2", "#fcf4f4"))

plot_boxplot(data_label = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA.pdf", width = 8, height = 12)
