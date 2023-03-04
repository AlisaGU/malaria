#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggpubr)
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

group_gene <- function(data = NULL, number = NULL) {
    data_sort <- data[order(data$AV, decreasing = TRUE), ]
    data_sort$label <- ceiling((1:nrow(data_sort)) / (nrow(data_sort) / number))
    data_sort$label <- factor(data_sort$label, levels = number:1)
    return(data_sort)
}

plot_boxplot <- function(data_label = NULL, label = NULL) {
    p <- ggplot(data_label) +
        geom_violin(aes(x = label, y = density_signal, color = label), position = position_dodge(), size = 2) +
        geom_boxplot(aes(x = label, y = density_signal, group = label), position = position_dodge(), outlier.shape = NA, alpha = 0.1) +
        geom_point(aes(x = label, y = density_signal, color = label), position = position_jitterdodge(), size = 1) +
        # scale_y_continuous(trans = "log10") +
        scale_color_manual(values = colors) +
        scale_x_discrete(labels = label) +
        labs(y = "6mA density signal", x = "Gene groups with different\nexpression levels") +
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

plot_boxplot1 <- function(data_for_plot = NULL, colors = NULL, labels = NULL) {
    p <- ggplot(data_for_plot, aes(x = label, y = density_signal)) +
        geom_violin(aes(color = label),
            position = position_dodge(), size = 2
        ) +
        geom_boxplot(aes(group = label),
            position = position_dodge(), outlier.shape = NA, alpha = 0.1
        ) +
        geom_point(aes(color = label),
            position = position_jitterdodge(), size = 1
        ) +
        scale_color_manual(values = colors) +
        scale_x_discrete(labels = labels) +
        scale_y_continuous(trans = "log10") +
        labs(y = "6mA density", x = "RNA expression") +
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
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_boxplot2 <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = expression_type, y = density_signal, alpha = exp, color = expression_type)) +
        geom_violin(position = position_dodge(width = 1), size = 1) +
        geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA, width = 0.3) +
        geom_point(position = position_dodge(width = 1)) +
        scale_color_manual(values = rev(c(
            rgb(220, 56, 56, maxColorValue = 255),
            rgb(178, 195, 37, maxColorValue = 255),
            rgb(66, 165, 217, maxColorValue = 255)
        )), guide = NULL) +
        # scale_fill_manual(values=c("white","#7c184f"))+
        scale_x_discrete(labels = c("Down", "No diff", "Up")) +
        scale_y_continuous(trans = "log10") +
        labs(y = "6mA density", x = "RNA expression") +
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
            # legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_boxplot2_WT <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = expression_type, y = density_signal, color = expression_type)) +
        geom_violin(position = position_dodge(), size = 1) +
        geom_boxplot(position = position_dodge(), width = 0.5) +
        geom_point(position = position_jitterdodge()) +
        scale_color_manual(values = rev(c(
            rgb(253, 9, 1, maxColorValue = 255),
            rgb(0, 1, 252, maxColorValue = 255)
        )), guide = NULL) +
        scale_x_discrete(labels = c("Down", "Up")) +
        scale_y_continuous(trans = "log10") +
        labs(y = "6mA density", x = "RNA expression") +
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
            # legend.position = "None",
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
TSS_peakpro_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/TSS_peakpro.txt"
upgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/up_gene"
downgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/down_gene"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
chip_info <- fread(chip_info_filename, stringsAsFactors = F, header = F)
input_info <- fread(input_info_filename, stringsAsFactors = F, header = F)
RNA <- fread(RNA_filename, stringsAsFactors = F, header = T)
TSS_peakpro <- fread(TSS_peakpro_filename, stringsAsFactors = FALSE, header = TRUE)
TSS_peakpro <- TSS_peakpro %>%
    group_by(chrom, TSS_start, TSS_end, gene_name, strand) %>%
    summarise(overlap_base = sum(overlap_base))
colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(RNA)[1] <- "gene"

density_signal_data <- get_6mA(
    all_genes_info = all_genes_info,
    chip_info = chip_info, input_info = input_info,
    chip_depth = chip_depth, input_depth = input_depth
)

data <- merge(density_signal_data, RNA, by = "gene")

TSS_peakpro$pro <- TSS_peakpro$overlap_base / (TSS_peakpro$TSS_end - TSS_peakpro$TSS_start)

data_peakpro <- merge(data, TSS_peakpro, by.x = "gene", by.y = "gene_name")[, c(1, 2, 5, 11)]
## 分成10份
data_label <- group_gene(data = data, number = 10)

colors <- rev(c("#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb", "#d0d1e6", "#ece2f0", "#f9ebf2", "#fcf4f4"))

plot_boxplot(data_label = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)), label = rev(c(paste(paste(seq(0, 90, length.out = 10), c(seq(10, 100, length.out = 10)), sep = "-"), "%", sep = ""))))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_.pdf", width = 8, height = 12)
## peak propor>10%,分成10份
data_label <- group_gene(data = data_peakpro %>% filter(pro >= 0.1), number = 10)

colors <- rev(c("#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb", "#d0d1e6", "#ece2f0", "#f9ebf2", "#fcf4f4"))

plot_boxplot(data_label = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)), label = rev(c(paste(paste(seq(0, 90, length.out = 10), c(seq(10, 100, length.out = 10)), sep = "-"), "%", sep = ""))))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_10part.pdf", width = 8, height = 12)
## peak propor>10%,分成3份(师兄选这个)
data_label <- group_gene(data = data_peakpro %>% filter(pro >= 0.1), number = 3)
t.test(data_label$AV[data_label$label == "1"], data_label$AV[data_label$label == "2"], alternative = "greater")
t.test(data_label$AV[data_label$label == "1"], data_label$AV[data_label$label == "3"], alternative = "greater")
t.test(data_label$AV[data_label$label == "2"], data_label$AV[data_label$label == "3"], alternative = "greater")

colors <- rev(c("#014636", "#67a9cf", "#ece2f0"))

plot_boxplot(data_label = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)), label = rev(c("High", "Medium", "Low")))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_3part.pdf", width = 8, height = 6)

plot_boxplot1(
    data_for_plot = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)), labels = c("Low", "Medium", "High"),
    colors = c(rgb(236, 114, 106, maxColorValue = 255), rgb(157, 43, 127, maxColorValue = 255), rgb(75, 43, 119, maxColorValue = 255))
)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_3part_1.pdf", width = 7, height = 10)

## peak propor>10%,分成2份
data_label <- group_gene(data = data_peakpro %>% filter(pro >= 0.1), number = 2)
t.test(data_label$AV[data_label$label == "1"], data_label$AV[data_label$label == "2"], alternative = "greater")

colors <- rev(c(rgb(75, 43, 119, maxColorValue = 255), rgb(157, 43, 127, maxColorValue = 255)))

plot_boxplot(data_label = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)), label = rev(c("High", "Low")))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_2part.pdf", width = 8, height = 6)

plot_boxplot1(
    data_for_plot = data_label %>% filter(!is.infinite(density_signal)) %>% filter(!is.nan(density_signal)),
    label = rev(c("High", "Low")), colors = colors
)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_2part_1.pdf", width = 8, height = 8)


## 差异表达基因
load("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/density_signal_data_KD.Rdata")
upgene <- as.character(unlist(read.table(upgene_filename, header = F, as.is = T)))
downgene <- as.character(unlist(read.table(downgene_filename, header = F, as.is = T)))
data_peakpro$expression_type <- sapply(data_peakpro$gene, function(x) {
    if (x %in% upgene) {
        return("up")
    } else if (x %in% downgene) {
        return("down")
    } else {
        return("nodiff")
    }
})
data_WT_KD <- merge(data_peakpro, density_signal_data_KD, by = "gene")
WT <- data_WT_KD[, c(1, 2, 3, 4, 5)]
KD <- data_WT_KD[, c(1, 6, 3, 4, 5)]
colnames(KD)[2] <- "density_signal"
data_for_plot <- rbind(data.frame(WT, exp = "WT"), data.frame(KD, exp = "KD"))
data_for_plot <- data_for_plot %>%
    filter(!is.infinite(density_signal)) %>%
    filter(!is.nan(density_signal))
plot_boxplot2(data_for_plot = data_for_plot)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_3part_RNAupDownn.pdf", width = 10, height = 8)

plot_boxplot2_WT(data_for_plot = data_for_plot %>% filter(exp == "WT") %>% filter(expression_type != "nodiff"))
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_pro10Perc_2part_RNAupDownn.pdf", width = 6, height = 8)


WT_down <- data_for_plot %>% filter(exp == "WT" & expression_type == "down")
WT_nodiff <- data_for_plot %>% filter(exp == "WT" & expression_type == "nodiff")
WT_up <- data_for_plot %>% filter(exp == "WT" & expression_type == "up")
KD_down <- data_for_plot %>% filter(exp == "KD" & expression_type == "down")
KD_nodiff <- data_for_plot %>% filter(exp == "KD" & expression_type == "nodiff")
KD_up <- data_for_plot %>% filter(exp == "KD" & expression_type == "up")

t.test(WT_down$density_signal, KD_down$density_signal, alternative = "greater")$p.value
t.test(WT_nodiff$density_signal, KD_nodiff$density_signal, alternative = "greater")$p.value
t.test(WT_up$density_signal, KD_up$density_signal, alternative = "greater")$p.value
