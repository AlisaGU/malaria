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


get_data_for_plot_facet <- function(dataname = NULL) {
    data <- fread(dataname, stringsAsFactors = F, sep = " ", select = 1:12)
    colnames(data) <- c("nopeak", "nopeak500", paste("w", 1:10, sep = ""))
    rownames(data) <- paste("chrom", 1:14, sep = "")

    data_for_plot <- data.frame(
        value = unlist(data),
        chrom = rep(rownames(data), times = ncol(data)),
        source = rep(colnames(data), each = nrow(data))
    )
    a <- gsub("chrom", "", data_for_plot$chrom)
    a[a != "1"] <- ""
    data_for_plot$label <- a
    data_for_plot$histone <- unlist(strsplit(dataname, "/"))[1]
    # data_for_plot$variant <- gsub("_mean_variant_count", "", unlist(strsplit(dataname, "/"))[length(unlist(strsplit(dataname, "/")))])
    # data_for_plot$variant <- gsub("_mean_variant_count_6mA", "", unlist(strsplit(dataname, "/"))[length(unlist(strsplit(dataname, "/")))])

    data_for_plot$variant <- gsub("_mean_variant_count_noK36me3", "", unlist(strsplit(dataname, "/"))[length(unlist(strsplit(dataname, "/")))])
    # data_for_plot$variant <- gsub("_mean_variant_count_noncore", "", unlist(strsplit(dataname, "/"))[length(unlist(strsplit(dataname, "/")))])

    return(data_for_plot)
}


plot_mutation_density <- function(data_for_plot = NULL, histone_type = NULL) {
    # stat.test <- compare_means(
    #     value ~ source,
    #     data = subset(data_for_plot, source == "nopeak" | source == "w1"),
    #     method = "t.test", paired = T
    # )
    # stat.test <- stat.test %>% mutate(y.position = 0.15)

    p <- ggplot(data_for_plot %>% filter(histone == histone_type), aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        facet_wrap(~variant_f, nrow = 2, ncol = 1, scale = "free") +
        # stat_pvalue_manual(data = stat.test, label = "p = {p.format}", size = 10, remove.bracket = T) +
        # stat_compare_means(
        #     method = "t.test", tip.length = 0, bracket.size = 0,
        #     comparisons = list(c("nopeak", "w1")),
        #     hide.ns = T, vjust = 0.5, show.legend = FALSE,
        #     label = "p = {p.format}", paired = T, size = 10
        # ) +
        ggtitle(histone_type) +
        scale_color_manual(values = color_value) +
        scale_x_discrete(labels = c("Control", "Control(500)", "Submit~5%", "5~10%", "10~15%", "15~20%", "20~25%", "25~30%", "30~35%", "35~40%", "40~45%", "45~50%")) +
        labs(y = "Mutation density") +
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


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
path <- Args[6]
popu_symbol <- Args[7]
# 4. variable setting of test module--------------------------------------- TODO:
# path<-"/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS"
# popu_symbol<-"ESEA_WSEA_OCE_SAM_SAS"
# 5. process -------------------------------------------------------------- TODO:
setwd(path)


cols <- brewer.pal(3, "YlOrRd")
pal <- colorRampPalette(cols)
mycolors <- rev(pal(10))
color_value <- c("#868686", mycolors)
names(color_value) <- c("control", paste("w", 1:10, sep = ""))
color_value_for_control <- c("nopeak" = "#868686", "nopeak_motif" = "#868686", "nopeak500_motif" = "#2f2f2f")
color_value_for_peak <- c("peak" = "#F03B20", "peak_motif" = "#F03B20")

variants_type <- c("snps", "indels")
histone_type <- c("6mA", "H3K36me3", "H3K9Ac", "H3K9me3", "H3K4me3")
# data_for_plot <- lapply(histone_type, function(histone) {
#     result <- lapply(variants_type, function(x) {
#         get_data_for_plot_facet(paste0(histone, "/sliding_window/relative.0.05/mutation_density/", x, "_mean_variant_count"))
#     })
#     return(do.call(rbind, result))
# })
# names(data_for_plot) <- histone_type

data_for_plot <- lapply("6mA", function(histone) {
    result <- lapply(variants_type, function(x) {
        get_data_for_plot_facet(paste0(histone, "/sliding_window/relative.0.05/mutation_density/", x, "_mean_variant_count_noK36me3"))
    })
    return(do.call(rbind, result))
})
data_for_plot_global <- do.call(rbind, data_for_plot)
data_for_plot_global$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", data_for_plot_global$variant))
data_for_plot_global$variant_f <- factor(data_for_plot_global$variant, levels = c("SNPs", "Indels"))
data_for_plot_global$source <- factor(data_for_plot_global$source, levels = c("nopeak", "nopeak500", paste("w", 1:10, sep = "")))

pdf(file = "histone_sliding_window.pdf", width = 20, height = 15)
for (histone in c("6mA", "H3K36me3", "H3K9Ac", "H3K9me3", "H3K4me3")) {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_no6mA.pdf", width = 15, height = 15)
for (histone_type in histone) {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_core.pdf", width = 15, height = 15)
for (histone_type in "H3K36me3") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_core_no6ma.pdf", width = 15, height = 15)
for (histone_type in "H3K36me3") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_no6ma_extend.pdf", width = 15, height = 15)
for (histone_type in "H3K36me3") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_noncore.pdf", width = 15, height = 15)
for (histone_type in "H3K36me3") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_6mA.pdf", width = 15, height = 15)
for (histone_type in "H3K36me3") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()

pdf(file = "histone_sliding_window_noK36me3.pdf", width = 15, height = 15)
for (histone_type in "6mA") {
    p <- plot_mutation_density(data_for_plot = data_for_plot_global, histone_type = histone_type)
    print(p)
}
dev.off()
