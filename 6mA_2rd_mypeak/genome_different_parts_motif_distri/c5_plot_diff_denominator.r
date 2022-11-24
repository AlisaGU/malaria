#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read.data <- function(filename = NULL) {
    data <- read.table(filename, stringsAsFactors = F, as.is = T)
    rownames(data) <- paste("chrom", 1:14, sep = "")
    colnames(data) <- paste(rep(c("intergenic", "TSS_2KB", "TTS_2KB", "genes", "exons", "introns"), each = 2), rep(c("motif_base_count", "region_len"), 6), sep = ".")
    motif_base_count <- data[, seq(1, ncol(data), by = 2)]
    region_len <- data[, seq(2, ncol(data), by = 2)]

    motif_avergae_base_count <- motif_base_count / region_len

    b <- data.frame(
        value = unlist(motif_avergae_base_count),
        chrom = rep(rownames(motif_avergae_base_count), ncol(motif_avergae_base_count)),
        source = rep(gsub(".motif_base_count", "", colnames(motif_avergae_base_count)),
            each = nrow(motif_avergae_base_count)
        ),
        class = "motif_base_count"
    )

    result <- as.data.frame(b)
    result$source <- factor(result$source, levels = c("intergenic", "TSS_2KB", "TTS_2KB", "genes", "exons", "introns"))
    return(result)
}

get_plot <- function(filename = NULL) {
    data_for_plot <- read.data(filename = filename)
    p <- ggplot(data_for_plot %>% filter(class == "motif_base_count") %>% filter(source != "intergenic") %>% filter(source != "genes"), aes(x = source, y = value)) +
        geom_boxplot(aes(color = source), width = 0.5, size = 1) +
        geom_point(aes(color = source), position = position_jitterdodge(), size = 4) +
        # stat_compare_means(
        #     comparisons = list(c("intergenic", "exons"), c("TSS_2KB", "exons"), c("TTS_2KB", "exons"), c("genes", "exons"), c("introns", "exons")),
        #     label = "p.signif", paired = T, size = 10, tip.length = 0.01, bracket.size = 0.1
        # ) +
        # facet_wrap(~class, scales = "free") +
        scale_color_manual(values = c(
            "intergenic" = rgb(96, 127, 155, maxColorValue = 255),
            "TSS_2KB" = rgb(125, 175, 76, maxColorValue = 255),
            "TTS_2KB" = rgb(91, 105, 79, maxColorValue = 255),
            "genes" = rgb(121, 84, 153, maxColorValue = 255),
            "exons" = rgb(192, 127, 59, maxColorValue = 255),
            "introns" = rgb(162, 79, 49, maxColorValue = 255)
        )) +
        # scale_x_discrete(labels = c("Intergenic", "TSS 2kb", "TTS 2kb", "Gene", "Exon", "Intron")) +
        # scale_x_discrete(labels = c("TSS 2kb", "TTS 2kb", "Gene", "Exon", "Intron")) +
        scale_x_discrete(labels = c("TSS 2kb", "TTS 2kb", "Exon", "Intron")) +
        labs(x = "", y = "Proportion of motif base") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
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
            # legend.position = "bottom", legend.direction = "horizontal",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

get_ratio_plot <- function() {
    data_for_plot_peak <- read.data(filename = "peak_motif_density")
    data_for_plot_nopeak <- read.data(filename = "nopeak_motif_density")
    data_for_plot_ratio <- data.frame(
        value = data_for_plot_peak$value / data_for_plot_nopeak$value,
        chrom = data_for_plot_peak$chrom, source = data_for_plot_peak$source,
        class = data_for_plot_peak$class
    )

    p <- ggplot(data_for_plot_ratio %>% filter(class == "motif_base_count") %>% filter(source != "intergenic") %>% filter(source != "genes"), aes(x = source, y = value)) +
        geom_boxplot(aes(color = source), width = 0.5, size = 1) +
        geom_point(aes(color = source), position = position_jitterdodge(), size = 4) +
        geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
        # stat_compare_means(
        #     comparisons = list(c("intergenic", "exons"), c("TSS_2KB", "exons"), c("TTS_2KB", "exons"), c("genes", "exons"), c("introns", "exons")),
        #     label = "p.signif", paired = T, size = 10, tip.length = 0.01, bracket.size = 0.1
        # ) +
        # facet_wrap(~class, scales = "free") +
        scale_color_manual(values = c(
            "intergenic" = rgb(96, 127, 155, maxColorValue = 255),
            "TSS_2KB" = rgb(125, 175, 76, maxColorValue = 255),
            "TTS_2KB" = rgb(91, 105, 79, maxColorValue = 255),
            "genes" = rgb(121, 84, 153, maxColorValue = 255),
            "exons" = rgb(192, 127, 59, maxColorValue = 255),
            "introns" = "black"
        )) +
        # scale_x_discrete(labels = c("Intergenic", "TSS 2kb", "TTS 2kb", "Gene", "Exon", "Intron")) +
        # scale_x_discrete(labels = c("TSS 2kb", "TTS 2kb", "Gene", "Exon", "Intron")) +
        scale_x_discrete(labels = c("TSS 2kb", "TTS 2kb", "Exon", "Intron")) +
        labs(x = "", y = "Motif proportion ratio\n(6mA vs control)") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            # legend.position = "bottom", legend.direction = "horizontal",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri/motif_density_diff_denominator_core_genome_two_outgroup")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak_p <- get_plot(filename = "peak_motif_density")
nopeak_p <- get_plot(filename = "nopeak_motif_density")
ratio_p <- get_ratio_plot()
ylab <- peak_p$labels$y
nopeak_p$labels$y <- ""
peak_p$labels$x <- "6mA"
nopeak_p$labels$x <- "Control"

# peak_p | nopeak_p | ratio_p
ratio_p
ggsave("peak_nopeak_ratio_motif_base_count_proportion_in_genome_parts.pdf", width = 10, height = 10)
