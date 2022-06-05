#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)
# 2. functions ------------------------------------------------------------ TODO:
get_data_for_plot <- function(popu_symbol = NULL, dataname = NULL) {
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/mutation_density")
    data <- read.table(paste0(path, "/", dataname, "_mean_variant_count"), stringsAsFactors = F, sep = " ")
    colnames(data) <- c(
        "nopeak", paste("relative", c(0.01, 0.05, 0.1, 0.2, 0.5), sep = ""),
        paste("absolute", c(10, 30, 50, 100, 200), sep = "")
    )
    rownames(data) <- paste("chrom", 1:14, sep = "")

    data_for_plot <- data.frame(
        value = unlist(data),
        chrom = rep(rownames(data), times = ncol(data)),
        source = rep(colnames(data), each = nrow(data))
    )
    a <- gsub("chrom", "", data_for_plot$chrom)
    a[a != "1"] <- ""
    data_for_plot$label <- a
    return(data_for_plot)
}

plot_mutation_density <- function(data_for_plot = NULL, sig_y = NULL, line_y = NULL, text_y = NULL, y_min = NULL, y_max = NULL, submit_type = NULL, ytitle = NULL, xlabel = NULL) {
    p <- ggboxplot(data = data_for_plot, x = "source", y = "value", color = "source", add = "point", add.params = list(size = 3)) +
        stat_compare_means(aes(group = source),
            ref.group = "nopeak", hide.ns = TRUE,
            label = "p.signif", paired = TRUE, label.y = sig_y, size = 9
        ) +
        geom_point(data = data_for_plot[data_for_plot$label == "1", ], aes(x = source, y = value), color = rgb(0, 149, 249, maxColorValue = 255), size = 3) +
        scale_color_manual(values = color_value) +
        scale_x_discrete(labels = c("control", "1%", "5%", "10%", "20%", "50%", "10bp", "30bp", "50bp", "100bp", "200bp")) +
        geom_segment(aes(x = 2.0, xend = 6, y = line_y, yend = line_y), col = "#41ab5d") +
        geom_segment(aes(x = 7.0, xend = 11, y = line_y, yend = line_y), col = "#fd8d3c") +
        # annotate("text", label = "Proportion summit", x = 4, y = text_y, size = 24 / .pt) +
        # annotate("text", label = "Absolute summit", x = 9, y = text_y, size = 24 / .pt) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(y_min, y_max)) +
        labs(y = "") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
            axis.text.y = element_text(size = 24, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 28, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 28, color = "black"),
            legend.position = "None",
            plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
        )
    if (!is.logical(submit_type)) {
        p <- p + annotate("text", label = "Proportion summit", x = 4, y = text_y, size = 24 / .pt) +
            annotate("text", label = "Absolute summit", x = 9, y = text_y, size = 24 / .pt)
    }
    if (!is.logical(submit_type)) {
        p <- p + labs(title = submit_type)
    }
    if (!is.logical(ytitle)) {
        p <- p + labs(y = ytitle)
    }
    if (xlabel) {
        p <- p + theme(axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7))
    }
    return(p)
}


# 3. input ---------------------------------------------------------------- TODO:
popu_symbol <- "OCE"
output_path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/mutation_density")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
color_value <- c(
    "nopeak" = "#868686", "relative0.01" = "#005a32", "relative0.05" = "#238443",
    "relative0.1" = "#41ab5d", "relative0.2" = "#78c679", "relative0.5" = "#addd8e",
    "absolute10" = "#a63603", "absolute30" = "#e6550d", "absolute50" = "#fd8d3c",
    "absolute100" = "#fdbe85", "absolute200" = "#feedde"
)

data <- list()
for (variant_type in c("allvar", "snps", "indels", "A2G", "G2A", "T2C", "C2T", "A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")) {
    data[[popu_symbol]][[variant_type]] <- get_data_for_plot(popu = popu_symbol, dataname = variant_type)
}

allvar <- plot_mutation_density(
    data_for_plot = data$OCE$allvar,
    sig_y = 0.123, line_y = 0.13, text_y = 0.135, y_min = 0.02, y_max = 0.14,
    submit_type = "OCE", ytitle = "All variants", xlabel = FALSE
)
snps <- plot_mutation_density(
    data_for_plot = data$OCE$snps,
    sig_y = 0.0285, line_y = 0.03, text_y = 0.0125, y_min = 0.004, y_max = 0.031,
    submit_type = FALSE, ytitle = "SNPs", xlabel = FALSE
)
indels <- plot_mutation_density(
    data_for_plot = data$OCE$indels,
    sig_y = 0.096, line_y = 0.102, text_y = 0.0105, y_min = 0.01, y_max = 0.11,
    submit_type = FALSE, ytitle = "Indels", xlabel = TRUE
)

(allvar / snps / indels) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 32))
ggsave(paste0(output_path, "/", "OCE.global.pdf"), width = 10, height = 26)
ggsave(paste0(output_path, "/", "OCE.global.jpg"), width = 10, height = 26)

A2G <- plot_mutation_density(
    data_for_plot = data$OCE$A2G,
    sig_y = 0.004, line_y = 0.0044, text_y = 0.0048, y_min = 0, y_max = 0.005,
    submit_type = "transition", ytitle = "A2G", xlabel = FALSE
)
G2A <- plot_mutation_density(
    data_for_plot = data$OCE$G2A,
    sig_y = 0.004, line_y = 0.0044, text_y = 0.0125, y_min = 0, y_max = 0.0046,
    submit_type = FALSE, ytitle = "G2A", xlabel = FALSE
)
T2C <- plot_mutation_density(
    data_for_plot = data$OCE$T2C,
    sig_y = 0.0029, line_y = 0.0031, text_y = 0.0105, y_min = 0, y_max = 0.0032,
    submit_type = FALSE, ytitle = "T2C", xlabel = FALSE
)
C2T <- plot_mutation_density(
    data_for_plot = data$OCE$C2T,
    sig_y = 0.0044, line_y = 0.0047, text_y = 0.0105, y_min = 0, y_max = 0.005,
    submit_type = FALSE, ytitle = "C2T", xlabel = TRUE
)
# TODO:
A2C <- plot_mutation_density(
    data_for_plot = data$OCE$A2C, sig_y = 0.00143, line_y = 0.0016, text_y = 0.0017, y_min = -0.0003, y_max = 0.0018,
    submit_type = "transversion", ytitle = "A2C", xlabel = FALSE
)
A2T <- plot_mutation_density(
    data_for_plot = data$OCE$A2T,
    sig_y = 0.0037, line_y = 0.0039, text_y = 0.011, y_min = 0, y_max = 0.0041,
    submit_type = FALSE, ytitle = "A2T", xlabel = FALSE
)
G2C <- plot_mutation_density(
    data_for_plot = data$OCE$G2C,
    sig_y = 0.00058, line_y = 0.00063, text_y = 0.0175, y_min = -0.00002, y_max = 0.0007,
    submit_type = FALSE, ytitle = "G2C", xlabel = FALSE
)
G2T <- plot_mutation_density(
    data_for_plot = data$OCE$G2T,
    sig_y = 0.002, line_y = 0.0022, text_y = 0.001, y_min = 0, y_max = 0.0023,
    submit_type = FALSE, ytitle = "G2T", xlabel = TRUE
)

C2A <- plot_mutation_density(
    data_for_plot = data$OCE$C2A, sig_y = 0.00165, line_y = 0.00177, text_y = 0.0019, y_min = 0, y_max = 0.002,
    submit_type = "transversion", ytitle = "C2A", xlabel = FALSE
)
C2G <- plot_mutation_density(
    data_for_plot = data$OCE$C2G,
    sig_y = 0.00099, line_y = 0.00105, text_y = 0.011, y_min = -0.00002, y_max = 0.0011,
    submit_type = FALSE, ytitle = "C2G", xlabel = FALSE
)
T2A <- plot_mutation_density(
    data_for_plot = data$OCE$T2A,
    sig_y = 0.0044, line_y = 0.0047, text_y = 0.003, y_min = 0, y_max = 0.005,
    submit_type = FALSE, ytitle = "T2A", xlabel = FALSE
)
T2G <- plot_mutation_density(
    data_for_plot = data$OCE$T2G,
    sig_y = 0.0024, line_y = 0.0026, text_y = 0.0078, y_min = 0, y_max = 0.0026,
    submit_type = FALSE, ytitle = "T2G", xlabel = TRUE
)

(A2G / G2A / T2C / C2T) | (A2C / A2T / G2C / G2T) | (C2A / C2G / T2A / T2G) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 32))
ggsave(paste0(output_path, "/", "OCE.mutation_direc.pdf"), width = 26, height = 34)
ggsave(paste0(output_path, "/", "OCE.mutation_direc.jpg"), width = 26, height = 34)
