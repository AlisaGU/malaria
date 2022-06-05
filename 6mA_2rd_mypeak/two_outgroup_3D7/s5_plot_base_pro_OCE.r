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
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/base_proporation_two_group_consistent_whole_chromsome")
    data <- read.table(paste0(path, "/", dataname, "_mean_proporation"), stringsAsFactors = F, sep = " ")
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
output_path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/base_proporation_two_group_consistent_whole_chromsome")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
color_value <- c(
    "nopeak" = "#868686", "relative0.01" = "#005a32", "relative0.05" = "#238443",
    "relative0.1" = "#41ab5d", "relative0.2" = "#78c679", "relative0.5" = "#addd8e",
    "absolute10" = "#a63603", "absolute30" = "#e6550d", "absolute50" = "#fd8d3c",
    "absolute100" = "#fdbe85", "absolute200" = "#feedde"
)

data <- list()
for (variant_type in c("A","T","C","G")) {
    data[[popu_symbol]][[variant_type]] <- get_data_for_plot(popu = popu_symbol, dataname = variant_type)
}

A <- plot_mutation_density(
    data_for_plot = data$OCE$A,
    sig_y = 0.45, line_y = 0.46, text_y = 0.47, y_min = 0.35, y_max = 0.48,
    submit_type = "OCE", ytitle = "A", xlabel = FALSE
)
T <- plot_mutation_density(
    data_for_plot = data$OCE$T,
    sig_y = 0.43, line_y = 0.44 ,text_y = 0.45, y_min = 0.32, y_max = 0.47,
    submit_type = "OCE", ytitle = "T", xlabel = FALSE
)
C <- plot_mutation_density(
    data_for_plot = data$OCE$C,
    sig_y = 0.123, line_y = 0.126, text_y = 0.135, y_min = 0.08, y_max = 0.128,
    submit_type =FALSE, ytitle = "C", xlabel = TRUE
)
G <- plot_mutation_density(
    data_for_plot = data$OCE$G,
    sig_y = 0.128, line_y = 0.132, text_y = 0.135, y_min = 0.08, y_max = 0.135,
    submit_type = FALSE, ytitle = "G", xlabel = TRUE
)

(A/C) | (T/G)  +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 32))
ggsave(paste0(output_path, "/", "base_pro.pdf"), width = 16, height = 16)
ggsave(paste0(output_path, "/", "base_pro.jpg"), width = 16, height = 16)
