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
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent")
    data <- read.table(paste0(path, "/", dataname, "_mean_variant_count"), stringsAsFactors = F, sep = " ")
    colnames(data) <- c("nopeak", "motif", paste("region", c(1, 2, 3, 4), sep = ""))
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
            ref.group = "nopeak", hide.ns = T,
            label = "p.signif", paired = T, label.y = sig_y, size = 9
        ) +
        geom_point(data = data_for_plot[data_for_plot$label == "1", ], aes(x = source, y = value), color = rgb(0, 149, 249, maxColorValue = 255), size = 3) +
        scale_color_manual(values = color_value) +
        scale_x_discrete(labels = c("control", "motif", "1st", "2rd", "3th", "4th")) +
        geom_segment(aes(x = 3.0, xend = 6, y = line_y, yend = line_y), col = "#35978f") +
        # annotate("text", label = "Motif adjacent areas", x = 4.5, y = text_y, size = 24 / .pt) +
        coord_cartesian(ylim = c(y_min, y_max)) +
        labs(y = "") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
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
        p <- p + annotate("text", label = "Motif adjacent areas", x = 4.5, y = text_y, size = 24 / .pt)
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
popu_symbol <- "WSEA"
output_path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
color_value <- c(
    "nopeak" = "#868686",
    "motif" = "#ebbe58",
    "region1" = "#024943",
    "region2" = "#2e8f87",
    "region3" = "#77beb3",
    "region4" = "#bfe0db"
)

data <- list()
for (variant_type in c("A2G", "G2A", "T2C", "C2T", "A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")) {
    data[[popu_symbol]][[variant_type]] <- get_data_for_plot(popu = popu_symbol, dataname = variant_type)
}

A2G <- plot_mutation_density(
    data_for_plot = data$WSEA$A2G, sig_y = 0.0163, line_y = 0.0176, text_y = 0.0185, y_min = 0.0019, y_max = 0.019,
    submit_type = "transition", ytitle = "A2G", xlabel = FALSE
)
G2A <- plot_mutation_density(
    data_for_plot = data$WSEA$G2A,
    sig_y = 0.034, line_y = 0.036, text_y = 0.011, y_min = 0.006, y_max = 0.04,
    submit_type = FALSE, ytitle = "G2A", xlabel = FALSE
)
T2C <- plot_mutation_density(
    data_for_plot = data$WSEA$T2C,
    sig_y = 0.0162, line_y = 0.017, text_y = 0.0175, y_min = 0.0015, y_max = 0.018,
    submit_type = FALSE, ytitle = "T2C", xlabel = FALSE
)
C2T <- plot_mutation_density(
    data_for_plot = data$WSEA$C2T,
    sig_y = 0.0265, line_y = 0.0277, text_y = 0.011, y_min = 0.006, y_max = 0.028,
    submit_type = FALSE, ytitle = "C2T", xlabel = TRUE
)
# TODO:
A2C <- plot_mutation_density(
    data_for_plot = data$WSEA$A2C, sig_y = 0.0047, line_y = 0.0051, text_y = 0.0055, y_min = 0, y_max = 0.0055,
    submit_type = "transversion", ytitle = "A2C", xlabel = FALSE
)
A2T <- plot_mutation_density(
    data_for_plot = data$WSEA$A2T,
    sig_y = 0.014, line_y = 0.0148, text_y = 0.011, y_min = 0.0009, y_max = 0.015,
    submit_type = FALSE, ytitle = "A2T", xlabel = FALSE
)
G2C <- plot_mutation_density(
    data_for_plot = data$WSEA$G2C,
    sig_y = 0.0087, line_y = 0.0093, text_y = 0.0175, y_min = -0.00002, y_max = 0.0095,
    submit_type = FALSE, ytitle = "G2C", xlabel = FALSE
)
G2T <- plot_mutation_density(
    data_for_plot = data$WSEA$G2T,
    sig_y = 0.0148, line_y = 0.0157, text_y = 0.011, y_min = 0.0018, y_max = 0.016,
    submit_type = FALSE, ytitle = "G2T", xlabel = TRUE
)

C2A <- plot_mutation_density(
    data_for_plot = data$WSEA$C2A, sig_y = 0.0115, line_y = 0.0123, text_y = 0.013, y_min = 0.0008, y_max = 0.014,
    submit_type = "transversion", ytitle = "C2A", xlabel = FALSE
)
C2G <- plot_mutation_density(
    data_for_plot = data$WSEA$C2G,
    sig_y = 0.017, line_y = 0.018, text_y = 0.011, y_min = 0.0006, y_max = 0.019,
    submit_type = FALSE, ytitle = "C2G", xlabel = FALSE
)
T2A <- plot_mutation_density(
    data_for_plot = data$WSEA$T2A,
    sig_y = 0.0093, line_y = 0.0099, text_y = 0.042, y_min = 0.0006, y_max = 0.01,
    submit_type = FALSE, ytitle = "T2A", xlabel = FALSE
)
T2G <- plot_mutation_density(
    data_for_plot = data$WSEA$T2G,
    sig_y = 0.0052, line_y = 0.0056, text_y = 0.0078, y_min = 0.0002, y_max = 0.006,
    submit_type = FALSE, ytitle = "T2G", xlabel = TRUE
)

# fixation <- plot_mutation_density(
#     data_for_plot = data$WSEA$fixation,
#     sig_y = 0.03, line_y = 0.0315, text_y = 0.033, y_min = 0.015, y_max = 0.033,
#     submit_type = "Fixation", ytitle = "Fixation", xlabel = TRUE
# )

(A2G / G2A / T2C / C2T) | (A2C / A2T / G2C / G2T) | (C2A / C2G / T2A / T2G) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 32))
ggsave(paste0(output_path, "/", "WSEA.mutation_direc.pdf"), width = 26, height = 34)
ggsave(paste0(output_path, "/", "WSEA.mutation_direc.jpg"), width = 26, height = 34)



# fixation
# ggsave(paste0(output_path, "/", "WSEA.fixation.pdf"), width = 10, height = 10)
# ggsave(paste0(output_path, "/", "WSEA.fixation.jpg"), width = 10, height = 10)
