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
    data <- read.table(paste0(path, "/", dataname), stringsAsFactors = F, sep = " ")
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
            ref.group = "nopeak",
            label = "p.signif", paired = T, label.y = sig_y, size = 9
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
    if (!is.logical(ytitle)) {
        p <- p + labs(y = ytitle)
    }
    if (xlabel) {
        p <- p + theme(axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7))
    }
    return(p)
}


# 3. input ---------------------------------------------------------------- TODO:
output_path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
color_value <- c(
    "nopeak" = "#868686", "relative0.01" = "#005a32", "relative0.05" = "#238443",
    "relative0.1" = "#41ab5d", "relative0.2" = "#78c679", "relative0.5" = "#addd8e",
    "absolute10" = "#a63603", "absolute30" = "#e6550d", "absolute50" = "#fd8d3c",
    "absolute100" = "#fdbe85", "absolute200" = "#feedde"
)
data <- list()
for (popu_symbol in c("SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE", "WAF_CAF_EAF", "WSEA", "WAF", "CAF", "EAF")) {
    for (variant_type in c("allvar_mean_variant_count", "snps_mean_variant_count", "indels_mean_variant_count")) {
        data[[popu_symbol]][[variant_type]] <- get_data_for_plot(popu = popu_symbol, dataname = variant_type)
    }
}

# 第一幅图
{
    world_all <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$allvar_mean_variant_count,
        sig_y = 0.663, line_y = 0.69, text_y = 0.7, y_min = 0.2, y_max = 0.731,
        submit_type = "World", ytitle = "All variants \nmutation density", xlabel = FALSE
    )
    world_snp <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$snps_mean_variant_count,
        sig_y = 0.475, line_y = 0.505, text_y = 0.52, y_min = 0.1, y_max = 0.54,
        submit_type = FALSE, ytitle = "SNPs \nmutation density", xlabel = FALSE
    )
    world_indel <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$indels_mean_variant_count,
        sig_y = 0.468, line_y = 0.497, text_y = 0.513, y_min = 0.1, y_max = 0.54,
        submit_type = FALSE, ytitle = "Indels \nmutation density", xlabel = TRUE
    )

    Africa_all <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$allvar_mean_variant_count,
        sig_y = 0.65, line_y = 0.68, text_y = 0.7, y_min = 0.2, y_max = 0.72,
        submit_type = "Africa", ytitle = FALSE, xlabel = FALSE
    )
    Africa_snp <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$snps_mean_variant_count,
        sig_y = 0.46, line_y = 0.49, text_y = 0.51, y_min = 0.1, y_max = 0.53,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    Africa_indel <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$indels_mean_variant_count,
        sig_y = 0.46, line_y = 0.49, text_y = 0.51, y_min = 0.1, y_max = 0.53,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    wsea_all <- plot_mutation_density(
        data_for_plot = data$WSEA$allvar_mean_variant_count,
        sig_y = 0.27, line_y = 0.29, text_y = 0.3, y_min = 0.05, y_max = 0.33,
        submit_type = "Western SE Asia", ytitle = FALSE, xlabel = FALSE
    )
    wsea_snp <- plot_mutation_density(
        data_for_plot = data$WSEA$snps_mean_variant_count,
        sig_y = 0.18, line_y = 0.197, text_y = 0.205, y_min = 0, y_max = 0.21,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    wsea_indel <- plot_mutation_density(
        data_for_plot = data$WSEA$indels_mean_variant_count,
        sig_y = 0.212, line_y = 0.23, text_y = 0.24, y_min = 0, y_max = 0.25,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    (world_all | Africa_all | wsea_all) / (world_snp | Africa_snp | wsea_snp) / (world_indel | Africa_indel | wsea_indel) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 32))
    ggsave(paste0(output_path, "/", "world_africa_wsea.mutation_density.pdf"), width = 26, height = 26)
}

# 第二幅图
{
    WAF_all <- plot_mutation_density(
        data_for_plot = data$WAF$allvar_mean_variant_count,
        sig_y = 0.575, line_y = 0.605, text_y = 0.625, y_min = 0.15, y_max = 0.65,
        submit_type = "West africa", ytitle = "All variants \nmutation density", xlabel = FALSE
    )
    WAF_snp <- plot_mutation_density(
        data_for_plot = data$WAF$snps_mean_variant_count,
        sig_y = 0.375, line_y = 0.405, text_y = 0.425, y_min = 0, y_max = 0.45,
        submit_type = FALSE, ytitle = "SNPs \nmutation density", xlabel = FALSE
    )
    WAF_indel <- plot_mutation_density(
        data_for_plot = data$WAF$indels_mean_variant_count,
        sig_y = 0.415, line_y = 0.445, text_y = 0.465, y_min = 0, y_max = 0.49,
        submit_type = FALSE, ytitle = "Indels \nmutation density", xlabel = TRUE
    )

    CAF_all <- plot_mutation_density(
        data_for_plot = data$CAF$allvar_mean_variant_count,
        sig_y = 0.37, line_y = 0.39, text_y = 0.403, y_min = 0.05, y_max = 0.415,
        submit_type = "Central africa", ytitle = FALSE, xlabel = FALSE
    )
    CAF_snp <- plot_mutation_density(
        data_for_plot = data$CAF$snps_mean_variant_count,
        sig_y = 0.255, line_y = 0.275, text_y = 0.288, y_min = 0, y_max = 0.305,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    CAF_indel <- plot_mutation_density(
        data_for_plot = data$CAF$indels_mean_variant_count,
        sig_y = 0.285, line_y = 0.305, text_y = 0.318, y_min = 0, y_max = 0.33,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    EAF_all <- plot_mutation_density(
        data_for_plot = data$EAF$allvar_mean_variant_count,
        sig_y = 0.44, line_y = 0.46, text_y = 0.473, y_min = 0.1, y_max = 0.485,
        submit_type = "East africa", ytitle = FALSE, xlabel = FALSE
    )
    EAF_snp <- plot_mutation_density(
        data_for_plot = data$EAF$snps_mean_variant_count,
        sig_y = 0.3, line_y = 0.32, text_y = 0.333, y_min = 0, y_max = 0.34,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    EAF_indel <- plot_mutation_density(
        data_for_plot = data$EAF$indels_mean_variant_count,
        sig_y = 0.33, line_y = 0.35, text_y = 0.363, y_min = 0, y_max = 0.38,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    (WAF_all | CAF_all | EAF_all) / (WAF_snp | CAF_snp | EAF_snp) / (WAF_indel | CAF_indel | EAF_indel) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 32))
    ggsave(paste0(output_path, "/", "WAF_CAF_EAF.mutation_density.pdf"), width = 26, height = 26)
}