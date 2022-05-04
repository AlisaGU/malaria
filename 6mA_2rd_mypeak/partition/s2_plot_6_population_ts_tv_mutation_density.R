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
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic/", popu_symbol, "/mutation_density")
    data <- read.table(paste0(path, "/", dataname), stringsAsFactors = F)
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
    inf_index <- which(is.infinite(data_for_plot$value))
    data_for_plot$value[inf_index] <- NA

    p <- ggboxplot(data = data_for_plot, x = "source", y = "value", color = "source", add = "point", add.params = list(size = 3)) +
        # stat_compare_means(aes(group = source),
        #     ref.group = "nopeak", hide.ns = T,method="wilcox.test",
        #     label = "p.signif", paired = T, label.y = sig_y, size = 9
        # ) +
        sapply(c(2:11), function(index) {
            annotate(geom = "text", x = index, y = sig_y, label = get_sigif(index = index, data_for_plot = data_for_plot), size = 9)
        }) +
        geom_point(data = data_for_plot[data_for_plot$label == "1", ], aes(x = source, y = value), color = rgb(0, 149, 249, maxColorValue = 255), size = 3) +
        scale_color_manual(values = color_value) +
        scale_x_discrete(labels = c("control", "1%", "5%", "10%", "20%", "50%", "10bp", "30bp", "50bp", "100bp", "200bp")) +
        geom_segment(aes(x = 2.0, xend = 6, y = line_y, yend = line_y), col = "#41ab5d") +
        geom_segment(aes(x = 7.0, xend = 11, y = line_y, yend = line_y), col = "#fd8d3c") +
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

get_sigif <- function(index = NULL, data_for_plot = NULL) {
    no <- data_for_plot[data_for_plot$source == "nopeak", 1]

    source <- c("nopeak", "relative0.01", "relative0.05", "relative0.1", "relative0.2", "relative0.5", "absolute10", "absolute30", "absolute50", "absolute100", "absolute200")
    index_value <- data_for_plot[data_for_plot$source == source[index], 1]
    if (all(is.na(index_value))) {
        return("")
    } else {
        no_na_inf_index <- intersect(
            intersect(
                which(!is.na(no)),
                which(!is.infinite(no))
            ),
            intersect(
                which(!is.na(index_value)),
                which(!is.infinite(index_value))
            )
        )

        p_value <- wilcox.test(no[no_na_inf_index], index_value[no_na_inf_index], paired = T)$p.value
        if (p_value > 0.05) {
            return("")
        } else if (p_value <= 0.05 & p_value > 0.01) {
            return("*")
        } else if (p_value <= 0.01 & p_value > 0.001) {
            return("**")
        } else if (p_value <= 0.001) {
            return("***")
        }
    }
}
# 3. input ---------------------------------------------------------------- TODO:
output_path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic"
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
    for (variant_type in c("transition_mean_variant_count", "transversion_mean_variant_count", "ts_tv_ratio.txt")) {
        data[[popu_symbol]][[variant_type]] <- get_data_for_plot(popu = popu_symbol, dataname = variant_type)
    }
}

# 第一幅图
{
    world_ts <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$transition_mean_variant_count,
        sig_y = 0.034, line_y = 0.037, text_y = 0.04, y_min = -0.002, y_max = 0.045,
        submit_type = "World", ytitle = "Transition", xlabel = FALSE
    )
    world_tv <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$transversion_mean_variant_count,
        sig_y = 0.04, line_y = 0.043, text_y = 0.045, y_min = -0.002, y_max = 0.05,
        submit_type = FALSE, ytitle = "Transversion", xlabel = FALSE
    )
    world_ratio <- plot_mutation_density(
        data_for_plot = data$SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE$ts_tv_ratio.txt,
        sig_y = 1.55, line_y = 1.62, text_y = 1.6, y_min = 0.5, y_max = 1.7,
        submit_type = FALSE, ytitle = "Ts/tv", xlabel = TRUE
    )

    Africa_ts <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$transition_mean_variant_count,
        sig_y = 0.033, line_y = 0.037, text_y = 0.04, y_min = -0.002, y_max = 0.045,
        submit_type = "Africa", ytitle = FALSE, xlabel = FALSE
    )
    Africa_tv <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$transversion_mean_variant_count,
        sig_y = 0.036, line_y = 0.04, text_y = 0.043, y_min = -0.002, y_max = 0.045,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    Africa_ratio <- plot_mutation_density(
        data_for_plot = data$WAF_CAF_EAF$ts_tv_ratio.txt,
        sig_y = 1.53, line_y = 1.6, text_y = 1.65, y_min = 0.5, y_max = 1.7,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    wsea_ts <- plot_mutation_density(
        data_for_plot = data$WSEA$transition_mean_variant_count,
        sig_y = 0.0062, line_y = 0.0068, text_y = 0.0071, y_min = -0.0001, y_max = 0.0075,
        submit_type = "Western SE Asia", ytitle = FALSE, xlabel = FALSE
    )
    wsea_tv <- plot_mutation_density(
        data_for_plot = data$WSEA$transversion_mean_variant_count,
        sig_y = 0.0065, line_y = 0.0071, text_y = 0.0071, y_min = -0.0001, y_max = 0.0075,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    wsea_ratio <- plot_mutation_density(
        data_for_plot = data$WSEA$ts_tv_ratio.txt,
        sig_y = 5.1, line_y = 5.3, text_y = 5.4, y_min = -0.1, y_max = 5.5,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    (world_ts | Africa_ts | wsea_ts) / (world_tv | Africa_tv | wsea_tv) / (world_ratio | Africa_ratio | wsea_ratio) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 32))
    ggsave(paste0(output_path, "/", "world_africa_wsea.ts_tv.mutation_density.pdf"), width = 26, height = 26)
}

# 第二幅图
{
    WAF_ts <- plot_mutation_density(
        data_for_plot = data$WAF$transition_mean_variant_count,
        sig_y = 0.025, line_y = 0.027, text_y = 0.029, y_min = -0.002, y_max = 0.031,
        submit_type = "West africa", ytitle = "Transition", xlabel = FALSE
    )
    WAF_tv <- plot_mutation_density(
        data_for_plot = data$WAF$transversion_mean_variant_count,
        sig_y = 0.028, line_y = 0.03, text_y = 0.032, y_min = -0.002, y_max = 0.032,
        submit_type = FALSE, ytitle = "Transversion", xlabel = FALSE
    )
    WAF_ratio <- plot_mutation_density(
        data_for_plot = data$WAF$ts_tv_ratio.txt,
        sig_y = 1.7, line_y = 1.8, text_y = 1.9, y_min = 0.5, y_max = 1.8,
        submit_type = FALSE, ytitle = "Ts/tv", xlabel = TRUE
    )

    CAF_ts <- plot_mutation_density(
        data_for_plot = data$CAF$transition_mean_variant_count,
        sig_y = 0.0125, line_y = 0.013, text_y = 0.014, y_min = -0.002, y_max = 0.015,
        submit_type = "Central africa", ytitle = FALSE, xlabel = FALSE
    )
    CAF_tv <- plot_mutation_density(
        data_for_plot = data$CAF$transversion_mean_variant_count,
        sig_y = 0.0074, line_y = 0.0076, text_y = 0.0078, y_min = -0.002, y_max = 0.008,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    CAF_ratio <- plot_mutation_density(
        data_for_plot = data$CAF$ts_tv_ratio.txt,
        sig_y = 5.1, line_y = 5.2, text_y = 5.3, y_min = 0, y_max = 5.5,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    EAF_ts <- plot_mutation_density(
        data_for_plot = data$EAF$transition_mean_variant_count,
        sig_y = 0.0147, line_y = 0.0151, text_y = 0.0157, y_min = -0.002, y_max = 0.017,
        submit_type = "East africa", ytitle = FALSE, xlabel = FALSE
    )
    EAF_tv <- plot_mutation_density(
        data_for_plot = data$EAF$transversion_mean_variant_count,
        sig_y = 0.012, line_y = 0.013, text_y = 0.014, y_min = -0.002, y_max = 0.015,
        submit_type = FALSE, ytitle = FALSE, xlabel = FALSE
    )
    EAF_ratio <- plot_mutation_density(
        data_for_plot = data$EAF$ts_tv_ratio.txt,
        sig_y = 1.77, line_y = 1.9, text_y = 1.9, y_min = 0, y_max = 2,
        submit_type = FALSE, ytitle = FALSE, xlabel = TRUE
    )

    (WAF_ts | CAF_ts | EAF_ts) / (WAF_tv | CAF_tv | EAF_tv) / (WAF_ratio | CAF_ratio | EAF_ratio) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 32))
    ggsave(paste0(output_path, "/", "WAF_CAF_EAF.ts_tv.mutation_density.pdf"), width = 26, height = 26)
}

CAF_ts <- plot_mutation_density(
    data_for_plot = data$CAF$transition_mean_variant_count,
    sig_y = 0.038, line_y = 0.04, text_y = 0.045, y_min = 0, y_max = 0.05,
    submit_type = "Central africa", ytitle = "Transition", xlabel = FALSE
)
CAF_tv <- plot_mutation_density(
    data_for_plot = data$CAF$transversion_mean_variant_count,
    sig_y = 0.074, line_y = 0.076, text_y = 0.0785, y_min = 0, y_max = 0.081,
    submit_type = FALSE, ytitle = "Transversion", xlabel = FALSE
)
CAF_ratio <- plot_mutation_density(
    data_for_plot = data$CAF$ts_tv_ratio.txt,
    sig_y = 1.35, line_y = 1.39, text_y = 1.39, y_min = 0.25, y_max = 1.42,
    submit_type = FALSE, ytitle = "Ts/Tv", xlabel = TRUE
)

CAF_ts / CAF_tv / CAF_ratio + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 32))
ggsave(paste0(output_path, "/", "CAF.ts_tv.mutation_density.pdf"), width = 10, height = 26)
ggsave(paste0(output_path, "/", "CAF.ts_tv.mutation_density.jpg"), width = 10, height = 26)

a <- data$CAF$transition_mean_variant_count
data_absolute <- a[grepl("absolute", a$source), ]
aov.abso <- aov(value ~ source, data = data_absolute)
summary(aov.abso)
data_relative <- a[grepl("relative", a$source), ]
aov.rela <- aov(value ~ source, data = data_relative)
summary(aov.rela)

a <- data$CAF$transversion_mean_variant_count
data_absolute <- a[grepl("absolute", a$source), ]
aov.abso <- aov(value ~ source, data = data_absolute)
summary(aov.abso)
data_relative <- a[grepl("relative", a$source), ]
aov.rela <- aov(value ~ source, data = data_relative)
summary(aov.rela)
