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

get_data_for_plot_3_control_global <- function(dataname = NULL) {
    data <- fread(dataname, stringsAsFactors = F, sep = " ", drop = c(1, 3))
    # colnames(data) <- c("nopeak", "nopeak_motif", "nopeak500_motif", "motif", paste("region", c(1, 2, 3, 4), sep = ""))
    colnames(data) <- c("nopeak_motif", "motif", paste("region", c(1, 2, 3, 4), sep = ""))

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

get_data_for_plot_3_control_N2N <- function(dataname = NULL) {
    data <- fread(dataname, stringsAsFactors = F, sep = " ", select = c(2, 4))
    colnames(data) <- c("nopeak_motif", "motif")

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

get_data_for_plot_3_control_N2N_threshold <- function(dataname = NULL) {
    data <- fread(dataname, stringsAsFactors = F, sep = " ", select = c(2, 4))
    data <- data[, 2] / data[, 1]
    colnames(data) <- c("peak/control")

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

plot_global <- function(data_for_all_direc = NULL, y_pos = NULL, nrow = NULL, ncol = NULL, label = NULL) {
    anno_data_for_all_direc <- compare_means(value ~ source, group.by = "variant", data = data_for_all_direc) %>% mutate(y_pos = y_pos)

    p <- ggplot(data_for_all_direc, aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge()) +
        geom_point(position = position_jitterdodge(), size = 4)
    if (label == T) {
        p <- p + geom_point(data = data_for_all_direc[data_for_all_direc$label == "1", ], aes(x = source, y = value), color = rgb(0, 149, 249, maxColorValue = 255), size = 4)
    }
    p <- p + facet_wrap(~variant_f, nrow = nrow, ncol = ncol, scales = "free_y") + stat_compare_means(aes(group = source),
        # ref.group = "nopeak_motif",
        comparisons = list(c("nopeak_motif", "motif")),
        hide.ns = T, tip.length = 0.01, vjust = 0.5,
        label = "p.signif", paired = T, size = 10
    ) +
        scale_color_manual(values = color_value) +
        # scale_x_discrete(labels = c("Control", "C(motif)", "C500(motif)", "P_motif", paste("P_region", c(1, 2, 3, 4), sep = ""))) +
        scale_x_discrete(labels = c("Control(motif)", "P_motif", paste("P_region", c(1, 2, 3, 4), sep = ""))) +

        theme_bw() +
        theme(
            strip.background = element_rect(fill = "#ebbe58"),
            strip.text.x = element_text(size = 30),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 24, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_N2N <- function(data_for_all_direc = NULL, label = NULL) {
    p <- ggplot(data_for_all_direc, aes(x = variant_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4)
    if (label == T) {
        p <- p + geom_point(data = data_for_all_direc[data_for_all_direc$label == "1", ], aes(x = variant_f, y = value), color = rgb(0, 149, 249, maxColorValue = 255), size = 4)
    }
    p <- p + scale_color_manual(values = color_value_N2N, labels = c("Control motif", "Peak motif")) +
        stat_compare_means(aes(group = source),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 16
        ) +
        facet_wrap(~oribase_f, nrow = 2, ncol = 2, scale = "free_x") +
        labs(y = "Mutation density") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 24, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_N2N_threshold <- function(data_for_all_direc = NULL, label = NULL) {
    p <- ggplot(data_for_all_direc, aes(x = variant_f, y = value)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 4) +
        facet_wrap(~oribase_f, nrow = 2, ncol = 2, scale = "free_x") +
        labs(y = "Peak/Control") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 24, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
path <- Args[6]
# 4. variable setting of test module--------------------------------------- TODO:
# path<-"/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation_jiang/WSEA/motif_2/mutation_dentisy_two_outgroup_consistent"
# 5. process -------------------------------------------------------------- TODO:
setwd(path)

color_value <- c(
    "nopeak" = "black",
    "nopeak_motif" = "black",
    "nopeak500_motif" = "#2f2f2f",
    "motif" = "#b2182b",
    "region1" = "#2166ac",
    "region2" = "#4393c3",
    "region3" = "#92c5de",
    "region4" = "#d1e5f0"
)

{
    variants_type <- c("allvar", "snps", "indels")
    variants_type_title <- c("All variants", "SNPs", "Indels")
    data_for_plot <- lapply(variants_type, function(x) {
        get_data_for_plot_3_control_global(paste0(x, "_mean_variant_count_3control"))
    })
    names(data_for_plot) <- variants_type


    data_for_all_direc <- do.call(rbind, data_for_plot)
    data_for_all_direc$variant <- rep(variants_type_title, each = nrow(data_for_plot[[1]]))
    data_for_all_direc$variant_f <- factor(data_for_all_direc$variant, levels = variants_type_title)
    # data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak", "nopeak_motif", "nopeak500_motif", "motif", paste("region", c(1, 2, 3, 4), sep = "")))
    data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak_motif", "motif", paste("region", c(1, 2, 3, 4), sep = "")))

    plot_global(data_for_all_direc = data_for_all_direc, y_pos = 0.2, nrow = 1, ncol = 3, label = FALSE)
    # ggsave("N2N_3control.pdf", width = 32, height = 30)
    # ggsave("N2N_3control.jpg", width = 32, height = 30)
    ggsave("global_variants_3control_single_del.pdf", width = 24, height = 8)
    ggsave("global_variants_3control_single_del.jpg", width = 24, height = 8)
    plot_global(data_for_all_direc = data_for_all_direc, y_pos = 0.2, nrow = 4, ncol = 4, label = TRUE)
    # ggsave("N2N_3control_chr1.pdf", width = 32, height = 30)
    # ggsave("N2N_3control_chr1.jpg", width = 32, height = 30)
    ggsave("global_variants_3control_chr1_single_del.pdf", width = 24, height = 8)
    ggsave("global_variants_3control_chr1_single_del.jpg", width = 24, height = 8)
}

color_value_N2N <- c(
    "nopeak_motif" = "black",
    "motif" = "#b2182b"
)

{
    variants_type <- c(
        "A2G", "A2T", "A2C", "A_single_del",
        "G2A", "G2T", "G2C", "G_single_del",
        "C2T", "C2A", "C2G", "C_single_del",
        "T2C", "T2A", "T2G", "T_single_del"
    )
    variant_oribase <- substr(variants_type, 1, 1)
    # variants_type_title <- sub("2|_", " to ", variants_type, perl = T)
    variants_type_title <- c(
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    )
    data_for_plot <- lapply(variants_type, function(x) {
        # get_data_for_plot_3_control(paste0(x, "_mean_variant_count_3control"))
        get_data_for_plot_3_control_N2N(paste0(x, "_mean_variant_count_3control"))
    })
    names(data_for_plot) <- variants_type


    data_for_all_direc <- do.call(rbind, data_for_plot)
    data_for_all_direc$variant <- rep(variants_type_title, each = nrow(data_for_plot[[1]]))
    data_for_all_direc$variant_f <- factor(data_for_all_direc$variant, levels = variants_type_title)
    data_for_all_direc$oribase <- rep(variant_oribase, each = nrow(data_for_plot[[1]]))
    data_for_all_direc$oribase_f <- factor(data_for_all_direc$oribase, levels = unique(variant_oribase))

    # data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak", "nopeak_motif", "nopeak500_motif", "motif", paste("region", c(1, 2, 3, 4), sep = "")))
    data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak_motif", "motif"))

    plot_N2N(data_for_all_direc = data_for_all_direc, label = FALSE)
    # ggsave("N2N_3control.pdf", width = 32, height = 30)
    # ggsave("N2N_3control.jpg", width = 32, height = 30)
    ggsave("N2N_3control_single_del_2_2.pdf", width = 16, height = 16)
    # ggsave("N2N_3control_single_del_2_2.jpg", width = 16, height = 16)
    plot_N2N(data_for_all_direc = data_for_all_direc, label = TRUE)
    # ggsave("N2N_3control_chr1.pdf", width = 32, height = 30)
    # ggsave("N2N_3control_chr1.jpg", width = 32, height = 30)
    ggsave("N2N_3control_single_del_2_2_chr1.pdf", width = 16, height = 16)
    # ggsave("N2N_3control_single_del_2_2_chr1.jpg", width = 16, height = 16)
}

{{ variants_type <- c(
    "A2G", "A2T", "A2C", "A_single_del",
    "G2A", "G2T", "G2C", "G_single_del",
    "C2T", "C2A", "C2G", "C_single_del",
    "T2C", "T2A", "T2G", "T_single_del"
)
variant_oribase <- substr(variants_type, 1, 1)
# variants_type_title <- sub("2|_", " to ", variants_type, perl = T)
variants_type_title <- c(
    "A to G", "A to T", "A to C", "A to sDel",
    "G to A", "G to T", "G to C", "G to sDel",
    "C to T", "C to A", "C to G", "C to sDel",
    "T to C", "T to A", "T to G", "T to sDel"
)
data_for_plot <- lapply(variants_type, function(x) {
    # get_data_for_plot_3_control(paste0(x, "_mean_variant_count_3control"))
    get_data_for_plot_3_control_N2N_threshold(paste0(x, "_mean_variant_count_3control"))
})
names(data_for_plot) <- variants_type


data_for_all_direc <- do.call(rbind, data_for_plot)
data_for_all_direc$variant <- rep(variants_type_title, each = nrow(data_for_plot[[1]]))
data_for_all_direc$variant_f <- factor(data_for_all_direc$variant, levels = variants_type_title)
data_for_all_direc$oribase <- rep(variant_oribase, each = nrow(data_for_plot[[1]]))
data_for_all_direc$oribase_f <- factor(data_for_all_direc$oribase, levels = unique(variant_oribase))

# data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak", "nopeak_motif", "nopeak500_motif", "motif", paste("region", c(1, 2, 3, 4), sep = "")))
# data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak_motif", "motif"))

plot_N2N_threshold(data_for_all_direc = data_for_all_direc, label = FALSE)
ggsave("N2N_3control_single_del_2_2_threshold.pdf", width = 16, height = 16) }}


{
    pdf("mutation_vs.pdf", width = 32, height = 20)
    vs_list <- list(c("A2T", "T2A"), c("A2G", "G2A"), c("A2C", "C2A"), c("T2G", "G2T"), c("T2C", "C2T"), c("C2G", "G2C"))
    for (vs in vs_list) {
        data_for_all_direc <- do.call(rbind, data_for_plot[vs])
        data_for_all_direc$variant <- rep(vs, each = 112)
        data_for_all_direc$source <- factor(data_for_all_direc$source, levels = c("nopeak", "nopeak_motif", "nopeak500_motif", "motif", paste("region", c(1, 2, 3, 4), sep = "")))
        p <- ggplot(data_for_all_direc, aes(x = variant, y = value, color = variant)) +
            geom_boxplot(position = position_dodge()) +
            geom_point(position = position_jitterdodge()) +
            stat_compare_means(
                hide.ns = T,
                label = "p.signif", paired = T, size = 10
            ) +
            facet_wrap(~source, nrow = 2, ncol = 4, scales = "free_y") +
            theme_bw() +
            theme(
                strip.background = element_rect(fill = "#ebbe58"),
                strip.text.x = element_text(size = 30),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
                axis.text.y = element_text(size = 30, color = "black"),
                plot.title = element_text(
                    colour = "black",
                    size = 24, vjust = 0.5, hjust = 0.5
                ),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "None",
                plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
                panel.spacing = unit(3, "lines")
            )
        print(p)
    }
    dev.off()
}
