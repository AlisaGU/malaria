#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(data.table)
library(ggridges)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        # data <- fread(paste0(var_type, "_mean_variant_count"), header = F, stringsAsFactors = F)
        data <- fread(paste0(var_type, "_mean_variant_count_peak"), header = F, stringsAsFactors = F, drop = c(7, 8))
        data <- as.data.frame(data)
        data <- data[, seq(1, ncol(data), 2)] / data[, seq(2, ncol(data), 2)]
        # data <- data.frame(value = unlist(data), chrom = rep(1:14, 4), source = rep(c("bp6", "bp9", "bp12", "bp15"), each = 14), var_type = var_type)
        data <- data.frame(value = unlist(data), chrom = rep(1:14, 3), source = rep(c("bp6", "bp9", "bp12"), each = 14), var_type = var_type)
    })
    result <- do.call(rbind, result)
    # result$source <- factor(result$source, levels = c("bp6", "bp9", "bp12", "bp15"))
    result$source <- factor(result$source, levels = c("bp6", "bp9", "bp12"))

    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))
    return(result)
}

read_peak_motif_scatter_each_chrom_peak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_all_possible_peak"), header = F, stringsAsFactors = F)
        data <- cbind(data[, 1], data[, 2], data[, 3] / data[, 4], var_type = var_type)
        colnames(data) <- c("chrom", "motif_length", "mutation_density", "var_type")
        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_chrom_nopeak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_all_possible_nopeak"), header = F, stringsAsFactors = F)
        data <- cbind(data[, 1], data[, 2], data[, 3] / data[, 4], var_type = var_type)
        colnames(data) <- c("chrom", "motif_length", "mutation_density", "var_type")
        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_len_peak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_all_possible_peak"), header = F, stringsAsFactors = F)
        data <- as.data.frame(data)
        colnames(data) <- c("chrom", "motif_length", "mutation_count", "region_len")

        data <- data[, -c(1)] %>%
            group_by(motif_length) %>%
            summarise(
                mutation_count = sum(mutation_count),
                region_len = sum(region_len)
            )
        data <- cbind(data[, 1], data[, 2] / data[, 3], var_type = var_type)
        colnames(data)[2] <- "mutation_density"
        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_len_nopeak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_all_possible_nopeak"), header = F, stringsAsFactors = F)
        data <- as.data.frame(data)
        colnames(data) <- c("chrom", "motif_length", "mutation_count", "region_len")

        data <- data[, -c(1)] %>%
            group_by(motif_length) %>%
            summarise(
                mutation_count = sum(mutation_count),
                region_len = sum(region_len)
            )
        data <- cbind(data[, 1], data[, 2] / data[, 3], var_type = var_type)
        colnames(data)[2] <- "mutation_density"
        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_len_total <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_all_possible_total"), header = F, stringsAsFactors = F)
        data <- as.data.frame(data)
        colnames(data) <- c("chrom", "motif_length", "mutation_count", "region_len")

        data <- data[, -c(1)] %>%
            group_by(motif_length) %>%
            summarise(
                mutation_count = sum(mutation_count),
                region_len = sum(region_len)
            )
        data <- cbind(data[, 1], data[, 2] / data[, 3], var_type = var_type)
        colnames(data)[2] <- "mutation_density"
        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_motif_peak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_each_motif_peak"), header = F, stringsAsFactors = F)
        data <- as.data.frame(data)
        colnames(data) <- c("chrom", "mutation_count", "motif_length")

        data <- data.frame(motif_length = data[, 3], mutation_density = data[, 2] / data[, 3], var_type = var_type)

        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

read_peak_motif_scatter_each_motif_nopeak <- function() {
    result <- lapply(c("snps", "indels"), function(var_type) {
        data <- fread(paste0(var_type, "_mean_variant_count_each_motif_nopeak"), header = F, stringsAsFactors = F)
        data <- as.data.frame(data)
        colnames(data) <- c("chrom", "mutation_count", "motif_length")

        data <- data.frame(motif_length = data[, 3], mutation_density = data[, 2] / data[, 3], var_type = var_type)

        return(data)
    })
    result <- do.call(rbind, result)
    result$var_type <- gsub("snps", "SNPs", gsub("indels", "Indels", result$var_type))
    result$var_type <- factor(result$var_type, levels = c("SNPs", "Indels"))

    return(result)
}

plot_boxplot <- function(data_for_plot = NULL) {
    p <- ggplot(data = data_for_plot, aes(x = source, y = value)) +
        # scale_y_continuous(expand = expansion(mult = 0.12)) +
        geom_boxplot(aes(color = source, group = source), outlier.shape = NA, show.legend = FALSE, lwd = 1.5) +
        geom_point(size = 8, aes(color = source)) +
        facet_wrap(vars(var_type), scales = "free_y") +
        stat_compare_means(
            # comparisons = list(c("bp15", "bp12"), c("bp15", "bp9"), c("bp15", "bp6")), aes(label = paste0("p = ", ..p.format..)), size = 12,
            comparisons = list(c("bp12", "bp9"), c("bp12", "bp6")), aes(label = paste0("p = ", ..p.format..)), size = 12,
            paired = T
        ) +
        scale_x_discrete(labels = c("6", "9", "12")) +
        scale_color_manual(values = c(
            "bp6" = rgb(128, 206, 242, maxColorValue = 255),
            "bp9" = rgb(110, 161, 214, maxColorValue = 255),
            "bp12" = rgb(82, 112, 180, maxColorValue = 255),
            "bp15" = rgb(65, 89, 165, maxColorValue = 255)
        )) +
        theme_bw() +
        labs(y = "Mutation density", x = "Merged motifs length (bp)") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
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


plot_scatter <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = motif_length, y = mutation_density)) +
        geom_point(color = "#8d8d8d", fill = "#78addc", size = 8, alpha = 0.5, shape = 21) +
        labs(y = "Mutation density", x = "Merged motif length") +
        facet_wrap(~var_type, scales = "free_y") +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_scatter_each_chrom <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = motif_length, y = mutation_density)) +
        geom_point(color = "#8d8d8d", fill = "#78addc", size = 9, alpha = 0.5, shape = 21) +
        labs(y = "Mutation density", x = "Merged motif length (bp)") +
        facet_wrap(~var_type, scales = "free_y") +
        coord_cartesian(ylim = c(0, 0.7)) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 15) +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_merged_motif_length_distri_nomatter_chrom <- function() {
    peak <- get_vector_nomatter_chrom(dataname = "snps_mean_variant_count_all_possible_peak")
    nonpeak <- get_vector_nomatter_chrom(dataname = "snps_mean_variant_count_all_possible_nopeak")
    data <- rbind(
        cbind(snp_peak, "peak"), cbind(snp_nonpeak, "nonpeak")
    )
    data <- data.frame(
        motif_len = as.numeric(data[, 1]), peak_type = data[, 2]
    )
    p <- ggplot(data = data, aes(x = peak_type, y = motif_len)) +
        geom_violin(fill = rgb(0, 72, 144, maxColorValue = 255)) +
        labs(x = "", y = "Merged motifs length (bp)") +
        scale_y_break(c(30, 85)) +
        scale_x_discrete(labels = c("Control", "Peak")) +
        coord_flip() +
        theme_bw() +
        theme(
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

get_vector_nomatter_chrom <- function(dataname = NULL) {
    snp_peak <- read.table(dataname, header = F, as.is = T) %>%
        group_by(V1, V2) %>%
        summarise(V4 = sum(V4))

    snp_peak_motif_count <- data.frame(motif_len = snp_peak[, 1], motif_count = snp_peak[, 2] / snp_peak[, 1])
    snp_peak_vector <- lapply(1:nrow(snp_peak_motif_count), function(i) {
        rep(snp_peak_motif_count[i, 1], snp_peak_motif_count[i, 2])
    })
    result <- do.call(c, snp_peak_vector)
    return(result)
}

plot_merged_motif_average_length_distri_each_chrom <- function() {
    peak <- get_vector_average_each_chrom(dataname = "snps_mean_variant_count_all_possible_peak")
    nonpeak <- get_vector_average_each_chrom(dataname = "snps_mean_variant_count_all_possible_nopeak")
    data <- rbind(cbind(peak, type = "Peak"), cbind(nonpeak, type = "Control"))
    p <- ggplot(data, aes(x = type, y = value, color = type)) +
        geom_point(position = position_jitterdodge(), size = 8, alpha = 0.5) +
        geom_boxplot(size = 3, show.legend = FALSE, outlier.shape = NA, alpha = 0.1, width = 0.5) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        scale_color_manual(values = c("Peak" = "#db161b", "Control" = "#1c4e6c")) +
        theme_bw() +
        labs(y = "Merged motifs \naverage length (bp)") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            # legend.position = "bottom",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

get_vector_average_each_chrom <- function(dataname = NULL) {
    snp_peak <- read.table(dataname, header = F, as.is = T)
    snp_peak_split <- split(snp_peak, f = snp_peak$V1, drop = T)
    mean_motif_len <- sapply(snp_peak_split, function(x) {
        sum(x[, 4]) / sum(x[, 4] / x[, 2])
    })
    result <- data.frame(chrom = names(mean_motif_len), value = mean_motif_len)
    return(result)
}

plot_merged_motif_median_length_distri_each_chrom <- function() {
    peak <- get_vector_median_each_chrom(dataname = "snps_mean_variant_count_all_possible_peak")
    nonpeak <- get_vector_median_each_chrom(dataname = "snps_mean_variant_count_all_possible_nopeak")
    data <- rbind(cbind(peak, type = "Peak"), cbind(nonpeak, type = "Control"))
    p <- ggplot(data, aes(x = type, y = value, color = type)) +
        geom_point(position = position_jitterdodge(), size = 8, alpha = 0.5) +
        geom_boxplot(size = 3, show.legend = FALSE, outlier.shape = NA, alpha = 0.1, width = 0.5) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        scale_color_manual(values = c("Peak" = "#db161b", "Control" = "#1c4e6c")) +
        theme_bw() +
        labs(y = "Merged motifs \naverage length (bp)") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            # legend.position = "bottom",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

get_vector_median_each_chrom <- function(dataname = NULL) {
    snp_peak <- read.table(dataname, header = F, as.is = T)
    snp_peak_split <- split(snp_peak, f = snp_peak$V1, drop = T)
    median_motif_len <- sapply(snp_peak_split, function(x) {
        a <- apply(x, 1, function(y) {
            count <- as.numeric(y[4]) / as.numeric(y[2])
            rep(as.numeric(unlist(y[2])), count)
        })
        a <- do.call(c, a)
        sort_a <- sort(a)
        return(median(sort_a))
    })
    result <- data.frame(chrom = names(mean_motif_len), value = mean_motif_len)
    return(result)
}

plot_chromosome_proportion_for_specific_length_motif <- function(data_for_plot = NULL) {
    data_for_plot_snp <- data_for_plot %>% filter(var_type == "SNPs")
    data_for_plot_indel <- data_for_plot %>% filter(var_type == "Indels")

    a <- data.frame(table(data_for_plot_snp$motif_length, data_for_plot_snp$chrom))
    a <- plyr::ddply(a, .(Var1), transform, percent = Freq / sum(Freq) * 100)
    a$label <- paste0(sprintf("%.1f", a$percent), "%")
    snp <- a

    a <- data.frame(table(data_for_plot_indel$motif_length, data_for_plot_indel$chrom))
    a <- plyr::ddply(a, .(Var1), transform, percent = Freq / sum(Freq) * 100)
    a$label <- paste0(sprintf("%.1f", a$percent), "%")
    indel <- a

    data_for_plot <- rbind(data.frame(snp, var_type = "SNPs"), data.frame(indel, var_type = "Indels"))
    data_for_plot$var_type <- factor(data_for_plot$var_type, levels = c("SNPs", "Indels"))
    p <- ggplot(data_for_plot, aes(Var1, percent, fill = Var2)) +
        geom_bar(stat = "identity", position = position_stack()) +
        facet_wrap(vars(var_type), nrow = 2, ncol = 1, scales = "free_x") +
        scale_fill_manual(values = c(
            "Pf3D7_01_v3" = rgb(133, 66, 100, maxColorValue = 255), "Pf3D7_02_v3" = rgb(82, 54, 98, maxColorValue = 255),
            "Pf3D7_03_v3" = rgb(128, 176, 127, maxColorValue = 255), "Pf3D7_04_v3" = rgb(37, 42, 84, maxColorValue = 255),
            "Pf3D7_05_v3" = rgb(229, 228, 173, maxColorValue = 255), "Pf3D7_06_v3" = rgb(25, 104, 123, maxColorValue = 255),
            "Pf3D7_07_v3" = rgb(185, 225, 216, maxColorValue = 255), "Pf3D7_08_v3" = rgb(86, 140, 103, maxColorValue = 255),
            "Pf3D7_09_v3" = rgb(248, 200, 200, maxColorValue = 255), "Pf3D7_10_v3" = rgb(184, 84, 157, maxColorValue = 255),
            "Pf3D7_11_v3" = rgb(57, 96, 174, maxColorValue = 255), "Pf3D7_12_v3" = rgb(160, 37, 94, maxColorValue = 255),
            "Pf3D7_13_v3" = rgb(231, 88, 13, maxColorValue = 255), "Pf3D7_14_v3" = rgb(97, 170, 196, maxColorValue = 255)
        )) +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.7, vjust = 0.7
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/different_length_overlap_motif")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
## boxplot
data_for_plot <- read_data_peak()
plot_boxplot(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_3part.pdf", width = 18, height = 12)


snp <- fread("snps_mean_variant_count", stringsAsFactors = F, select = c(2, 4, 6, 8))
snp_count <- t(apply(snp, 1, function(x) {
    x / c(6, 9, 12, 15)
}))
rownames(snp_count) <- paste("chrom", 1:14, sep = "")
colnames(snp_count) <- c("6bp", "9bp", "12bp", "15bp")
## 合并motif，同一长度平均密度，每个点代表一个染色体，散点图
data_for_plot <- read_peak_motif_scatter_each_chrom_peak()
plot_scatter_each_chrom(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_chrom_peak.pdf", width = 16, height = 8)

plot_chromosome_proportion_for_specific_length_motif(data_for_plot = data_for_plot)
ggsave("chromosome_proportion_for_specific_length_motif_peak.pdf", width = 18, height = 16)


data_for_plot <- read_peak_motif_scatter_each_chrom_nopeak()
plot_scatter_each_chrom(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_chrom_nopeak.pdf", width = 16, height = 8)

ggplot(data_for_plot, aes(x = motif_length, y = chrom))
## 合并motif，同一长度平均密度，散点图

data_for_plot <- read_peak_motif_scatter_each_len_peak()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_len_peak.pdf", width = 18, height = 9)


data_for_plot <- read_peak_motif_scatter_each_len_nopeak()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_len_nopeak.pdf", width = 18, height = 9)

data_for_plot <- read_peak_motif_scatter_each_len_total()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_len_total.pdf", width = 18, height = 9)
## 合并motif，不算平均密度，散点图


data_for_plot <- read_peak_motif_scatter_each_motif_peak()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_motif_peak.pdf", width = 18, height = 9)


data_for_plot <- read_peak_motif_scatter_each_motif_nopeak()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_motif_nopeak.pdf", width = 18, height = 9)

data_for_plot <- read_peak_motif_scatter_each_len_total()
plot_scatter(data_for_plot = data_for_plot)
ggsave("different_length_overlap_motif_scatter_each_len_total.pdf", width = 18, height = 9)
## 统计信号区和非信号区的合并motif平均长度分布
plot_merged_motif_average_length_distri_each_chrom()
ggsave("merged_motifs_average_length.pdf", height = 7, width = 7)
