#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_data_motif_nonmotif <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)

        data <- data.frame(
            rbind(
                cbind(data[, 1], data[, 2] / data[, 3], "motif"),
                cbind(data[, 1], data[, 4] / data[, 5], "nonmotif")
            ), gene_set
        )
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "mutation_density_value", "peak_type", "gene_set_type")
    result <- as.data.frame(result)
    result$mutation_density_value <- as.numeric(as.character(result$mutation_density_value))
    result$peak_type <- factor(result$peak_type, levels = c("nonmotif", "motif"))
    return(result)
}

read_data_Mergedmotif_nonmotif <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count_Merged"), header = F, as.is = T)

        data <- data.frame(
            rbind(
                cbind(data[, 1], data[, 2] / data[, 3], "motif"),
                cbind(data[, 1], data[, 4] / data[, 5], "nonmotif")
            ), gene_set
        )
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "mutation_density_value", "peak_type", "gene_set_type")
    result <- as.data.frame(result)
    result$mutation_density_value <- as.numeric(as.character(result$mutation_density_value))
    result$peak_type <- factor(result$peak_type, levels = c("nonmotif", "motif"))
    return(result)
}

read_data_peakmotif_nonpeakmotif_other <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count_peakmotif_nonpeakmotif_other"), header = F, as.is = T)

        data <- data.frame(
            rbind(
                cbind(data[, 1], data[, 2] / data[, 3], "peak_motif"),
                cbind(data[, 1], data[, 4] / data[, 5], "peak_nonmotif"),
                cbind(data[, 1], data[, 6] / data[, 7], "control_motif"),
                cbind(data[, 1], data[, 8] / data[, 9], "control_nonmotif")
            ), gene_set
        )
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "mutation_density_value", "peak_type", "gene_set_type")
    result <- as.data.frame(result)
    result$mutation_density_value <- as.numeric(as.character(result$mutation_density_value))
    result$peak_type <- factor(result$peak_type, levels = c("peak_motif", "peak_nonmotif", "control_motif", "control_nonmotif"))
    return(result)
}

read_base_content <- function() {
    data <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_exclude_pseudo/", gene_set, "/", gene_set, ".base_content"),
            sep = "\t", header = F, as.is = T, row.names = 1
        )
        colnames(data) <- c("A", "T", "C", "G")
        data <- data.frame(
            value = unlist(data), base = rep(c("A", "T", "C", "G"), each = nrow(data)),
            gene_set = gene_set, gene_name = rep(rownames(data), ncol(data))
        )
    })
    data <- do.call(rbind, data)
    data$base <- factor(data$base, levels = c("A", "T", "C", "G"))
    return(data)
}

plot_motif <- function(data_for_plot = NULL) {
    p <- ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value)) +
        scale_y_continuous(expand = expansion(mult = 0.12)) +
        facet_wrap(vars(gene_set_type), scales = "free_y") +
        geom_boxplot(aes(color = peak_type, group = peak_type), outlier.shape = NA, show.legend = FALSE, lwd = 1.5) +
        geom_point(size = 8, aes(color = peak_type)) +
        stat_compare_means(comparisons = list(c("nonmotif", "motif")), aes(label = paste0("p = ", ..p.format..)), size = 12, paired = T) +
        scale_x_discrete(labels = c("Non-motif", "Motif")) +
        scale_color_manual(values = c(
            "nonmotif" = "grey", "motif" = "red"
        )) +
        theme_bw() +
        labs(y = "Mutation density") +
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
            axis.title.x = element_blank(),
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

plot_base_content <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = base, y = value, fill = base)) +
        geom_boxplot() +
        geom_point(shape = 21, size = 6, position = position_jitterdodge()) +
        facet_wrap(~gene_set) +
        theme_bw() +
        labs(y = "% Base") +
        scale_fill_manual(values = c(
            "A" = "#bad7e8", "T" = "#3079b3",
            "C" = "#d5ed8d", "G" = "#32a720"
        )) +
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

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_6_sets")
gene_sets <- c("DR", "HDR", "RNA_translation", "STEVOR", "VAR", "RIF")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read_data_motif_nonmotif()
plot_motif(data_for_plot = data_for_plot)
ggsave("gene_6_sets.pdf", width = 18, height = 12)

data_for_plot <- read_data_Mergedmotif_nonmotif()
plot_motif(data_for_plot = data_for_plot)
ggsave("gene_6_sets_Merged_motifs.pdf", width = 18, height = 12)


data_for_plot <- read_data_peakmotif_nonpeakmotif_other()
ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value)) +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    facet_wrap(vars(gene_set_type), scales = "free_y") +
    geom_boxplot(aes(color = peak_type, group = peak_type), outlier.shape = NA, show.legend = FALSE, lwd = 1.5) +
    geom_point(size = 8, aes(color = peak_type)) +
    # stat_compare_means(comparisons = list(c("nonmotif", "motif")), aes(label = paste0("p = ", ..p.format..)), size = 12, paired = T) +
    scale_x_discrete(labels = c("Peak\nmotif", "Peak\nnonmotif", "Control\nmotif", "Control\nnonmotif")) +
    scale_color_manual(values = c(
        "control_nonmotif" = "grey", "control_motif" = "#575757", "peak_motif" = "red", "peak_nonmotif" = "#f06d6d"
    )) +
    theme_bw() +
    labs(y = "Mutation density") +
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
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        # legend.position = "bottom", legend.direction = "horizontal",
        legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
ggsave("gene_6_sets_peakmotif_nonpeakmotif_other.pdf", width = 27, height = 15)




## 探究为什么基因内的peak motif突变密度不高: motif区间长度
library(data.table)
result <- lapply(gene_sets, function(gene_set) {
    data <- fread(paste0(gene_set, "/snps_mean_variant_count_peakmotif_nonpeakmotif_other"), stringsAsFactors = F, select = 3)
    return(cbind(data, gene_set))
})
result <- do.call(rbind, result)
colnames(result) <- c("peak_motif_length", "gene_set")
ggplot(result, aes(x = peak_motif_length, color = gene_set)) +
    geom_density() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("peak_motif_region_len_of_each_gene_set.pdf", width = 10, height = 7)

## A,T,C,G占序列的比例
data_for_plot <- read_base_content()
plot_base_content(data_for_plot = data_for_plot)
ggsave("peak_motif_base_content_of_each_gene_set.pdf", width = 18, height = 14)
