#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function() {
    mutation_density <- read.table(paste0(mutation_density_global_dir, "/snps_mean_variant_count_flank"), as.is = T, stringsAsFactors = F)
    colnames(mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    mutation_density_CDS <- data.frame(
        gene_name = mutation_density$gene_name,
        mutation_density = mutation_density$mutation_count / mutation_density$cds_length
    )


    peak_length <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_flank"), as.is = T, stringsAsFactors = F)
    colnames(peak_length) <- c("gene_name", "peak_length", "control_length")
    peak_proportion <- data.frame(gene_name = peak_length$gene_name, peak_proportion = apply(peak_length[, 2:3], 1, function(x) {
        x[1] / sum(x)
    }))

    result <- merge(mutation_density_CDS, peak_proportion, by = "gene_name")
    return(result)
}

read_data_with_peak_control_mutation_density <- function(filter = NULL) {
    data <- read.table("snps_mean_variant_count_peak_nopeak_flank", header = F, as.is = T)
    if (filter) {
        data <- data[!is.na(data[, 2]), ]
    }
    data <- data.frame(
        rbind(
            cbind(data[, 1], data[, 2] / data[, 3], "peak"),
            cbind(data[, 1], data[, 4] / data[, 5], "control")
        )
    )
    colnames(data) <- c("gene_name", "mutation_density_value", "peak_type")
    data$mutation_density_value <- as.numeric(as.character(data$mutation_density_value))

    gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(data, gene_class, by = "gene_name")
    return(result)
}

read_data_with_fold_enrichment_and_peak_proportion <- function(filter = NULL) {
    data <- read.table("snps_mean_variant_count_peak_nopeak_flank", header = F, as.is = T)
    if (filter) {
        data <- data[!is.na(data[, 2]), ]
    }
    data <- data.frame(data[, 1], data[, 2] / data[, 3], data[, 4] / data[, 5], data[, 6])
    colnames(data) <- c("gene_name", "peak_mutation_density", "control_mutation_density", "peak_FE")
    data_peak_proportion <- read_peak_proportion_and_mutation_density()
    data <- merge(data, data_peak_proportion, by = "gene_name")

    data <- data.frame(rbind(
        cbind(gene_name = data[, 1], part_mutation_density = data[, 2], data[, 4:6], peak_type = "peak"),
        cbind(gene_name = data[, 1], part_mutation_density = data[, 3], data[, 4:6], peak_type = "nopeak")
    ))

    gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(data, gene_class, by = "gene_name")
    return(result)
}

plot_graph <- function(data_for_plot = NULL, title = NULL) {
    p <- ggplot(
        data_for_plot,
        aes(x = peak_type, y = part_mutation_density, color = peak_type)
    ) +
        geom_boxplot(width = 0.7, lwd = 2) +
        # geom_point(position = position_jitterdodge(), size = 3) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        ggtitle(title) +
        labs(y = "Mutation density") +
        scale_color_manual(values = c("nopeak" = rgb(137, 157, 159, maxColorValue = 255), "peak" = rgb(236, 174, 34, maxColorValue = 255))) +
        scale_x_discrete(labels = c("Control", "6mA")) +
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
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome")
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_filter <- read_data_with_fold_enrichment_and_peak_proportion(filter = TRUE)

w <- plot_graph(data_for_plot = data_filter %>% filter(peak_proportion > 0.1), title = "Whole genome")
l <- plot_graph(data_for_plot = data_filter %>% filter(gene_class == "core") %>% filter(peak_proportion > 0.1), title = "Core genome")
h <- plot_graph(data_for_plot = data_filter %>% filter(gene_class != "core") %>% filter(peak_proportion > 0.1), title = "Non-core genome")

l$labels$y <- ""
h$labels$y <- ""
l | h
ggsave("whole_genome_peak_nopeak_mutation_density.pdf", width = 16, height = 8)
