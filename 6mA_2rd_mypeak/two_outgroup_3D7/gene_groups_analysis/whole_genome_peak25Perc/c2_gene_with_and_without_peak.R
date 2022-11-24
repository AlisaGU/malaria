#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(grid)
library(ggsignif)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function(gene_class_criterion = NULL, flank_length = NULL) {
    snp_mutation_density <- read.table(paste0(mutation_density_global_dir, "/snps_mean_variant_count_flank", flank_length), as.is = T, stringsAsFactors = F)
    colnames(snp_mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    indel_mutation_density <- read.table(paste0(mutation_density_global_dir, "/indels_mean_variant_count_flank", flank_length), as.is = T, stringsAsFactors = F)
    colnames(indel_mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    mutation_density <- rbind(data.frame(snp_mutation_density, type = "SNP"), data.frame(indel_mutation_density, type = "Indel"))

    mutation_density_CDS <- data.frame(
        gene_name = mutation_density$gene_name,
        mutation_density = mutation_density$mutation_count / mutation_density$cds_length,
        type = mutation_density$type
    )


    peak_length_25Perc <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_flank", flank_length, "_peak25Perc"), as.is = T, stringsAsFactors = F)
    colnames(peak_length_25Perc) <- c("gene_name", "peak_length_25Perc", "control_length_25Perc", "peak_count_25Perc")
    peak_proportion_25Perc <- data.frame(
        gene_name = peak_length_25Perc$gene_name,
        peak_proportion_25Perc = apply(peak_length_25Perc[, 2:3], 1, function(x) {
            x[1] / sum(x)
        }),
        peak_count_25Perc = peak_length_25Perc$peak_count_25Perc
    )

    peak_length_wholePeak <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_flank", flank_length, "_wholePeak"), as.is = T, stringsAsFactors = F)
    colnames(peak_length_wholePeak) <- c("gene_name", "peak_length_wholePeak", "control_length_wholePeak")
    peak_proportion_wholePeak <- data.frame(
        gene_name = peak_length_wholePeak$gene_name,
        peak_proportion_wholePeak = apply(peak_length_wholePeak[, 2:3], 1, function(x) {
            x[1] / sum(x)
        })
    )

    gene_class <- ""
    if (gene_class_criterion == "h3k9me3") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class_h3k9me3.bed"), as.is = T, stringsAsFactors = F)
    } else if (gene_class_criterion == "core") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    }
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(merge(merge(mutation_density_CDS, peak_proportion_25Perc, by = "gene_name"), peak_proportion_wholePeak, by = "gene_name"), gene_class, by = "gene_name")

    result$class <- sapply(1:nrow(result), function(x) {
        x <- result[x, ]
        if (x["peak_proportion_25Perc"] == 0 & x["peak_proportion_wholePeak"] == 0) {
            return("Control")
        } else if (x["peak_proportion_25Perc"] == 0 & x["peak_proportion_wholePeak"] != 0) {
            return("Mixed")
        } else if (x["peak_proportion_25Perc"] != 0) {
            return("Peak")
        }
    })
    return(result)
}

plot_boxplot_youwu6mA_3parts <- function(data_for_plot = NULL, title = NULL) {
    position1 <- max(data_for_plot$mutation_density)
    position2 <- position1 + 0.2 * (max(data_for_plot$mutation_density) - min(data_for_plot$mutation_density))
    position3 <- position1 + 0.1 * (max(data_for_plot$mutation_density) - min(data_for_plot$mutation_density))
    p <- ggplot(data_for_plot, aes(x = class, y = mutation_density)) +
        geom_boxplot(aes(color = class), outlier.shape = NA, size = 2) +
        geom_point(aes(color = class), size = 4) +
        geom_signif(
            comparisons = list(c("Control", "Mixed"), c("Control", "Peak"), c("Mixed", "Peak")),
            y_position = c(position1, position2, position3)
        ) +
        labs(y = "Mutation density") +
        theme_bw() +
        ggtitle(title) +
        scale_color_manual(values = c(
            "Control" = rgb(190, 229, 161, maxColorValue = 255),
            "Mixed" = rgb(84, 174, 171, maxColorValue = 255),
            "Peak" = rgb(233, 94, 70, maxColorValue = 255)
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
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )

    return(p)
}

plot_boxplot_youwu6mA_2parts <- function(data_for_plot = NULL, title = NULL) {
    position1 <- max(data_for_plot$mutation_density)
    p <- ggplot(data_for_plot, aes(x = class, y = mutation_density)) +
        geom_boxplot(aes(color = class), outlier.shape = NA, size = 2) +
        geom_point(aes(color = class), size = 4) +
        # geom_signif(
        #     # test = "wilcox.test",
        #     test = "t.test",
        #     comparisons = list(c("Control", "Peak")),
        #     y_position = position1, textsize = 12
        # ) +
        stat_signif(comparisons = list(c("Control", "Peak")), test = "t.test", test.args = list(alternative = "less")) +
        labs(y = "Mutation density") +
        theme_bw() +
        ggtitle(title) +
        scale_color_manual(values = c(
            "Control" = rgb(137, 157, 159, maxColorValue = 255),
            "Peak" = rgb(236, 174, 34, maxColorValue = 255)
        )) +
        scale_y_continuous(limits = c(0, position1 * 1.2)) +
        scale_x_discrete(labels = c("Control", "6mA")) +
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
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/wholeGenome_peak25Perc/3D7"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

for (gene_class_criterion in c("core")) {
    for (flank_length in c("0kb", "1kb", "2kb")) {
        data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
        core_snp <- plot_boxplot_youwu6mA_3parts(data_for_plot = data %>% filter(gene_class == "core") %>% filter(type == "SNP"), title = paste(c("Core", "SNP", paste0("flank", flank_length)), collapse = "_"))
        core_indel <- plot_boxplot_youwu6mA_3parts(data_for_plot = data %>% filter(gene_class == "core") %>% filter(type == "Indel"), title = paste(c("Core", "Indel", paste0("flank", flank_length)), collapse = "_"))
        noncore_snp <- plot_boxplot_youwu6mA_3parts(data_for_plot = data %>% filter(gene_class == "noncore") %>% filter(type == "SNP"), title = paste(c("Noncore", "SNP", paste0("flank", flank_length)), collapse = "_"))
        noncore_indel <- plot_boxplot_youwu6mA_3parts(data_for_plot = data %>% filter(gene_class == "noncore") %>% filter(type == "Indel"), title = paste(c("Noncore", "Indel", paste0("flank", flank_length)), collapse = "_"))
        (core_snp | core_indel) / (noncore_snp | noncore_indel)
        ggsave(paste0(mutation_density_global_dir, "/comparing_gene_with_and_without_peak_3parts_flank", flank_length, ".pdf"), width = 16, height = 16)
    }
}


for (gene_class_criterion in c("core")) {
    for (flank_length in c("0kb", "1kb", "2kb")) {
        data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
        core_snp <- plot_boxplot_youwu6mA_2parts(data_for_plot = data %>% filter(class !=
            "Mixed") %>% filter(gene_class == "core") %>% filter(type == "SNP"), title = paste(c("Core", "SNP", paste0("flank", flank_length)), collapse = "_"))
        core_indel <- plot_boxplot_youwu6mA_2parts(data_for_plot = data %>% filter(class !=
            "Mixed") %>% filter(gene_class == "core") %>% filter(type == "Indel"), title = paste(c("Core", "Indel", paste0("flank", flank_length)), collapse = "_"))
        noncore_snp <- plot_boxplot_youwu6mA_2parts(data_for_plot = data %>% filter(class !=
            "Mixed") %>% filter(gene_class == "noncore") %>% filter(type == "SNP"), title = paste(c("Noncore", "SNP", paste0("flank", flank_length)), collapse = "_"))
        noncore_indel <- plot_boxplot_youwu6mA_2parts(data_for_plot = data %>% filter(class !=
            "Mixed") %>% filter(gene_class == "noncore") %>% filter(type == "Indel"), title = paste(c("Noncore", "Indel", paste0("flank", flank_length)), collapse = "_"))
        (core_snp | core_indel) / (noncore_snp | noncore_indel)
        ggsave(paste0(mutation_density_global_dir, "/comparing_gene_with_and_without_peak_2parts_flank", flank_length, "_t.test.pdf"), width = 16, height = 16)
    }
}
