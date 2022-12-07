#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function(gene_class_criterion = NULL) {
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

    gene_class <- ""
    if (gene_class_criterion == "h3k9me3") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class_h3k9me3.bed"), as.is = T, stringsAsFactors = F)
    } else if (gene_class_criterion == "core") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    }
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(merge(mutation_density_CDS, peak_proportion, by = "gene_name"), gene_class, by = "gene_name")
    return(result)
}

read_6mA_site_density <- function() {
    site_info <- read.table(paste0(peak_proportion_global_dir, "/gene_flank_6mA_site_info"), as.is = T, stringsAsFactors = F)
    data <- split(site_info, f = site_info$V4)

    result <- lapply(data, function(data_subset) {
        leng <- unique(data_subset$V3 - data_subset$V2)
        site_count <- NROW(data_subset)
        return(site_count / leng)
    })

    result <- do.call(c, result)
    result <- data.frame(gene_name = names(result), site_count_density = result)
    return(result)
}

plot_scatter <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density)) +
        geom_point(size = 4, color = rgb(0, 0, 139, maxColorValue = 255), alpha = 0.3) +
        # geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, color = "red") +
        labs(x = "Proportion of peak\nin gene and flank 2kb", y = "Mutation density") +
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

plot_scatter1 <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = site_count_density, y = mutation_density)) +
        geom_point(size = 2) +
        # geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12, color = "red") +
        labs(x = "6mA count density\nin gene and flank 2kb", y = "Mutation density") +
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

plot_boxplot_youwu6mA <- function(data_for_plot = NULL, title = NULL) {
    data_for_plot$has_peak <- data_for_plot$peak_proportion == 0
    pvalue <- t.test(
        data_for_plot$mutation_density[data_for_plot$has_peak == TRUE],
        data_for_plot$mutation_density[data_for_plot$has_peak == FALSE]
    )$p.value
    pvalue <- format(pvalue, digits = 2)

    ypos <- max(data_for_plot$mutation_density) * 0.8

    p <- ggplot(data_for_plot, aes(x = has_peak, y = mutation_density)) +
        geom_boxplot(aes(color = has_peak)) +
        ggtitle(title) +
        labs(y = "Mutation density") +
        annotate("text", x = 1.5, y = ypos, label = pvalue, size = 12) +
        # stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method = "t.test", size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        scale_x_discrete(labels = c("Without peak", "With peak")) +
        scale_color_manual(values = c(rgb(137, 157, 159, maxColorValue = 255), rgb(236, 174, 34, maxColorValue = 255))) +
        scale_y_log10() +
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

plot_boxplot <- function(data_for_plot = NULL) {
    # data_for_plot$site_count_class <- sapply(data_for_plot$site_count_density, function(x) {
    #     unit <- 0.00025
    #     return(x)
    # })

    unit <- 0.0001

    data_for_plot$site_count_class <- ceiling(data_for_plot$site_count_density / unit)
    # data_for_plot$site_count_class[data_for_plot$site_count_class > 12] <- 12
    # data_for_plot$site_count_class[data_for_plot$site_count_class <=3] <- 3
    data_for_plot$site_count_class <- factor(data_for_plot$site_count_class,
        levels = seq(1, max(data_for_plot$site_count_class))
    )
    p <- ggplot(data_for_plot, aes(x = site_count_class, y = mutation_density)) +
        geom_boxplot() +
        geom_point(aes(color = gene_class)) +
        # scale_x_discrete(labels = c("[0,0.0005]", "(0.0005,0.001]", "(0.001,0.0015]", "(0.0015,0.002]", "(0.002,0.0025]", "(0.0025,+∞]")) +
        labs(x = "Site count density", y = "Mutation density") +
        # scale_y_continuous(trans = "log10") +
        theme_bw()
    # theme(
    #     strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
    #     strip.text.x = element_text(size = 30, color = "white"),
    #     panel.border = element_blank(),
    #     panel.grid = element_blank(),
    #     axis.line = element_line(colour = "black"),
    #     axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
    #     # axis.text.x = element_text(size = 36, color = "black"),
    #     axis.text.y = element_text(size = 36, color = "black"),
    #     plot.title = element_text(
    #         colour = "black",
    #         size = 40, vjust = 0.5, hjust = 0.5
    #     ),
    #     axis.title.x = element_text(size = 40, color = "black"),
    #     axis.title.y = element_text(size = 40, color = "black"),
    #     legend.title = element_blank(),
    #     legend.text = element_text(size = 24, color = "black"),
    #     legend.position = "bottom",
    #     plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
    #     panel.spacing = unit(3, "lines")
    # )
    return(p)
}

plot_boxplot_50part <- function(data_for_plot = NULL) {
    data_for_plot_ordered <- data_for_plot[order(data_for_plot$site_count_density, decreasing = F), ]
    unit <- ceiling(nrow(data_for_plot_ordered) / 20)
    data_for_plot_ordered$site_count_class <- (1:nrow(data_for_plot_ordered)) %/% unit + 1
    data_for_plot_ordered$site_count_class <- factor(data_for_plot_ordered$site_count_class,
        levels = c(1:max(data_for_plot_ordered$site_count_class))
    )
    p <- ggplot(data_for_plot_ordered, aes(x = site_count_class, y = mutation_density)) +
        geom_boxplot() +
        geom_point(aes(color = gene_class)) +
        # scale_x_discrete(labels = c("[0,0.0005]", "(0.0005,0.001]", "(0.001,0.0015]", "(0.0015,0.002]", "(0.002,0.0025]", "(0.0025,+∞]")) +
        labs(x = "Site count class", y = "Mutation density") +
        # scale_y_continuous(trans = "log10") +
        theme_bw()
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
for (gene_class_criterion in c("core", "h3k9me3")) {
    data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion)

    data_core <- data %>% filter(gene_class == "core")
    data_noncore <- data %>% filter(gene_class == "noncore")

    core_p <- plot_scatter(data_for_plot = data_core %>% filter(peak_proportion != 0) %>% filter(length >= 300)) +
        ggtitle(ifelse(gene_class_criterion == "core", "Core genome", "Non-H3K9me3 region"))
    noncore_p <- plot_scatter(data_for_plot = data_noncore %>% filter(peak_proportion != 0) %>% filter(length >= 300)) +
        ggtitle(ifelse(gene_class_criterion == "core", "Non-core genome", "H3K9me3 region"))

    noncore_p$labels$y <- ""
    core_p | noncore_p
    ggsave(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome/correlation_between_peak_proportion_with_MD_whole_genome_flank_", gene_class_criterion, ".pdf"), width = 16, height = 8)
}

for (gene_class_criterion in c("core", "h3k9me3")) {
    data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion)
    data_core <- data %>% filter(gene_class == "core")
    data_noncore <- data %>% filter(gene_class == "noncore")
    whole_p_youwu6mA <- plot_boxplot_youwu6mA(data_for_plot = data, title = "Whole genome")
    core_p_youwu6mA <- plot_boxplot_youwu6mA(data_for_plot = data_core, title = ifelse(gene_class_criterion == "core", "Core genome", "Non-H3K9me3 region"))
    noncore_p_youwu6mA <- plot_boxplot_youwu6mA(data_for_plot = data_noncore, title = ifelse(gene_class_criterion == "core", "Non-core genome", "H3K9me3 region"))
    core_p_youwu6mA$labels$y <- ""
    noncore_p_youwu6mA$labels$y <- ""

    whole_p_youwu6mA | core_p_youwu6mA | noncore_p_youwu6mA
    ggsave(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome/comparing_gene_with_and_without_peak_", gene_class_criterion, ".pdf"), height = 8, width = 24)
}

## mitochondria_1里的基因数要比mitochondria_2里的多很多
mito_gene_list <- unique(unlist(read.table(paste0(peak_proportion_global_dir, "/mitochondria_1"), as.is = T, header = F)))
# mito_gene_list <- unique(unlist(read.table(paste0(peak_proportion_global_dir, "/mitochondria_2"), as.is = T, header = F)))
export_gene_list <- unique(unlist(read.table(paste0(peak_proportion_global_dir, "/export_gene_list"), as.is = T, header = F)))

invasion_gene_list <- unique(unlist(read.table(paste0(peak_proportion_global_dir, "/invasion"), as.is = T, header = F)))

mito_data <- data[match(mito_gene_list, data$gene_name), ]
export_data <- data[match(export_gene_list, data$gene_name), ]
invasion_data <- data[match(invasion_gene_list, data$gene_name), ]
mito <- plot_scatter(mito_data) + ggtitle("Mitochondria") + scale_x_continuous(limits = c(0, 0.6))
export <- plot_scatter(export_data) + ggtitle("Export protein") + scale_x_continuous(limits = c(0, 0.6))
invasion <- plot_scatter(invasion_data) + ggtitle("Invasion") + scale_x_continuous(limits = c(0, 0.6))
mito | invasion
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome/correlation_between_peak_proportion_with_MD_mitochodria_invasion_flank.pdf", width = 16, height = 8)
write.table(mito_data, "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/mito_mutation_density.txt", quote = F, sep = "\t", row.names = F)
write.table(invasion_data, "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/invasion_mutation_density.txt", quote = F, sep = "\t", row.names = F)


## 3代数据6mA位点密度与基因突变密度的关系
site_density <- read_6mA_site_density()
data_with_site_density <- merge(data, site_density, by = "gene_name")
w <- plot_scatter1(data_for_plot = data_with_site_density) + ggtitle("Whole genome")
h <- plot_scatter1(data_for_plot = data_with_site_density %>% filter(gene_class == "noncore")) + ggtitle("heterochromatin") + scale_x_continuous(limits = c(0, 0.0025))
e <- plot_scatter1(data_for_plot = data_with_site_density %>% filter(gene_class == "core")) + ggtitle("Euchromatin")


plot_boxplot(data_for_plot = data_with_site_density %>% filter(gene_class == "noncore")) + ggtitle("heterochromatin")
plot_boxplot(data_for_plot = data_with_site_density) + ggtitle("Whole genome")
plot_boxplot(data_for_plot = data_with_site_density %>% filter(gene_class == "core")) + ggtitle("Euchromatin")


w <- plot_boxplot_50part(data_for_plot = data_with_site_density) + ggtitle("Whole genome")
w_limits <- plot_boxplot_50part(data_for_plot = data_with_site_density) + ggtitle("Whole genome") + scale_y_continuous(limits = c(0.015, 0.03))
w | w_limits

e <- plot_boxplot_50part(data_for_plot = data_with_site_density %>% filter(gene_class == "core")) + ggtitle("Euchromatin")
e_limits <- plot_boxplot_50part(data_for_plot = data_with_site_density %>% filter(gene_class == "core")) + ggtitle("Euchromatin") + scale_y_continuous(limits = c(0.015, 0.03))
e | e_limits

plot_boxplot_50part(data_for_plot = data_with_site_density %>% filter(gene_class == "noncore")) + ggtitle("heterochromatin")
