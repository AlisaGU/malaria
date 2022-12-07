#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        mutation_density <- read.table(paste0(mutation_density_global_dir, "/", gene_set, "/snps_mean_variant_count_flank"), as.is = T, stringsAsFactors = F)
        colnames(mutation_density) <- c("gene_name", "mutation_count", "cds_length")
        mutation_density_CDS <- data.frame(
            gene_name = mutation_density$gene_name,
            mutation_density = mutation_density$mutation_count / mutation_density$cds_length,
            gene_set = gene_set
        )


        peak_length <- read.table(paste0(peak_proportion_global_dir, "/", gene_set, "/", gene_set, "_peak_control_length_in_flank"), as.is = T, stringsAsFactors = F)
        colnames(peak_length) <- c("gene_name", "peak_length", "control_length")
        peak_proportion <- data.frame(gene_name = peak_length$gene_name, peak_proportion = apply(peak_length[, 2:3], 1, function(x) {
            x[1] / sum(x)
        }))

        mutation_density_gene <- read.table(paste0(mutation_density_global_dir, "/", gene_set, "/snps_mean_variant_count_gene"), as.is = T, stringsAsFactors = F)
        peak_fold_enrichment <- mutation_density_gene[, c(1, 6)]
        colnames(peak_fold_enrichment) <- c("gene_name", "fold_enrichment")
        result <- merge(merge(mutation_density_CDS, peak_proportion, by = "gene_name"), peak_fold_enrichment, by = "gene_name")
        result$fold_enrichment[is.na(result$fold_enrichment)] <- 0
        return(result)
    })
    result <- do.call(rbind, result)
    return(result)
}

plot_scatter <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density)) +
        geom_point(aes(color = gene_set), size = 8) +
        # geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12, color = "red") +
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

plot_scatter_nocolor <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density)) +
        geom_point(size = 8, color = rgb(250, 217, 126, maxColorValue = 255), alpha = 0.5) +
        # geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12, color = "red") +
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
# 3. input ---------------------------------------------------------------- TODO:
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_exclude_pseudo"
gene_sets <- c("HDR", "RIF", "RNA_translation", "STEVOR", "VAR")

mito_data <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/mito_mutation_density.txt", header = T, as.is = T)
invasion_data <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/invasion_mutation_density.txt", header = T, as.is = T)

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read_peak_proportion_and_mutation_density()

data_core <- data %>% filter(gene_set == "HDR" | gene_set == "DR" | gene_set == "RNA_translation")
data_noncore <- data %>% filter(gene_set == "RIF" | gene_set == "STEVOR" | gene_set == "VAR")

core_p <- plot_scatter(data_for_plot = data_core)
noncore_p <- plot_scatter(data_for_plot = data_noncore)
noncore_p$labels$y <- ""
core_p | noncore_p
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/correlation_between_peak_proportion_with_MD_flank.pdf", width = 16, height = 8)


noncore_p_nocolor <- plot_scatter_nocolor(data_for_plot = data_noncore)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/correlation_between_peak_proportion_with_MD_flank_var_rifin_stevor.pdf", width = 8, height = 8)
# plot_scatter(data_for_plot = data %>% filter(gene_set == "VAR")) + scale_color_manual(values = c("VAR" = "#faca5d"))
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/correlation_between_peak_proportion_with_MD_flank_VAR.pdf", width = 8, height = 8)


data <- rbind(
    data[, c("gene_name", "mutation_density", "gene_set", "peak_proportion")],
    data.frame(gene_name = mito_data$gene_name, mutation_density = mito_data$mutation_density, gene_set = "Mitochondria", peak_proportion = mito_data$peak_proportion),
    data.frame(gene_name = invasion_data$gene_name, mutation_density = invasion_data$mutation_density, gene_set = "Invasion", peak_proportion = invasion_data$peak_proportion)
)
ggplot(data, aes(x = peak_proportion, y = mutation_density, group = gene_set, color = gene_set)) +
    geom_point(size = 2, alpha = 0.3) +
    stat_smooth(aes(color = gene_set), method = "lm", se = FALSE) +
    stat_cor(aes(color = gene_set), method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12) +
    scale_color_manual() +
    labs(x = "Proportion of peak\nin gene and flank 2kb", y = "Mutation density") +
    theme_bw() +
    # scale_y_continuous(trans = "log10") +
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
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets/each_set_correlation_log.pdf", height = 12, width = 12)
