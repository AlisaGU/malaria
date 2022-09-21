#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        mutation_density <- read.table(paste0(mutation_density_global_dir, "/", gene_set, "/snps_mean_variant_count_gene"), as.is = T, stringsAsFactors = F)
        colnames(mutation_density) <- c("gene_name", "mutation_count", "cds_length")
        mutation_density_CDS <- data.frame(
            gene_name = mutation_density$gene_name,
            mutation_density = mutation_density$mutation_count / mutation_density$cds_length,
            gene_set = gene_set
        )


        peak_length <- read.table(paste0(peak_proportion_global_dir, "/", gene_set, "/", gene_set, "_peak_control_length_in_gene"), as.is = T, stringsAsFactors = F)
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
        geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12, color = "red") +
        labs(x = "Proportion of peak in gene", y = "Mutation density") +
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
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read_peak_proportion_and_mutation_density()

data_core <- data %>% filter(gene_set == "HDR" | gene_set == "DR" | gene_set == "RNA_translation")
data_noncore <- data %>% filter(gene_set == "RIF" | gene_set == "STEVOR" | gene_set == "VAR")

core_p <- plot_scatter(data_for_plot = data_core)
noncore_p <- plot_scatter(data_for_plot = data_noncore)
noncore_p$labels$y <- ""
core_p | noncore_p
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets_exclude_pseudo_with_fold_enrichment_including_RIF/correlation_between_peak_proportion_with_MD_gene.pdf", width = 16, height = 8)
