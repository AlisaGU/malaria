#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function() {
    mutation_density <- read.table(paste0(mutation_density_global_dir, "/snps_mean_variant_count_CDS"), as.is = T, stringsAsFactors = F)
    colnames(mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    mutation_density_CDS <- data.frame(
        gene_name = mutation_density$gene_name,
        mutation_density = mutation_density$mutation_count / mutation_density$cds_length
    )


    peak_length <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_CDS"), as.is = T, stringsAsFactors = F)
    colnames(peak_length) <- c("gene_name", "peak_length", "control_length")
    peak_proportion <- data.frame(gene_name = peak_length$gene_name, peak_proportion = apply(peak_length[, 2:3], 1, function(x) {
        x[1] / sum(x)
    }))

    gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)[, c(4, 6)]
    colnames(gene_class) <- c("gene_name", "gene_class")
    result <- merge(merge(mutation_density_CDS, peak_proportion, by = "gene_name"), gene_class, by = "gene_name")
    return(result)
}

plot_scatter <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density)) +
        geom_point(size = 2) +
        # geom_point(size = 8) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12, color = "red") +
        labs(x = "Proportion of peak in CDS", y = "Mutation density") +
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
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read_peak_proportion_and_mutation_density()

data_core <- data %>% filter(gene_class == "core")
data_noncore <- data %>% filter(gene_class == "noncore")

core_p <- plot_scatter(data_for_plot = data_core)
noncore_p <- plot_scatter(data_for_plot = data_noncore)
noncore_p$labels$y <- ""
core_p | noncore_p
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/wholeGenome/correlation_between_peak_proportion_with_MD_whole_genome_CDS.pdf", width = 16, height = 8)
