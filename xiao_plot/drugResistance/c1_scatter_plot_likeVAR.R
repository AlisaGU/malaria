#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(data.table)
library(patchwork)
library(grid)
library(gg.gap)
library(ggbreak)
library(ggrepel)
library(dplyr)

# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function(flank_length = NULL) {
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


    peak_length <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_flank", flank_length, "_peak25Perc"), as.is = T, stringsAsFactors = F)
    colnames(peak_length) <- c("gene_name", "peak_length", "control_length")
    peak_proportion <- data.frame(
        gene_name = peak_length$gene_name,
        peak_proportion = apply(peak_length[, 2:3], 1, function(x) {
            x[1] / sum(x)
        })
    )

    gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(merge(mutation_density_CDS, peak_proportion, by = "gene_name"), gene_class, by = "gene_name")
    return(result)
}

plot_scatter_just_by_mutationDensity <- function(data = NULL) {
    data$index <- 1:nrow(data)
    p <- ggplot(data, aes(x = index, y = mutation_density, color = peak_proportion)) +
        geom_text_repel(aes(label = gene_name)) +
        geom_point(size = 5) +
        scale_color_gradient(
            low = "grey", high = "#aed477"
        ) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            # legend.title = element_blank(),
            # legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_scatter_just_by_MD_and_gene_6mAdensity <- function(data = NULL, title = NULL) {
    formula <- y ~ poly(x, 2, raw = TRUE)
    p <- ggplot(data, aes(x = gene_density_signal, y = mutation_density, color = peak_proportion)) +
        geom_text_repel(aes(label = gene_name)) +
        geom_point(size = 5) +
        labs(title = title) +
        # geom_smooth(
        #    method = 'lm', formula = y~poly(x, 2),se=FALSE
        # ) +
        # stat_cor( method = "pearson") +
        stat_poly_line(formula = formula,se=FALSE) +
        stat_poly_eq(
            mapping = use_label(c("eq", "adj.R2")),
            formula = formula
        ) +
        scale_color_gradient(
            low = "grey", high = "#aed477"
        ) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            # legend.title = element_blank(),
            # legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_scatter_just_by_MD_and_gene_flank2kb_6mAdensity <- function(data = NULL) {
    p <- ggplot(data, aes(x = gene_flank2kb_density_signal, y = mutation_density, color = peak_proportion)) +
        geom_text_repel(aes(label = gene_name)) +
        geom_point(size = 5) +
        scale_color_gradient(
            low = "grey", high = "#aed477"
        ) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            # legend.title = element_blank(),
            # legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/mutation_density_distribution/wholeGenome/ESEA_WSEA_OCE_SAM_SAS/allgenes/3D7"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/mutation_density_distribution/wholeGenome/ESEA_WSEA_OCE_SAM_SAS/allgenes/3D7"

WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"

WT_chip_depth <- 36971023
WT_input_depth <- 21500156


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance")
gene_list <- read.table("drugResistanceMarker", stringsAsFactors = F, sep = "\t", header = F)
gene_6mAdensity <- read.table("drug_resistance_gene_6mADensity.txt", as.is = T, header = T, sep = "\t")
gb_flank2kb_6mAdensity <- read.table("drug_resistance_gene_flank2kb_6mADensity.txt", as.is = T, header = T, sep = "\t")
colnames(gb_flank2kb_6mAdensity)[1] <- colnames(gene_6mAdensity)[1] <- colnames(gene_list) <- "gene_name"

data <- read_peak_proportion_and_mutation_density(flank_length = "0kb")
data_for_plot <- merge(merge(merge(data, gene_list, by = "gene_name"), gene_6mAdensity, by = "gene_name"), gb_flank2kb_6mAdensity, by = "gene_name")


snp_data <- data_for_plot %>% filter(type == "SNP")
snp_data_order <- snp_data[order(snp_data$mutation_density, decreasing = T), ]
plot_scatter_just_by_mutationDensity(data = snp_data_order)
plot_scatter_just_by_MD_and_gene_6mAdensity(data = snp_data_order, title = "SNP")
plot_scatter_just_by_MD_and_gene_flank2kb_6mAdensity(data = snp_data_order)

indel_data <- data_for_plot %>% filter(type == "Indel")
indel_data_order <- indel_data[order(indel_data$mutation_density, decreasing = T), ]
plot_scatter_just_by_mutationDensity(data = indel_data_order)
plot_scatter_just_by_MD_and_gene_6mAdensity(data = indel_data_order, title = "Indel")
plot_scatter_just_by_MD_and_gene_flank2kb_6mAdensity(data = indel_data_order)
