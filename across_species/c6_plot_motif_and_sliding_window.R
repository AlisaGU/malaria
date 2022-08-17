#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# 2. functions ------------------------------------------------------------ TODO:
read_motif_data <- function() {
    files <- "snps_mean_variant_count"
    data <- read.table(files, as.is = T)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    colnames(data1) <- c("peak_motif", "whole_peak_region", "control_motif", "whole_control_region")
    rownames(data1) <- paste("chrom", 1:14, sep = "")

    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    a <- gsub("chrom", "", data_for_plot$chrom)
    a[a != "1"] <- ""
    data_for_plot$label <- a
    data_for_plot$variant <- "SNPs"
    data_for_plot$source <- factor(data_for_plot$source, levels = c("peak_motif", "whole_peak_region", "control_motif", "whole_control_region"))
    return(data_for_plot)
}

plot_motif <- function(data_for_plot = NULL) {
    color_value_N2N <- c("control_motif" = "black", "peak_motif" = "#b2182b", "whole_peak_region" = "#d75665", "whole_control_region" = "grey")
    p <- ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 8) +
        scale_color_manual(values = color_value_N2N, labels = c("Peak motif", "Whole peak region", "Control motif", "Whole control region")) +
        scale_x_discrete(labels = c("Peak motif", "Whole peak region", "Control motif", "Whole control region")) +
        stat_compare_means(
            comparisons = list(c("peak_motif", "whole_peak_region"), c("peak_motif", "control_motif"), c("control_motif", "whole_control_region"), c("peak_motif", "whole_control_region")),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 14
        ) +
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
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

read_sliding_window_data <- function() {
    # result <- lapply(c("snps", "indels"), function(variant) {
    result <- lapply(c("indels"), function(variant) {
        files <- paste0(variant, "_sliding_window_mean_variant_count_LI_2")
        # files <- paste0(variant, "_sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_",variant)

        data <- read.table(files, as.is = T)
        data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
        colnames(data1) <- c(paste("w", 1:10, sep = ""), "control")
        rownames(data1) <- paste("chrom", 1:14, sep = "")


        data_for_plot <- data.frame(
            value = unlist(data1),
            chrom = rep(rownames(data1), times = ncol(data1)),
            source = rep(colnames(data1), each = nrow(data1))
        )
        a <- gsub("chrom", "", data_for_plot$chrom)
        a[a != "1"] <- ""
        data_for_plot$label <- a
        data_for_plot$variant <- ifelse(variant == "snps", "SNPs", "Indels")
        data_for_plot$source <- factor(data_for_plot$source, levels = c("control", paste("w", 1:10, sep = "")))
        return(data_for_plot)
    })
    result <- do.call(rbind, result)
    result$variant <- factor(result$variant, levels = c("SNPs", "Indels"))
    return(result)
}


plot_sliding_window <- function(data_for_plot = NULL, title = "") {
    cols <- brewer.pal(3, "YlOrRd")
    pal <- colorRampPalette(cols)
    mycolors <- rev(pal(10))
    color_value <- c("#868686", mycolors)
    names(color_value) <- c("control", paste("w", 1:10, sep = ""))
    p <- ggplot(data_for_plot, aes(x = source, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        facet_wrap(~variant, scales = "free_y") +
        scale_color_manual(values = color_value) +
        scale_x_discrete(labels = c("Control", "Submit~5%", "5~10%", "10~15%", "15~20%", "20~25%", "25~30%", "30~35%", "35~40%", "40~45%", "45~50%")) +
        # stat_compare_means(
        #     # ref = "control",
        #     comparisons = list(c("control","w1")),
        #      size = 14, paired = T,
        #     hide.ns = T, vjust = 0.5, show.legend = FALSE,
        #     label = "p.signif",
        #     # aes(label = paste0("p = ", ..p.format..)),
        #     # paired = T, size = 3
        # ) +
        stat_compare_means(
            method = "t.test",
            tip.length = 0, bracket.size = 0,
            comparisons = list(c("control", "w1")),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p = {p.format}", paired = T, size = 10
        ) +
        scale_y_continuous(expand = expansion(mult = 0.12)) +
        labs(y = "Mutation density", title = title) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black", face = "italic",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_reichenowi/mutation_density")
data_for_plot <- read_motif_data()
plot_motif(data_for_plot = data_for_plot)
ggsave("motif_3_parts_comparison.pdf", width = 16, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_reichenowi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot)
ggsave("species3_sliding_window_comparison.pdf", width = 10, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/P_falciparum_only/P_praefalciparum/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "P.falciparum - P.praefalciparum")
ggsave("P_praefalciparum_sliding_window_comparison.pdf", width = 10, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/P_falciparum_only/P_reichenowi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "P.falciparum  - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison.pdf", width = 20, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison.pdf", width = 20, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock.pdf", width = 20, height = 7)

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock_bedFirst.pdf", width = 20, height = 7)

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock_bedFirst_precise.pdf", width = 20, height = 7)

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data() # _sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_indels
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock_bedFirst_precise_indel_excluding_P_bill_gap.pdf", width = 20, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data()
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock_bedFirst_precise_indel_two_outgroup_nongap.pdf", width = 20, height = 7)

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data() # _sliding_window_mean_variant_count_noFilter_include_indel_bedFirst_precise_indel_twoOutgroup_nongap_indels_both
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_comparison_noFilterBlock_bedFirst_precise_indel_two_outgroup_nongap_both.pdf", width = 10, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data() # _sliding_window_mean_variant_count_LI
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_LI.pdf", width = 10, height = 7)


setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data() # _sliding_window_mean_variant_count_LI_1
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_LI_1.pdf", width = 10, height = 7)

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/mutation_density")
data_for_plot <- read_sliding_window_data() # _sliding_window_mean_variant_count_LI_2
plot_sliding_window(data_for_plot = data_for_plot, title = "(P.falciparum + P.billcollinsi) - P.reichenowi")
ggsave("P_reichenowi_sliding_window_LI_2.pdf", width = 10, height = 7)
