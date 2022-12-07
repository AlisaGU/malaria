#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

# 2. functions ------------------------------------------------------------ TODO:
plot_motif_length_ratio <- function(peak = NULL, nopeak = NULL) {
    peak_index <- which(peak$window_size == window_size)
    peak_50 <- peak[peak_index, ]
    nopeak_index <- which(nopeak$window_size == window_size)
    nopeak_50 <- nopeak[nopeak_index, ]

    peak_motif_base_count_in_window <- as.data.frame(table(peak_50$base_count_in_motif))
    peak_motif_base_count_in_window$proportion <- peak_motif_base_count_in_window$Freq / sum(peak_motif_base_count_in_window$Freq)

    nopeak_motif_base_count_in_window <- as.data.frame(table(nopeak_50$base_count_in_motif))
    nopeak_motif_base_count_in_window$proportion <- nopeak_motif_base_count_in_window$Freq / sum(nopeak_motif_base_count_in_window$Freq)

    merged_data <- merge(peak_motif_base_count_in_window, nopeak_motif_base_count_in_window, by = "Var1", all = T)
    colnames(merged_data) <- c("motif_len_in_window", "Freq.peak", "proportion.peak", "Freq.nopeak", "proportion.nopeak")

    ratio <- merged_data$proportion.peak / merged_data$proportion.nopeak
    names(ratio) <- seq(0, window_size)
    data_for_plot <- data.frame("motif_len_in_window" = 0:window_size, "ratio" = ratio)
    data_for_plot$class <- sapply(data_for_plot$ratio, function(x) {
        if (is.na(x)) {
            return("peak")
        } else {
            return("both")
        }
    })
    data_for_plot$ratio <- sapply(data_for_plot$ratio, function(x) {
        if (is.na(x)) {
            return(max(data_for_plot$ratio, na.rm = TRUE) + 10)
        } else {
            return(x)
        }
    })
    p <- ggplot(data_for_plot %>% filter(class == "both"), aes(x = motif_len_in_window, y = ratio)) +
        geom_point(size = 4, color = rgb(204, 161, 47, maxColorValue = 255)) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.25, label.y.npc = 0.1, size = 12) +
        scale_x_continuous(limits = c(0, 50)) +
        scale_y_continuous(trans = "log10") +
        labs(x = "Motif length per 50 bp", y = "Window proportion ratio\n(6mA vs control)") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            # axis.title.x = element_blank(),

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

plot_motif_mutation_count <- function(peak = NULL) {
    peak_index <- which(peak$window_size == window_size)
    peak_50 <- peak[peak_index, ]

    snp <- peak_50 %>%
        group_by(base_count_in_motif) %>%
        summarise(mutation_count = mean(snp_mutation_count))

    indel <- peak_50 %>%
        group_by(base_count_in_motif) %>%
        summarise(mutation_count = mean(indel_mutation_count))

    data_for_plot <- rbind(data.frame(snp, type = "SNP"), data.frame(indel, type = "Indel"))
    data_for_plot$mutation_density <- data_for_plot$mutation_count / 50
    data_for_plot$type <- factor(data_for_plot$type, levels = c("SNP", "Indel"))

    p_snp <- ggplot(data_for_plot %>% filter(type == "SNP"), aes(x = base_count_in_motif, y = mutation_density, color = type, group = type)) +
        geom_point(size = 4) +
        scale_x_continuous(limits = c(0, 50)) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.9, size = 12) +
        scale_color_manual(values = c("SNP" = rgb(190, 66, 65, maxColorValue = 255), "Indel" = rgb(29, 113, 169, maxColorValue = 255))) +
        labs(x = "Motif length per 50 bp", y = "Mutation density") +
        theme_bw() +
        ggtitle("SNP") +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_text(size = 40, color = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            # axis.title.x = element_blank(),

            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            # legend.position = "bottom", legend.direction = "horizontal",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )

    p_indel <- ggplot(data_for_plot %>% filter(type == "Indel"), aes(x = base_count_in_motif, y = mutation_density, color = type, group = type)) +
        geom_point(size = 4) +
        scale_x_continuous(limits = c(0, 50)) +
        stat_smooth(method = "lm", col = "black", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.9, size = 12) +
        scale_color_manual(values = c("SNP" = rgb(190, 66, 65, maxColorValue = 255), "Indel" = rgb(29, 113, 169, maxColorValue = 255))) +
        labs(x = "Motif length per 50 bp", y = "Mutation density") +
        theme_bw() +
        ggtitle("Indel") +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_text(size = 40, color = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            # axis.title.x = element_blank(),

            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            # legend.position = "bottom", legend.direction = "horizontal",
            legend.position = "none",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    result <- p_snp | p_indel
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
window_size <- 50
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_50_core_genome_two_outgroup")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak <- fread("allChrom.window.peak.filter.summary.gz", header = F, stringsAsFactors = F, sep = "\t")
nopeak <- fread("allChrom.window.nopeak.filter.summary.gz", header = F, stringsAsFactors = F, sep = "\t")
colnames(peak) <- colnames(nopeak) <- c(
    "chrom", "start_minus_1", "end",
    "base_count_in_motif", "snp_mutation_count", "indel_mutation_count"
)

peak$window_size <- peak$end - peak$start_minus_1
nopeak$window_size <- nopeak$end - nopeak$start_minus_1

ratio <- plot_motif_length_ratio(peak = peak, nopeak = nopeak)
mutation <- plot_motif_mutation_count(peak = peak)
ratio + mutation + plot_layout(widths = c(1, 2))
ggsave("window_proportion_ratio_and_mutation_density.pdf", width = 27, height = 10)
