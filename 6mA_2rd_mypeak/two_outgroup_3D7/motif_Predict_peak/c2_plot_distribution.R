#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
window_size <- 50
setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_", window_size))

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak <- fread("allChrom.window.peak.filter.summary.gz", header = F, stringsAsFactors = F)
nopeak <- fread("allChrom.window.nopeak.awayFromPeak.filter.summary.gz", header = F, stringsAsFactors = F)
colnames(peak) <- colnames(nopeak) <- c("chrom", "start_minus_1", "end", "base_count_in_motif", "motif_count_in_window", "maxFE", "meanFE", "medianFE")

{
    a <- peak$base_count_in_motif
}



peak_motif_base_count_in_per_base <- peak$base_count_in_motif / (peak$end - peak$start_minus_1)
nopeak_motif_base_count_in_per_base <- nopeak$base_count_in_motif / (nopeak$end - nopeak$start_minus_1)
data_for_plot <- as.data.frame(rbind(cbind(value = peak_motif_base_count_in_per_base, type = "peak"), cbind(nopeak_motif_base_count_in_per_base, type = "nopeak")))
data_for_plot$value <- as.numeric(data_for_plot$value)

ggplot(data_for_plot, aes(x = value, group = type, color = type)) +
    geom_vline(xintercept = 6 / window_size, color = "blue") +
    geom_vline(xintercept = 9 / window_size, color = "blue") +
    geom_vline(xintercept = 12 / window_size, color = "blue") +
    geom_vline(xintercept = 15 / window_size, color = "blue") +
    geom_vline(xintercept = 18 / window_size, color = "blue") +
    geom_vline(xintercept = 21 / window_size, color = "blue") +
    geom_vline(xintercept = 24 / window_size, color = "blue") +
    geom_vline(xintercept = 27 / window_size, color = "blue") +
    geom_density(size = 0.5) +
    scale_x_continuous(
        limits = c(0, 42 / window_size), labels = c("6", "9", "12", "15", "18", "21", "24", "27"),
        # breaks = c(0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135),
        breaks = c(6, 9, 12, 15, 18, 21, 24, 27) / window_size
    ) +
    labs(x = "Average motif base count in per site") +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        panel.border = element_blank(),
        # panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_text(size = 36, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
        axis.text.x = element_text(size = 36, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_text(size = 40, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        legend.position = "bottom", legend.direction = "horizontal",
        # legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave("motif_base_count_distribution.pdf", width = 10, height = 7)





# peak_split <- split(peak, f = peak$chrom, drop = T)
# mean(sapply(peak_split, function(x) {
#     sum(x$base_count_in_motif) / sum(x$end - x$start_minus_1)
# }))
# mean(sapply(peak_split, function(x) {
#     sum(x$motif_count_in_window) / sum(x$end - x$start_minus_1)
# }))

# nopeak_split <- split(nopeak, f = nopeak$chrom, drop = T)
# mean(sapply(nopeak_split, function(x) {
#     sum(x$base_count_in_motif) / sum(x$end - x$start_minus_1)
# }))
# mean(sapply(nopeak_split, function(x) {
#     sum(x$motif_count_in_window) / sum(x$end - x$start_minus_1)
# }))

peak_motif_count_in_per_base <- peak$motif_count_in_window / (peak$end - peak$start_minus_1)
nopeak_motif_count_in_per_base <- nopeak$motif_count_in_window / (nopeak$end - nopeak$start_minus_1)
data_for_plot <- as.data.frame(rbind(cbind(value = peak_motif_count_in_per_base, type = "peak"), cbind(nopeak_motif_count_in_per_base, type = "nopeak")))
data_for_plot$value <- as.numeric(data_for_plot$value)
ggplot(data_for_plot, aes(x = value, group = type, color = type)) +
    geom_vline(xintercept = 1 / window_size, color = "blue") +
    geom_vline(xintercept = 2 / window_size, color = "blue") +
    geom_vline(xintercept = 3 / window_size, color = "blue") +
    geom_vline(xintercept = 4 / window_size, color = "blue") +
    geom_vline(xintercept = 5 / window_size, color = "blue") +
    geom_vline(xintercept = 6 / window_size, color = "blue") +
    geom_vline(xintercept = 7 / window_size, color = "blue") +
    geom_density(size = 0.2, alpha = 0.3) +
    scale_x_continuous(
        limits = c(0, 9 / window_size),
        breaks = c(1, 2, 3, 4, 5, 6, 7) / window_size,
        labels = c(1, 2, 3, 4, 5, 6, 7)
    ) +
    labs(x = "Average motif count in per site") +
    theme_bw() +
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
        axis.title.x = element_text(size = 40, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, color = "black"),
        legend.position = "bottom", legend.direction = "horizontal",
        # legend.position = "none",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave("motif_count_distribution.pdf", width = 10, height = 7)
