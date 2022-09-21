#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
plot_motif_length <- function(peak = NULL, nopeak = NULL) {
    peak_motif_base_count_in_per_base <- peak$base_count_in_motif / (peak$end - peak$start_minus_1)
    nopeak_motif_base_count_in_per_base <- nopeak$base_count_in_motif / (nopeak$end - nopeak$start_minus_1)
    data_for_plot <- as.data.frame(rbind(cbind(value = peak_motif_base_count_in_per_base, type = "peak"), cbind(nopeak_motif_base_count_in_per_base, type = "nopeak")))
    data_for_plot$value <- as.numeric(data_for_plot$value)

    p <- ggplot(data_for_plot, aes(x = value, group = type, color = type)) +
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
    return(p)
}

plot_motif_count <- function(peak = NULL, nopeak = NULL) {
    peak_motif_count_in_per_base <- peak$motif_count_in_window / (peak$end - peak$start_minus_1)
    nopeak_motif_count_in_per_base <- nopeak$motif_count_in_window / (nopeak$end - nopeak$start_minus_1)
    data_for_plot <- as.data.frame(rbind(cbind(value = peak_motif_count_in_per_base, type = "peak"), cbind(nopeak_motif_count_in_per_base, type = "nopeak")))
    data_for_plot$value <- as.numeric(data_for_plot$value)
    p <- ggplot(data_for_plot, aes(x = value, group = type, color = type)) +
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
    return(p)
}

plot_motif_length_ratio <- function(peak = NULL, nopeak = NULL, outfileName = NULL) {
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
    p <- ggplot(data_for_plot, aes(x = motif_len_in_window, y = ratio, color = class)) +
        geom_line(aes(group = 1), color = "black") +
        geom_point(size = 4) +
        scale_color_manual(values = c("both" = "black", "peak" = "red")) +
        labs(x = "Motif length per window", y = "Window proportion ratio\n(6mA vs control)") +
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
    # p1 <- p + scale_y_continuous(trans = "log10") + theme(axis.title.y = element_blank())


    # xlab <- p$labels$x
    # p1$labels$x <- p$labels$x <- " "

    # pdf(file = outfileName, width = 15, height = 7)
    # p | p1
    # grid::grid.draw(grid::textGrob(xlab, x = 0.6, y = 0.08, gp = grid::gpar(fontsize = 40, col = "black")))
    # dev.off()

    p1 <- p + scale_y_continuous(trans = "log10")

    ggsave(outfileName, width = 10, height = 10)
}

plot_motif_count_ratio <- function(peak = NULL, nopeak = NULL, outfileName = NULL) {
    peak_index <- which(peak$window_size == window_size)
    peak_50 <- peak[peak_index, ]
    nopeak_index <- which(nopeak$window_size == window_size)
    nopeak_50 <- nopeak[nopeak_index, ]

    peak_motif_count_in_window <- as.data.frame(table(peak_50$motif_count_in_window))
    peak_motif_count_in_window$proportion <- peak_motif_count_in_window$Freq / sum(peak_motif_count_in_window$Freq)

    nopeak_motif_count_in_window <- as.data.frame(table(nopeak_50$motif_count_in_window))
    nopeak_motif_count_in_window$proportion <- nopeak_motif_count_in_window$Freq / sum(nopeak_motif_count_in_window$Freq)

    merged_data <- merge(peak_motif_count_in_window, nopeak_motif_count_in_window, by = "Var1", all = T)
    colnames(merged_data) <- c("motif_count_in_window", "Freq.peak", "proportion.peak", "Freq.nopeak", "proportion.nopeak")

    ratio <- merged_data$proportion.peak / merged_data$proportion.nopeak
    names(ratio) <- merged_data$motif_count_in_window
    data_for_plot <- data.frame("motif_count_in_window" = as.numeric(names(ratio)), "ratio" = ratio)
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
    p <- ggplot(data_for_plot, aes(x = motif_count_in_window, y = ratio, color = class)) +
        geom_line(aes(group = 1), color = "black") +
        geom_point(size = 4) +
        labs(x = "Motif count per 50 bp", y = "Window proportion ratio of\n6mA over control") +
        scale_color_manual(values = c("both" = "black", "peak" = "red")) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            # axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
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
    # p1 <- p + scale_y_continuous(trans = "log10") + theme(axis.title.y = element_blank())


    # xlab <- p$labels$x
    # p1$labels$x <- p$labels$x <- " "

    # pdf(file = outfileName, width = 15, height = 7)
    # p | p1
    # grid::grid.draw(grid::textGrob(xlab, x = 0.6, y = 0.08, gp = grid::gpar(fontsize = 40, col = "black")))
    # dev.off()

    p1 <- p + scale_y_continuous(trans = "log10")
    ggsave(outfileName, width = 10, height = 10)
}
# 3. input ---------------------------------------------------------------- TODO:
window_size <- 50
setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_", window_size))

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak <- fread("allChrom.window.peak.filter.summary.gz", header = F, stringsAsFactors = F)
nopeak <- fread("allChrom.window.nopeak.awayFromPeak.filter.summary.gz", header = F, stringsAsFactors = F)
colnames(peak) <- colnames(nopeak) <- c("chrom", "start_minus_1", "end", "base_count_in_motif", "motif_count_in_window", "maxFE", "meanFE", "medianFE")

peak$window_size <- peak$end - peak$start_minus_1
nopeak$window_size <- nopeak$end - nopeak$start_minus_1


# plot_motif_length(peak = peak, nopeak = nopeak)
# ggsave("motif_base_count_distribution.pdf", width = 10, height = 7)

# plot_motif_count(peak = peak, nopeak = nopeak)
# ggsave("motif_count_distribution.pdf", width = 10, height = 7)


plot_motif_length_ratio(peak = peak, nopeak = nopeak, outfileName = "motif_length_ratio.pdf")

# plot_motif_count_ratio(peak = peak, nopeak = nopeak, outfileName = "motif_Count_ratio.pdf")
