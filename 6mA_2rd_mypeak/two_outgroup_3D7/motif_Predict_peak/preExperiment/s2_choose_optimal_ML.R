#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
# 2. functions ------------------------------------------------------------ TODO:
classfication_using_peak_base_count <- function(data = NULL, base_count_in_peak_threshold = NULL, base_count_in_motif_threshold = NULL) {
    data$in_peak <- as.numeric(data$base_count_in_peak >= base_count_in_peak_threshold)
    data$motif_base_satisfying_need <- as.numeric(data$base_count_in_motif >= base_count_in_motif_threshold)


    # result <- table(data$in_peak, data$motif_base_satisfying_need) / nrow(data)
    # # rownames(result) <- paste("peak_base_lg", base_count_in_peak_threshold, "_", rownames(result), sep = "")
    # # colnames(result) <- paste("motif_base_lg", base_count_in_motif_threshold, "_", colnames(result), sep = "")
    # rownames(result) <- c("real_nopeak", "real_peak")
    # colnames(result) <- c("predict_nopeak", "predict_peak")
    # return(result)
    ##                  predict_nopeak       predict_peak
    ## real_nopeak          a(TN)               b(FP)
    ## real_peak            c(FN)               d(TP)
    a <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 0))
    b <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 1))
    c <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 0))
    d <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 1))
    return(c(a, b, c, d))
}

classfication_using_peak_base_count_and_mean_FE <- function(data = NULL, base_count_in_peak_threshold = NULL, base_count_in_motif_threshold = NULL, meanFE_threshold = NULL) {
    data$in_peak <- as.numeric(data$base_count_in_peak >= base_count_in_peak_threshold & data$meanFE >= meanFE_threshold)
    data$motif_base_satisfying_need <- as.numeric(data$base_count_in_motif >= base_count_in_motif_threshold)
    ##                  predict_nopeak       predict_peak
    ## real_nopeak          a(TN)               b(FP)
    ## real_peak            c(FN)               d(TP)
    a <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 0))
    b <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 1))
    c <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 0))
    d <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 1))
    return(c(a, b, c, d))
}

classfication_using_peak_base_count_and_max_FE <- function(data = NULL, base_count_in_peak_threshold = NULL, base_count_in_motif_threshold = NULL, maxFE_threshold = NULL) {
    data$in_peak <- as.numeric(data$base_count_in_peak >= base_count_in_peak_threshold & data$maxFE >= maxFE_threshold)
    data$motif_base_satisfying_need <- as.numeric(data$base_count_in_motif >= base_count_in_motif_threshold)
    ##                  predict_nopeak       predict_peak
    ## real_nopeak          a(TN)               b(FP)
    ## real_peak            c(FN)               d(TP)
    a <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 0))
    b <- length(which(data$in_peak == 0 & data$motif_base_satisfying_need == 1))
    c <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 0))
    d <- length(which(data$in_peak == 1 & data$motif_base_satisfying_need == 1))
    return(c(a, b, c, d))
}

compute_TPR_FPR_using_peak_base_count <- function(data = NULL, base_count_in_peak_thresholds = NULL, base_count_in_motif_thresholds = NULL) {
    result <- lapply(base_count_in_peak_thresholds, function(base_count_in_peak_threshold) {
        a <- t(sapply(base_count_in_motif_thresholds, function(base_count_in_motif_threshold) {
            values <- classfication_using_peak_base_count(data = data, base_count_in_peak_threshold = base_count_in_peak_threshold, base_count_in_motif_threshold = base_count_in_motif_threshold)
            TP <- values[4]
            FN <- values[3]
            FP <- values[2]
            TN <- values[1]
            TPR <- TP / (TP + FN)
            FPR <- FP / (FP + TN)
            return(c(base_count_in_peak_threshold, base_count_in_motif_threshold, FPR, TPR))
        }))
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("base_count_in_peak_threshold", "base_count_in_motif_threshold", "FPR", "TPR")
    return(result)
}

compute_TPR_FPR_using_peak_base_count_and_mean_FE <- function(data = NULL, base_count_in_peak_thresholds = NULL, base_count_in_motif_thresholds = NULL, meanFE_thresholds = NULL) {
    result <- lapply(base_count_in_peak_thresholds, function(base_count_in_peak_threshold) {
        a <- lapply(base_count_in_motif_thresholds, function(base_count_in_motif_threshold) {
            a <- lapply(meanFE_thresholds, function(meanFE_threshold) {
                values <- classfication_using_peak_base_count_and_mean_FE(data = data, base_count_in_peak_threshold = base_count_in_peak_threshold, base_count_in_motif_threshold = base_count_in_motif_threshold, meanFE_threshold = meanFE_threshold)
                TP <- values[4]
                FN <- values[3]
                FP <- values[2]
                TN <- values[1]
                TPR <- TP / (TP + FN)
                FPR <- FP / (FP + TN)
                return(c(base_count_in_peak_threshold, base_count_in_motif_threshold, meanFE_threshold, FPR, TPR))
            })
            a <- do.call(rbind, a)
            return(a)
        })
        a <- do.call(rbind, a)
        return(a)
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("base_count_in_peak_threshold", "base_count_in_motif_threshold", "meanFE_threshold", "FPR", "TPR")
    return(result)
}

compute_TPR_FPR_using_peak_base_count_and_max_FE <- function(data = NULL, base_count_in_peak_thresholds = NULL, base_count_in_motif_thresholds = NULL, maxFE_thresholds = NULL) {
    result <- lapply(base_count_in_peak_thresholds, function(base_count_in_peak_threshold) {
        a <- lapply(base_count_in_motif_thresholds, function(base_count_in_motif_threshold) {
            a <- lapply(maxFE_thresholds, function(maxFE_threshold) {
                values <- classfication_using_peak_base_count_and_max_FE(data = data, base_count_in_peak_threshold = base_count_in_peak_threshold, base_count_in_motif_threshold = base_count_in_motif_threshold, maxFE_threshold = maxFE_threshold)
                TP <- values[4]
                FN <- values[3]
                FP <- values[2]
                TN <- values[1]
                TPR <- TP / (TP + FN)
                FPR <- FP / (FP + TN)
                return(c(base_count_in_peak_threshold, base_count_in_motif_threshold, maxFE_threshold, FPR, TPR))
            })
            a <- do.call(rbind, a)
            return(a)
        })
        a <- do.call(rbind, a)
        return(a)
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("base_count_in_peak_threshold", "base_count_in_motif_threshold", "maxFE_threshold", "FPR", "TPR")
    return(result)
}


plot_data_using_peak_base_count <- function(data_for_plot = NULL) {
    data_for_plot <- as.data.frame(data_for_plot)

    data_for_plot$base_count_in_motif_threshold <- factor(data_for_plot$base_count_in_motif_threshold,
        levels = sort(unique(data_for_plot$base_count_in_motif_threshold))
    )
    data_for_plot$base_count_in_peak_threshold <- factor(data_for_plot$base_count_in_peak_threshold,
        levels = sort(unique(data_for_plot$base_count_in_peak_threshold))
    )
    p <- ggplot(
        data_for_plot, aes(
            x = FPR, y = TPR,
            group = base_count_in_peak_threshold
        )
    ) +
        geom_line(aes(color = base_count_in_peak_threshold)) +
        geom_point(aes(shape = base_count_in_motif_threshold)) +
        # scale_x_continuous(limits = c(0, 1)) +
        # scale_y_continuous(limits = c(0, 1)) +
        # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") +
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

plot_data_using_peak_base_count_and_mean_FE <- function(data_for_plot = NULL) {
    data_for_plot <- as.data.frame(data_for_plot)

    data_for_plot$base_count_in_motif_threshold <- factor(data_for_plot$base_count_in_motif_threshold,
        levels = sort(unique(data_for_plot$base_count_in_motif_threshold))
    )
    data_for_plot$base_count_in_peak_threshold <- factor(data_for_plot$base_count_in_peak_threshold,
        levels = sort(unique(data_for_plot$base_count_in_peak_threshold))
    )
    data_for_plot$meanFE_threshold <- factor(data_for_plot$meanFE_threshold,
        levels = sort(unique(data_for_plot$meanFE_threshold))
    )
    p <- ggplot(
        data_for_plot, aes(
            x = FPR, y = TPR,
            group = base_count_in_peak_threshold
        )
    ) +
        geom_line(aes(color = base_count_in_peak_threshold)) +
        geom_point(aes(shape = base_count_in_motif_threshold, size = meanFE_threshold), alpha = 0.3) +
        # scale_x_continuous(limits = c(0, 1)) +
        # scale_y_continuous(limits = c(0, 1)) +
        # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") +
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

plot_data_using_peak_base_count_and_max_FE <- function(data_for_plot = NULL) {
    data_for_plot <- as.data.frame(data_for_plot)

    data_for_plot$base_count_in_motif_threshold <- factor(data_for_plot$base_count_in_motif_threshold,
        levels = sort(unique(data_for_plot$base_count_in_motif_threshold))
    )
    data_for_plot$base_count_in_peak_threshold <- factor(data_for_plot$base_count_in_peak_threshold,
        levels = sort(unique(data_for_plot$base_count_in_peak_threshold))
    )
    data_for_plot$maxFE_threshold <- factor(data_for_plot$maxFE_threshold,
        levels = sort(unique(data_for_plot$maxFE_threshold))
    )
    p <- ggplot(
        data_for_plot, aes(
            x = FPR, y = TPR,
            group = base_count_in_peak_threshold
        )
    ) +
        geom_line(aes(color = base_count_in_peak_threshold)) +
        geom_point(aes(shape = base_count_in_motif_threshold, size = maxFE_threshold), alpha = 0.3) +
        # scale_x_continuous(limits = c(0, 1)) +
        # scale_y_continuous(limits = c(0, 1)) +
        # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") +
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
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
path <- Args[6]

# 4. variable setting of test module--------------------------------------- TODO:
# path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/300_20"

# 5. process -------------------------------------------------------------- TODO:
setwd(path)
data <- fread("allChrom.window.filter.summary.gz", header = F, stringsAsFactors = F)
colnames(data) <- c("chrom", "start_minus_1", "end", "base_count_in_peak", "base_count_in_motif", "motif_count_in_window", "maxFE", "meanFE", "medianFE")

data_for_plot_using_peak_base_count <- compute_TPR_FPR_using_peak_base_count(data = data, base_count_in_peak_thresholds = c(50, 60, 70, 80, 90, 100) * 3, base_count_in_motif_thresholds = c(6, 9, 12, 15, 18, 21))
plot_data_using_peak_base_count(data_for_plot = data_for_plot_using_peak_base_count)
ggsave("fpr_tpr_using_peak_base_count.pdf", width = 8, height = 7)

data_for_plot_using_peak_base_count_and_mean_FE <- compute_TPR_FPR_using_peak_base_count_and_mean_FE(data = data, base_count_in_peak_thresholds = c(50, 60, 70, 80, 90, 100) * 3, base_count_in_motif_thresholds = c(6, 9, 12, 15, 18, 21), meanFE_thresholds = c(1.2, 1.5, 1.7, 1.8, 2))
plot_data_using_peak_base_count_and_mean_FE(data_for_plot = data_for_plot_using_peak_base_count_and_mean_FE)
ggsave("fpr_tpr_using_peak_base_count_and_mean_FE.pdf", width = 8, height = 7)

data_for_plot_using_peak_base_count_and_max_FE <- compute_TPR_FPR_using_peak_base_count_and_max_FE(data = data, base_count_in_peak_thresholds = c(50, 60, 70, 80, 90, 100) * 3, base_count_in_motif_thresholds = c(6, 9, 12, 15, 18, 21), maxFE_thresholds = c(1.2, 1.5, 1.7, 1.8, 2))
plot_data_using_peak_base_count_and_max_FE(data_for_plot = data_for_plot_using_peak_base_count_and_max_FE)
ggsave("fpr_tpr_using_peak_base_count_and_max_FE.pdf", width = 8, height = 7)
