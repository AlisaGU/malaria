#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:
read_sliding_window_data <- function(dataname = NULL) {
    data <- read.table(dataname, as.is = T, header = F)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    rownames(data1) <- paste("chrom", 1:14, sep = "")
    colnames(data1) <- colnames(data) <- c("nopeak", "nopeak500", paste("w", 1:10, sep = ""))
    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c("nopeak", "nopeak500", paste("w", 1:10, sep = "")))
    return(data_for_plot)
}


read_motif_flank_data <- function(dataname = NULL) {
    data <- read.table(dataname, as.is = T, header = F)
    data1 <- data[, seq(1, ncol(data), by = 2)] / data[, seq(2, ncol(data), by = 2)]
    rownames(data1) <- paste("chrom", 1:14, sep = "")
    colnames(data1) <- c(paste("site", 1:6, sep = ""), paste("flank", 1:10, sep = ""))
    data_for_plot <- data.frame(
        value = unlist(data1),
        chrom = rep(rownames(data1), times = ncol(data1)),
        source = rep(colnames(data1), each = nrow(data1))
    )
    data_for_plot$source <- factor(data_for_plot$source, levels = c(paste("site", 1:6, sep = ""), paste("flank", 1:10, sep = "")))
    data_for_plot$class <- NA
    data_for_plot$class[grep("site", data_for_plot$source)] <- "motif"
    data_for_plot$class[grep("flank", data_for_plot$source)] <- "motif_flank"
    return(data_for_plot)
}
# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
sliding_window <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/sliding_window/sliding_window_6mA_distri_FRACgt20_COVgt25", header = F, as.is = T)
motif_flank <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern/motif_flank_pattern_6mA_distri_FRACgt20_COVgt25", header = F, as.is = T)


sliding_window_frac <- sliding_window[, seq(1, ncol(sliding_window), by = 2)] / sliding_window[, seq(2, ncol(sliding_window), by = 2)]
rownames(sliding_window) <- 1:14

motif_flank_frac <- motif_flank[, seq(1, ncol(motif_flank), by = 2)] / motif_flank[, seq(2, ncol(motif_flank), by = 2)]
data <- rbind(sliding_window_frac, motif_flank_frac[, 1:3])
ggplot(data) +
    geom_boxplot(aes(x = source, y = value),
        position = position_dodge(), outlier.shape = NA
    ) +
    geom_point(aes(x = source, y = value), size = 2) +
    scale_x_discrete(labels = c("Control", "Control(500)", "Submit~5%", "5~10%", "10~15%", "15~20%", "20~25%", "25~30%", "30~35%", "35~40%", "40~45%", "45~50%", "Site1", "Site2", "Site3", "Site4", "Site5", "Site6", "Flank1-10")) +
    labs(y = "6mA sites density")+
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
            size = 30, color = "black",
            angle = 45, hjust = 0.5, vjust = 0.5
        ),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.title = element_blank(),
        legend.position = "None",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/6mA_frac_comparision_between_sliding_window_and_motif_flank.pdf", width = 20, height = 7)
