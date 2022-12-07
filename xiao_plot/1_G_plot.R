#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function(prefix = NULL, array_1_label = NULL, array_2_label = NULL) {
    array_1_data <- read.csv(paste0(prefix, "_array_1.csv"), header = T, row.names = 1, as.is = T)
    array_1_data_mean <- apply(array_1_data, 2, mean)

    array_2_data <- read.csv(paste0(prefix, "_array_2.csv"), header = T, row.names = 1, as.is = T)
    array_2_data_mean <- apply(array_2_data, 2, mean)

    data_for_plot <- rbind(
        data.frame(value = array_1_data_mean, label = array_1_label, position = 1:length(array_1_data_mean)),
        data.frame(value = array_2_data_mean, label = array_2_label, position = 1:length(array_2_data_mean))
    )

    data_for_plot$label <- factor(data_for_plot$label, levels = c(array_1_label, array_2_label))
    return(data_for_plot)
}

plot_as_wen <- function(data_for_plot = NULL, array_color_value = NULL) {
    p <- ggplot(data = data_for_plot, aes(x = position, y = value)) +
        geom_line(aes(color = label), size = 4) +
        geom_area(aes(fill = label), alpha = 0.1, position = "identity") +
        scale_fill_manual(values = array_color_value, guide = "none") +
        scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
        scale_color_manual(values = array_color_value) +
        labs(y = "6mA signal density") +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = c(0.85, 0.9),
            legend.text = element_text(size = 36, color = "black"),
            legend.key.width = unit(4, "cm"),
            plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_bar_line <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = position, y = value, group = label)) +
        geom_col(aes(fill = label, color = label), position = "identity", alpha = 0.01) +
        geom_line(aes(color = label)) +
        scale_color_manual(values = c("WT" = "#a09999", "KD" = "red")) +
        scale_fill_manual(values = c("WT" = "white", "KD" = "red")) +
        scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
        labs(y = "6mA signal density") +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = c(0.85, 0.9),
            legend.text = element_text(size = 36, color = "black"),
            legend.key.width = unit(4, "cm"),
            plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_two_direction <- function(data_for_plot = NULL) {
    data_for_plot$value[data_for_plot$label == "WT"] <- -data_for_plot$value[data_for_plot$label == "WT"]
    ggplot(data = data_for_plot, aes(x = position, y = value)) +
        geom_line(aes(color = label), size = 4) +
        geom_area(aes(fill = label), position = "identity") +
        scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
        scale_y_continuous(labels = c(1.0, 0.5, 0.0, 0.5, 1.0), breaks = c(-1.0, -0.5, 0.0, 0.5, 1.0)) +
        labs(y = "6mA signal density") +
        coord_flip() +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = NULL,
            legend.text = element_text(size = 36, color = "black"),
            legend.key.width = unit(4, "cm"),
            plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
}

plot_point_and_line <- function(data_for_plot = NULL) {
    p <- ggplot(data_for_plot, aes(x = position, y = value, color = label, shape = label)) +
        geom_line() +
        geom_point(size = 3) +
        scale_color_manual(values = c("WT" = "#a09999", "KD" = "red")) +
        scale_fill_manual(values = c("WT" = "white", "KD" = "red")) +
        scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
        labs(y = "6mA signal density") +
        theme_classic() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            # panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(
                size = 36, color = "black"
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.position = c(0.85, 0.9),
            legend.text = element_text(size = 36, color = "black"),
            legend.key.width = unit(4, "cm"),
            plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot")
prefix <- "CT300-2-500-tss_gss"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read_data(prefix = prefix, array_1_label = "KD", array_2_label = "WT")
plot_as_wen(data_for_plot, array_color_value = c("KD" = "red", "WT" = "black"))
