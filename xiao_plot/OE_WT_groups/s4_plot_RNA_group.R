#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))

library(ggplot2)
library(pheatmap)
# 2. functions ------------------------------------------------------------ TODO:

# 3. input ---------------------------------------------------------------- TODO:

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/")
load("RNA_class1/RNA_class1.Rdata")
class1_data <- data_for_plot
load("RNA_class5/RNA_class5.Rdata")
class5_data <- data_for_plot
data_for_plot <- rbind(class1_data, class5_data)
data_for_plot$type <- paste(data_for_plot$group, data_for_plot$prefix, sep = "_")

pdf("RNA_class1_5.OE.pdf", width = 10, height = 7)
ggplot(data_for_plot, aes(x = order, y = value, group = type, color = type)) +
    # geom_line() +
    # geom_point() +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, size = 3) +
    scale_color_manual(values = c(
        "OE_RNA_class5" = "#d96e6e", "WT_RNA_class5" = "#4d4747",
        "OE_RNA_class1" = "red", "WT_RNA_class1" = "black"
    ), guide = "none") +
    scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
    labs(y = "6mA signal density") +
    theme_classic() +
    theme(
        axis.ticks.length = unit(.25, "cm"),
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
        legend.text = element_text(size = 36, color = "black"),
        legend.key.width = unit(2, "cm"),
        plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
dev.off()
