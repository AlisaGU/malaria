#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/up_down_WT_KD_T3")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
load("data_for_plot_KD.Rdata")
load("data_for_plot_wt.Rdata")

data_for_plot <- rbind(data_for_plot_wt, data_for_plot_KD)

ggplot(data_for_plot, aes(x = order, y = value, color = group, linetype = type)) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, size = 3) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue"), guide = "none") +
    scale_x_continuous(breaks = c(0, 100, 250, 400, 500), labels = c("2kb\nupstream", "Start", "gene body", "Stop", "2kb\ndownstream")) +
    scale_linetype_manual(values = c("KD" = "dashed", "WT" = "solid"), guide = "none") +
    # coord_cartesian(ylim = c(0.5, 2.6)) +
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
        legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 36, color = "black"),
        legend.key.width = unit(2, "cm"),
        plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
ggsave("combine.pdf", width = 10, height = 7)
