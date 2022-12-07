#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(data.table)

# 2. functions ------------------------------------------------------------


# 3. input ----------------------------------------------------------------
method1_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA/sliding_window/relative0.05/mutation_density_by_whole_chromosome_with_submit"
method2_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA/sliding_window_noCollapseWindow/relative0.05/mutation_density"
method1_all <- read.table(paste0(method1_dir, "/allvar_mean_variant_count"), header = F, stringsAsFactors = F)
method2_all <- read.table(paste0(method2_dir, "/allvar_mean_variant_count"), header = F, stringsAsFactors = F)

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
method1_chrom_mean <- apply(method1_all[, 2:11], 2, mean)
method2_chrom_mean <- apply(method2_all[, 2:11], 2, mean)
wilcox.test(method1_chrom_mean, method2_chrom_mean, paired = T)
mean_data <- data.frame(
    value = c(method1_chrom_mean, method2_chrom_mean),
    method = c(rep("method1", 10), rep("method2", 10))
)

mean_p <- ggboxplot(data = mean_data, x = "method", y = "value", color = "method", add = "jitter", add.params = list(size = 3)) +
    stat_compare_means(aes(group = method),
        label = "p.signif", paired = T, size = 10
    ) +
    labs(y = "", x = "", title = "mean") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(size = 24, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 28, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 28, color = "black"),
        legend.position = "None",
        plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
    )


method1_chrom_var <- apply(method1_all[, 2:11], 2, var)
method2_chrom_var <- apply(method2_all[, 2:11], 2, var)
wilcox.test(method1_chrom_var, method1_chrom_var, paired = T)
var_data <- data.frame(
    value = c(method1_chrom_var, method2_chrom_var),
    method = c(rep("method1", 10), rep("method2", 10)),
    source = c(paste("w", 1:10, sep = ""), paste("w", 1:10, sep = ""))+
    labs(y = "", x = "", title = "var") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(size = 24, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 28, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 28, color = "black"),
        legend.position = "None",
        plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
    )
)
ggplot(data = var_data,aes(x=method,y=value)) +
    geom_boxplot()+
    geom_point(aes(color=source),size=3)+
    labs(y = "", x = "", title = "var") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(size = 24, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 28, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 28, color = "black"),
        plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
    )

var_p <- ggboxplot(data = var_data, x = "method", y = "value", color = "method", add = "jitter", add.params = list(size = 3)) +
    labs(y = "", x = "", title = "var") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(size = 24, color = "black"),
        plot.title = element_text(
            colour = "black",
            size = 28, vjust = 0.5, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 28, color = "black"),
        legend.position = "None",
        plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
    )
mean_p | var_p
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA/compare_two_mutation_density_method/compare.pdf", width = 14, height = 7)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA/compare_two_mutation_density_method/compare.jpg", width = 14, height = 7)