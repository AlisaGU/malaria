#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
# 2. functions ------------------------------------------------------------ TODO:
my_test <- function(data_for_plot = NULL, group1 = NULL, group2 = NULL) {
    value1 <- data_for_plot$peak.pro[data_for_plot$newGroup == group1]
    value2 <- data_for_plot$peak.pro[data_for_plot$newGroup == group2]
    result <- list()
    result$wilcox <- wilcox.test(value1, value2)$p.value
    result$t <- t.test(value1, value2)$p.value
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
all_gene_6mA_proportion <- read.table("/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult1/mutation_density_distribution/wholeGenome/ESEA_WSEA_OCE_SAM_SAS/allgenes/3D7/peak_control_length_in_flank0kb_wholePeak")
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/7groups")
gene_groups <- read.table("7Groups_list.txt", header = F, as.is = T)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(all_gene_6mA_proportion) <- c("gene", "gene.peak.len", "gene.control.len")
colnames(gene_groups) <- c("gene", "groups")

all_gene_6mA_proportion$peak.pro <- apply(all_gene_6mA_proportion, 1, function(x) {
    as.numeric(unlist(x[2])) / (as.numeric(unlist(x[2])) + as.numeric(unlist(x[3])))
})

gene_groups <- gene_groups %>%
    # filter(groups != "DNA_repair") %>%
    filter(groups != "merged")

gene_groups$newGroup <- sapply(gene_groups$groups, function(x) {
    return(ifelse(x == "eba" | x == "ama1", "eba_ama1", x))
})


data_for_plot <- left_join(gene_groups, all_gene_6mA_proportion, by = "gene")
data_for_plot$newGroup <- factor(data_for_plot$newGroup, levels = c("RNA_translation", "DNA_repair", "var", "stevor", "rifin", "rhoptry", "msp", "eba_ama1", "rh"))

Horizontal_value <- median(data_for_plot$peak.pro[data_for_plot$newGroup == "RNA_translation"])

ggplot(data_for_plot, aes(x = newGroup, y = peak.pro)) +
    geom_boxplot(size = 1.5, fill = "white", outlier.shape = NA, color = "#2f2f2f") +
    geom_point(size = 5, color = "#939393", position = position_jitter(w = 0.2, h = 0)) +

    # geom_jitter(size = 5, color = "#d78825", width = 0.2, alpha = 0.5) +
    # geom_violin(size = 1.5, color = "black", fill = NA) +
    scale_x_discrete(labels = c("RNA translation", "DNA repair", "var", "stevor", "rifin", "rhoptry", "msp", "eba and ama1", "rh")) +
    geom_hline(yintercept = Horizontal_value, size = 1, color = "#8d3727", linetype = "dashed") +
    labs(y = "6mA proportion\nin gene body") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
            size = 36, color = "black", ,
            angle = 45, hjust = 1, vjust = 1
        ),
        axis.text.y = element_text(size = 36, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.position = "none"
    )
ggsave("9groups.6mAproportion.pdf", width = 8, height = 8)


# my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "eba_ama1")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "var")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "stevor")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "rifin")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "rhoptry")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "msp")


value1 <- data_for_plot$peak.pro[data_for_plot$newGroup == "RNA_translation"]
value2 <- data_for_plot$peak.pro[data_for_plot$newGroup == "eba_ama1" | data_for_plot$newGroup == "rh"]
result <- list()
result$wilcox <- wilcox.test(value1, value2)$p.value
result$t <- t.test(value1, value2)$p.value
result
