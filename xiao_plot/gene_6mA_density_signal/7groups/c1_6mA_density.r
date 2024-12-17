#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(ggbreak)))

# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(chip_info_name = NULL, input_info_name = NULL, chip_depth = NULL, input_depth = NULL) {
    chip_info <- fread(chip_info_name, stringsAsFactors = F, header = F)
    input_info <- fread(input_info_name, stringsAsFactors = F, header = F)
    colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

    chip_info <- chip_info %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")
    input_info <- input_info %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")

    chip_read_count <- chip_info %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))

    input_read_count <- input_info %>%
        group_by(gene_name) %>%
        summarise(depth = sum(depth))

    intersect_gene <- intersect(chip_read_count$gene_name, input_read_count$gene_name)
    chip_read_count_intersect <- chip_read_count[match(intersect_gene, chip_read_count$gene_name), ]
    input_read_count_intersect <- input_read_count[match(intersect_gene, input_read_count$gene_name), ]

    if (all(chip_read_count_intersect$gene_name == input_read_count_intersect$gene_name)) {
        density_signal <- (chip_read_count_intersect$depth / chip_depth) / (input_read_count_intersect$depth / input_depth)
        data <- data.frame(gene = chip_read_count_intersect$gene_name, density_signal = density_signal)
        colnames(data)[2] <- "gene_density_signal"
        return(data)
    }
}

my_test <- function(data_for_plot = NULL, group1 = NULL, group2 = NULL) {
    value1 <- data_for_plot$WT[data_for_plot$newGroup == group1]
    value2 <- data_for_plot$WT[data_for_plot$newGroup == group2]
    result <- list()
    result$wilcox <- wilcox.test(value1, value2)
    result$t <- t.test(value1, value2)
    return(result)
}

test_with_1 <- function(group = NULL) {
    value <- data_for_plot$WT[data_for_plot$newGroup == group]
    p <- wilcox.test(value, mu = 1, alternative = "greater")$p.value
    return(p)
}

my_test_median_oneside <- function(group2 = NULL) {
    value1 <- data_for_plot$WT[data_for_plot$newGroup == "RNA_translation"]
    value2 <- data_for_plot$WT[data_for_plot$newGroup == group2]
    result <- wilcox.test(value1, value2, alternative = "less")$p.value
    return(result)
}
# 3. input ---------------------------------------------------------------- TODO:
WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"

WT_chip_depth <- 36971023
WT_input_depth <- 21500156

setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/7groups")
gene_groups <- read.table("7Groups_list.txt", header = F, as.is = T)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(gene_groups) <- c("gene", "groups")
gene_groups <- gene_groups %>%
    filter(groups != "DNA_repair") %>%
    filter(groups != "merged")

gene_groups$newGroup <- sapply(gene_groups$groups, function(x) {
    return(ifelse(x == "eba" | x == "ama1" | x == "rh", "eba_ama1_rh", x))
})

gene_WT_6mA_density <- get_6mA(
    chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename,
    chip_depth = WT_chip_depth, input_depth = WT_input_depth
)
colnames(gene_WT_6mA_density)[2] <- "WT"


data_for_plot <- left_join(gene_groups, gene_WT_6mA_density, by = "gene")
# data_for_plot$newGroup <- factor(data_for_plot$newGroup, levels = c("RNA_translation", "var", "stevor", "rifin", "rhoptry", "msp", "eba_ama1_rh"))
data_for_plot$newGroup <- factor(data_for_plot$newGroup, levels = c("RNA_translation", "stevor", "rifin", "var", "rhoptry", "msp", "eba_ama1_rh"))

Horizontal_value <- median(data_for_plot$WT[data_for_plot$newGroup == "RNA_translation"])

ggplot(data_for_plot, aes(x = newGroup, y = WT)) +
    geom_boxplot(size = 1.5, fill = "white", outlier.shape = NA, color = "#2f2f2f") +
    geom_point(size = 5, color = "#939393", position = position_jitter(w = 0.2, h = 0)) +

    # geom_jitter(size = 5, color = "#d78825", width = 0.2, alpha = 0.5) +
    # geom_violin(size = 1.5, color = "black", fill = NA) +
    scale_x_discrete(labels = c("RNA translation", "stevor", "rifin", "var", "rhoptry", "msp", "eba, ama1, rh")) +
    scale_y_continuous(labels = c(1, 3, 5), breaks = c(1, 3, 5)) +
    geom_hline(yintercept = 1, size = 1, color = "#8d3727", linetype = "dashed") +
    # geom_hline(yintercept = Horizontal_value, size = 1, color = "#8d3727", linetype = "dashed") +
    labs(y = "6mA density in WT") +
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
# ggsave("7groups.pdf", width = 8, height = 8)
ggsave("7groups_changeOrder.pdf", width = 8, height = 8)

## compare between groups
# my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "eba_ama1")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "var")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "stevor")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "rifin")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "rhoptry")
my_test(data_for_plot = data_for_plot, group1 = "RNA_translation", group2 = "msp")

value1 <- data_for_plot$WT[data_for_plot$newGroup == "RNA_translation"]
value2 <- data_for_plot$WT[data_for_plot$newGroup == "eba_ama1" | data_for_plot$newGroup == "rh"]
result <- list()
result$wilcox <- wilcox.test(value1, value2, alternative = "less")$p.value
result$t <- t.test(value1, value2)$p.value
result

options(scipen = 10)
sapply(names(table(data_for_plot$newGroup)), my_test_median_oneside)
value1 <- data_for_plot$WT[data_for_plot$newGroup == "RNA_translation"]
value2 <- data_for_plot$WT[data_for_plot$newGroup == "eba_ama1" | data_for_plot$newGroup == "rh"]
wilcox.test(value1, value2, alternative = "less")$p.value


## compare with 1
sapply(names(table(data_for_plot$newGroup)), test_with_1)


## RNA_translation
ggplot(data_for_plot %>% filter(newGroup == "RNA_translation"), aes(x = newGroup, y = WT, label = gene)) +
    geom_point(size = 5, color = "red", alpha = 0.2) +
    geom_label_repel(aes(label = ifelse(WT > 3, as.character(gene), "")),
        box.padding = 1,
        point.padding = 0.5,
        segment.color = "grey50"
    )


RNA_tranlation_order <- data_for_plot %>%
    filter(newGroup == "RNA_translation") %>%
    arrange(desc(WT)) %>%
    select(c("gene", "WT"))
write.table(RNA_tranlation_order, "RNA_translation.ordered.txt", quote = F, sep = "\t", col.names = T, row.names = F)
