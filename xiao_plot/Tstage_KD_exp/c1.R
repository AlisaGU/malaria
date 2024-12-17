#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(chip_info_name = NULL, input_info_name = NULL, chip_depth = NULL, input_depth = NULL) {
    chip_info <- fread(chip_info_name, stringsAsFactors = F, header = F)
    input_info <- fread(input_info_name, stringsAsFactors = F, header = F)
    colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "strand", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

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
        colnames(data)[2] <- "gene_flank2kb_density_signal"
        return(data)
    }
}


get_WT_KD_6mA <- function() {
    WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Tstage"
    WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Tstage"
    WT_chip_depth <- 36971023
    WT_input_depth <- 21500156

    KD_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/chip.gene_flank2kb.inter"
    KD_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/input.gene_flank2kb.inter"
    KD_chip_depth <- 39189761
    KD_input_depth <- 27773718

    gene_WT_6mA_density <- get_6mA(chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename, chip_depth = WT_chip_depth, input_depth = WT_input_depth)
    gene_KD_6mA_density <- get_6mA(chip_info_name = KD_chip_info_filename, input_info_name = KD_input_info_filename, chip_depth = KD_chip_depth, input_depth = KD_input_depth)
    colnames(gene_WT_6mA_density)[2] <- "WT"
    colnames(gene_KD_6mA_density)[2] <- "KD"

    result <- dplyr::left_join(gene_WT_6mA_density, gene_KD_6mA_density, by = "gene") %>% na.omit()
    colnames(result) <- c("gene", "methy_WT", "methy_KD")
    result$methy_KDsubtractWT <- result$methy_KD - result$methy_WT
    result$methy_KDdivideWT <- result$methy_KD / result$methy_WT

    WT_gt_KD_index <- which(result$methy_WT > result$methy_KD)
    result$methy_diff <- NA
    result$methy_diff[WT_gt_KD_index] <- result$methy_WT[WT_gt_KD_index] / result$methy_KD[WT_gt_KD_index]
    result$methy_diff[-WT_gt_KD_index] <- result$methy_KD[-WT_gt_KD_index] / result$methy_WT[-WT_gt_KD_index]
    return(result)
}


get_WT_KD_exp <- function() {
    exp_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/Tstage_KD_exp/T_6mAKDvs3D7_dotplotdata.txt"
    exp <- fread(exp_filename, header = T, stringsAsFactors = F, sep = "\t")
    exp1 <- data.frame(
        gene = exp$gene_id,
        WT = apply(exp[, c("TPM_3D7_T2", "TPM_3D7_T10")], 1, mean),
        KD = apply(exp[, c("TPM_6mAKD_T2", "TPM_6mAKD_T8")], 1, mean),
        significant = exp$significant
    )

    colnames(exp1) <- c("gene", "exp_WT", "exp_KD", "exp_significant")
    exp1$exp_KDsubtractWT <- exp1$exp_KD - exp1$exp_WT
    exp1$exp_KDdivideWT <- exp1$exp_KD / exp1$exp_WT

    WT_gt_KD_index <- which(exp1$exp_WT > exp1$exp_KD)
    exp1$exp_diff <- NA
    exp1$exp_diff[WT_gt_KD_index] <- exp1$exp_WT[WT_gt_KD_index] / exp1$exp_KD[WT_gt_KD_index]
    exp1$exp_diff[-WT_gt_KD_index] <- exp1$exp_KD[-WT_gt_KD_index] / exp1$exp_WT[-WT_gt_KD_index]

    return(exp1)
}
# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/Tstage_KD_exp")

# 5. process -------------------------------------------------------------- TODO:
methy <- get_WT_KD_6mA()
exp <- get_WT_KD_exp()
data_for_plot <- dplyr::left_join(methy, exp, by = "gene")


#####
data_for_plot <- data_for_plot %>% arrange(methy_WT)
data_for_plot$methy_WT_level <- ceiling(1:nrow(data_for_plot) / (nrow(data_for_plot) / 5)) # level 1的6mA最低
data_for_plot$methy_WT_level <- as.factor(data_for_plot$methy_WT_level)



# down组
down <- data_for_plot$gene[data_for_plot$exp_significant == "down"]
sapply(1:5, function(x) {
    chisq_data <- table(data_for_plot$exp_significant == "down", data_for_plot$methy_WT_level == x)
    chisq.test(chisq_data)
})

ggplot(data = data_for_plot, aes(x = methy_WT_level, y = -methy_KDsubtractWT)) +
    geom_boxplot(size = 1, outlier.shape = NA, color = "#edb126") +
    geom_point(data = data_for_plot %>% filter(exp_significant == "down"), aes(x = methy_WT_level, y = -methy_KDsubtractWT, ), color = "#67aec7", size = 1, position = position_jitter(w = 0.2, h = 0)) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-0.5, 0.6)) +
    annotate("text", x = 1, y = 0.6, label = "25") +
    annotate("text", x = 2, y = 0.6, label = "52") +
    annotate("text", x = 3, y = 0.6, label = "60") +
    annotate("text", x = 4, y = 0.6, label = "100") +
    annotate("text", x = 5, y = 0.6, label = "117") +
    labs(y = "Difference of total 6mA\n(WT-KD)", x = "Gene groups") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none"
    )
ggsave("WTsubtractKD.pdf", height = 4, width = 3.5)



# up组
up <- data_for_plot$gene[data_for_plot$exp_significant == "up"]

sapply(1:5, function(x) {
    chisq_data <- table(data_for_plot$exp_significant == "up", data_for_plot$methy_WT_level == x)
    chisq.test(chisq_data)
})



ggplot(data = data_for_plot, aes(x = methy_WT_level, y = -methy_KDsubtractWT)) +
    geom_boxplot(size = 1, outlier.shape = NA, color = "#edb126") +
    geom_point(data = data_for_plot %>% filter(exp_significant == "up"), aes(x = methy_WT_level, y = -methy_KDsubtractWT, ), color = "#D6764F", size = 1, position = position_jitter(w = 0.2, h = 0)) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-0.5, 0.6)) +
    annotate("text", x = 1, y = 0.6, label = "48") +
    annotate("text", x = 2, y = 0.6, label = "34") +
    annotate("text", x = 3, y = 0.6, label = "19") +
    annotate("text", x = 4, y = 0.6, label = "19") +
    annotate("text", x = 5, y = 0.6, label = "22") +
    labs(y = "Difference of total 6mA\n(WT-KD)", x = "Gene groups") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none"
    )
ggsave("WTsubtractKD-up.pdf", height = 4, width = 3.5)


# up和down
ggplot(data = data_for_plot, aes(x = methy_WT_level, y = -methy_KDsubtractWT)) +
    geom_boxplot(size = 1, outlier.shape = NA, color = "#edb126") +
    geom_point(data = data_for_plot %>% filter(exp_significant == "up"), aes(x = methy_WT_level, y = -methy_KDsubtractWT, ), color = "#d6764f", size = 1, position = position_jitter(w = 0.2, h = 0)) +
    geom_point(data = data_for_plot %>% filter(exp_significant == "down"), aes(x = methy_WT_level, y = -methy_KDsubtractWT, ), color = "#67aec7", size = 1, position = position_jitter(w = 0.2, h = 0)) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-0.5, 0.6)) +
    annotate("text", x = 0.64, y = 0.62, label = "Up:", color = "#d6764f") +
    annotate("text", x = 0.64, y = 0.56, label = "Down:", color = "#67aec7") +
    annotate("text", x = 1, y = 0.62, label = "48", color = "#d6764f") +
    annotate("text", x = 2, y = 0.62, label = "34", color = "#d6764f") +
    annotate("text", x = 3, y = 0.62, label = "19", color = "#d6764f") +
    annotate("text", x = 4, y = 0.62, label = "19", color = "#d6764f") +
    annotate("text", x = 5, y = 0.62, label = "22", color = "#d6764f") +
    annotate("text", x = 1, y = 0.56, label = "25", color = "#67aec7") +
    annotate("text", x = 2, y = 0.56, label = "52", color = "#67aec7") +
    annotate("text", x = 3, y = 0.56, label = "60", color = "#67aec7") +
    annotate("text", x = 4, y = 0.56, label = "100", color = "#67aec7") +
    annotate("text", x = 5, y = 0.56, label = "117", color = "#67aec7") +
    labs(y = "Difference of total 6mA\n(WT-KD)", x = "Gene groups") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none"
    )
ggsave("WTsubtractKD-upAndDown.pdf", height = 4, width = 6)
#####
data_filtered <- data_for_plot %>%
    filter(exp_significant != "LowTPM") %>%
    filter(exp_KDdivideWT != 0) %>%
    filter(!is.infinite(exp_KDdivideWT))


ggplot(data_filtered %>% filter(exp_significant == "stable"), aes(x = methy_diff, y = exp_diff, color = exp_significant, group = 1)) +
    geom_point() +
    # scale_x_continuous(trans = "log2", limits=c(max(data_filtered$methy_diff), min(data_filtered$methy_diff))) +
    scale_x_continuous(trans = "reverse") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = c("down" = "blue", "up" = "red", "LowTPM" = "grey", "stable" = "grey")) +
    geom_smooth(method = "lm", se = T) +
    labs(x = "methy_diff", y = "exp_diff") +
    stat_cor(method = "pearson") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.position = "none"
    )


#####
ggplot(data_for_plot, aes(x = methy_WT, y = exp_WT)) +
    geom_point() +
    stat_cor(method = "spearman") +
    geom_smooth(method = "lm", se = T)
