#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggpubr)
library(patchwork)
library(ggridges)
# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(chip_info = NULL, input_info = NULL, chip_depth = NULL, input_depth = NULL) {
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
        return(data)
    }
}



group_gene <- function(data = NULL, number = NULL) {
    data_sort <- data[order(data$AV, decreasing = TRUE), ]
    data_sort$label <- ceiling((1:nrow(data_sort)) / (nrow(data_sort) / number))
    data_sort$label <- factor(data_sort$label, levels = number:1)
    return(data_sort)
}

group_gene_by_geneSMRT <- function(data = NULL, number = NULL) {
    data_sort <- data[order(data$gene_smrt_pro, decreasing = TRUE), ]
    data_sort$geneSMRTLabel <- ceiling((1:nrow(data_sort)) / (nrow(data_sort) / number))
    data_sort$geneSMRTLabel <- factor(data_sort$geneSMRTLabel, levels = number:1)
    return(data_sort)
}

group_gene_by_coluname <- function(data = NULL, colnumname = NULL, number = NULL) {
    data_sort <- data[order(data[, colnumname], decreasing = TRUE), ]
    data_sort[, paste0(colnumname, "_type")] <- ceiling((1:nrow(data_sort)) / (nrow(data_sort) / number))
    data_sort[, paste0(colnumname, "_type")] <- factor(data_sort[, paste0(colnumname, "_type")], levels = number:1)
    return(data_sort)
}

group_gene_by_coluname1 <- function(data = NULL, colnumname = NULL, number = NULL) {
    data_sort <- data[order(data[, colnumname], decreasing = TRUE), ]
    high_level_index <- 1:ceiling(nrow(data_sort) * 0.4)
    medium_level_index <- (ceiling(nrow(data_sort) * 0.4) + 1):(ceiling(nrow(data_sort) * 0.9))
    lower_level_index <- (ceiling(nrow(data_sort) * 0.9) + 1):nrow(data_sort)
    data_sort[high_level_index, paste0(colnumname, "_type")] <- 1
    data_sort[medium_level_index, paste0(colnumname, "_type")] <- 2
    data_sort[lower_level_index, paste0(colnumname, "_type")] <- 3

    data_sort[, paste0(colnumname, "_type")] <- factor(data_sort[, paste0(colnumname, "_type")], levels = 3:1)
    return(data_sort)
}

subsample_smrt <- function(expression_type = NULL, data = NULL, color = NULL) {
    observed <- as.data.frame(table(data$tss_smrt[data$RNA_expression_type == expression_type]))
    observed$Freq <- observed$Freq / sum(observed$Freq)
    observed$Var1 <- as.numeric(as.character(observed$Var1))
    observed$group <- "observed"

    gene_count <- length(which(data$RNA_expression_type == expression_type))

    simulation <- lapply(1:1000, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data$tss_smrt[index]
        data2 <- as.data.frame(table(data1))
        colnames(data2) <- c("Var1", "Freq")
        data2$Freq <- data2$Freq / sum(data2$Freq)
        data2$Var1 <- as.numeric(as.character(data2$Var1))
        data2$group <- paste0("simulation", x)
        return(data2)
    })
    simulation <- do.call(rbind, simulation)
    data_for_plot <- rbind(observed, simulation)



    p <- ggplot(data_for_plot, aes(x = Var1, y = Freq, group = group)) +
        geom_line(data = function(x) {
            x[x$group %in% paste("simulation", 1:100, sep = ""), ]
        }, color = "grey", alpha = 0.2, size = 2) +
        geom_line(data = function(x) {
            x[x$group %in% "observed", ]
        }, color = color, size = 2) +
        labs(y = "Fraction of genes", x = "Total starting 6mA in gene") +
        scale_y_continuous(limits = c(0, 0.3)) +
        scale_x_continuous(trans = "log10") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

subsample_smrt_density <- function(expression_type = NULL, data = NULL, color = NULL) {
    observed <- data.frame(value = data$tss_smrt_pro[data$RNA_expression_type == expression_type])
    observed$group <- "observed"

    gene_count <- length(which(data$RNA_expression_type == expression_type))

    simulation <- lapply(1:1000, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data.frame(value = data$tss_smrt_pro[index])

        data1$group <- paste0("simulation", x)
        return(data1)
    })
    simulation <- do.call(rbind, simulation)
    data_for_plot <- rbind(observed, simulation)



    # p <- ggplot(data_for_plot %>% filter(Var1 != 0), aes(x = Var1, y = Freq, group = group)) +
    p <- ggplot(data_for_plot, aes(x = value, group = group)) +
        geom_density(data = function(x) {
            x[x$group %in% paste("simulation", 1:100, sep = ""), ]
        }, color = "grey", alpha = 0.2, size = 2) +
        geom_density(data = function(x) {
            x[x$group %in% "observed", ]
        }, color = color, size = 2) +
        labs(y = "Density", x = "TSS smrt average density") +
        # scale_y_continuous(limits = c(0, 600)) +
        scale_x_continuous(trans = "log10") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

subsample_bam <- function(expression_type = NULL, data = NULL, color = NULL) {
    observed <- data.frame(value = data$tss_density_signal[data$RNA_expression_type == expression_type])
    observed$group <- "observed"

    gene_count <- length(which(data$RNA_expression_type == expression_type))

    simulation <- lapply(1:1000, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data.frame(value = data$tss_density_signal[index])

        data1$group <- paste0("simulation", x)
        return(data1)
    })
    simulation <- do.call(rbind, simulation)
    data_for_plot <- rbind(observed, simulation)



    # p <- ggplot(data_for_plot %>% filter(Var1 != 0), aes(x = Var1, y = Freq, group = group)) +
    p <- ggplot(data_for_plot, aes(x = value, group = group)) +
        geom_density(data = function(x) {
            x[x$group %in% paste("simulation", 1:100, sep = ""), ]
        }, color = "grey", alpha = 0.2, size = 2) +
        geom_density(data = function(x) {
            x[x$group %in% "observed", ]
        }, color = color, size = 2) +
        labs(y = "Density", x = "TSS bam average density") +
        scale_y_continuous(limits = c(0, 2.3)) +
        scale_x_continuous(trans = "log10") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

subsample_bam_gene_density <- function(expression_type = NULL, data = NULL, color = NULL) {
    observed <- data.frame(value = data$gene_flank2kb_density_signal[data$RNA_expression_type == expression_type])
    observed$group <- "observed"

    gene_count <- length(which(data$RNA_expression_type == expression_type))

    simulation <- lapply(1:100, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data.frame(value = data$gene_flank2kb_density_signal[index])

        data1$group <- paste0("simulation", x)
        return(data1)
    })
    simulation <- do.call(rbind, simulation)
    data_for_plot <- rbind(observed, simulation)



    # p <- ggplot(data_for_plot %>% filter(Var1 != 0), aes(x = Var1, y = Freq, group = group)) +
    p <- ggplot(data_for_plot, aes(x = value, group = group)) +
        geom_density(data = function(x) {
            x[x$group %in% paste("simulation", 1:100, sep = ""), ]
        }, color = "grey", alpha = 0.2, size = 2) +
        geom_density(data = function(x) {
            x[x$group %in% "observed", ]
        }, color = color, size = 2) +
        labs(y = "Density", x = "gene +-2kb bam average density") +
        # scale_y_continuous(limits = c(0, 2.3)) +
        scale_x_continuous(trans = "log10") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/gene_flank2kb.bed"



## T stage
chip_depth <- 36971023
input_depth <- 21500156
gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Tstage"
gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Tstage"
RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Tstage_new.txt" # tpm
# RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Tstage_fpkm.txt" # fpkm


## R stage
# chip_depth <- 54719200
# input_depth <- 50938409
# gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Rstage"
# gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Rstage"
# RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Rstage_new.txt"

## S stage
# chip_depth <- 43520326
# input_depth <- 26404329
# gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Sstage"
# gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Sstage"
# RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Sstage_new.txt"


upgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/up_gene"
downgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/down_gene"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
gene_chip_info <- fread(gene_chip_info_filename, stringsAsFactors = F, header = F)
gene_input_info <- fread(gene_input_info_filename, stringsAsFactors = F, header = F)

RNA <- fread(RNA_filename, stringsAsFactors = F, header = T) %>% na.omit()

up_gene <- as.character(unlist(read.table(upgene_filename, stringsAsFactors = F, header = F)))
down_gene <- as.character(unlist(read.table(downgene_filename, stringsAsFactors = F, header = F)))


colnames(all_genes_info) <- c("chrom", "start_minus_1", "end", "gene_name", "strand")
colnames(gene_chip_info) <- colnames(gene_input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "strand", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(RNA)[1] <- "gene"
## 整合数据
gene_density_signal_data <- get_6mA(
    chip_info = gene_chip_info, input_info = gene_input_info,
    chip_depth = chip_depth, input_depth = input_depth
)
colnames(gene_density_signal_data)[2] <- "gene_flank2kb_density_signal"
colnames(all_genes_info)[4] <- "gene"

data <- merge(merge(gene_density_signal_data, all_genes_info, by = "gene", all = TRUE), RNA, by = "gene", all = TRUE)

data <- data %>% filter(!is.na(chrom))
data <- data[, c(1, 2, 9)]

data$RNA_expression_type <- sapply(data$gene, function(x) {
    if (x %in% up_gene) {
        return("Up")
    } else if (x %in% down_gene) {
        return("Down")
    } else {
        return("Nodiff")
    }
})



## 上下调基因
up <- subsample_bam_gene_density(expression_type = "Up", data = data, color = rgb(240, 0, 0, maxColorValue = 255))
down <- subsample_bam_gene_density(expression_type = "Down", data = data, color = rgb(0, 0, 252, maxColorValue = 255))
up | down
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/up_down_gene+-2kb_density_distri.pdf", width = 14, height = 5)
## 所有基因 6mA density 和基因表达的关系
ggplot(data %>% filter(!is.infinite(gene_flank2kb_density_signal) & !is.na(gene_flank2kb_density_signal)), aes(x = gene_flank2kb_density_signal, y = AV)) +
    geom_point() +
    # scale_y_continuous(trans = "log10") +
    stat_smooth(aes(group = 1), method = "lm", col = "black", se = FALSE, show.legend = F) +
    stat_cor(aes(group = 1), method = "spearman", label.x.npc = 0.05, label.y.npc = 0.9, size = 12)
## 根据6mA density将基因分组
data1 <- group_gene_by_coluname(
    data = data %>% filter(!is.infinite(gene_flank2kb_density_signal) & !is.na(gene_flank2kb_density_signal)),
    colnumname = "gene_flank2kb_density_signal", number = 5
)
write.table(
    data1$gene[data1$gene_flank2kb_density_signal_type == "1"],
    "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/RNA_class1/RNA_class1.geneList.txt",
    row.names = F, col.names = F, quote = F
)
write.table(
    data1$gene[data1$gene_flank2kb_density_signal_type == "5"],
    "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/RNA_class5/RNA_class5.geneList.txt",
    row.names = F, col.names = F, quote = F
)
# data1 <- group_gene_by_coluname1(
#     data = data %>% filter(!is.infinite(gene_flank2kb_density_signal) & !is.na(gene_flank2kb_density_signal)),
#     colnumname = "gene_flank2kb_density_signal"
# )
# colnames(data1) <- c("gene", "6mA density in gene and +-2kb", "Average of TPM", "RNA expression type", "class of 6mA density")
# write.table(data1, "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/Tstage_gene_group_by_6mA_density.txt", sep = "\t", quote = F, row.names = F, col.names = T)
# wilcox.test(data1$AV[data1$gene_flank2kb_density_signal_type == "1"], data1$AV[data1$gene_flank2kb_density_signal_type == "2"], alternative = "greater")
ggplot(data1, aes(y = gene_flank2kb_density_signal_type, x = AV, fill = gene_flank2kb_density_signal_type)) +
    # stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
    # geom_density_ridges(alpha=0.5)+

    #     theme_ridges()+
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.5) +
    geom_text(
        data = data1 %>% group_by(gene_flank2kb_density_signal_type) %>%
            summarise(AV = median(AV)),
        aes(label = sprintf("%1.1f", AV)),
        position = position_nudge(y = -0.1, x = 20), colour = "black", size = 3.5
    ) +
    coord_cartesian(xlim = c(0, 500)) +
    scale_fill_manual(values = c("1" = "#8c2d04", "2" = "#cc4c02", "3" = "#ec7014", "4" = "#fe9929", "5" = "#fec44f"))



ggplot(data1, aes(x = gene_flank2kb_density_signal_type, y = AV, fill = gene_flank2kb_density_signal_type)) +
    scale_y_log10() +

    # geom_boxplot(position = position_dodge(), width = 0.5, outlier.shape = NA, size = 1) +
    # geom_point(position = position_jitter(w = 0.2, h = 0), alpha = 0.1) +
    geom_violin() +
    # stat_summary(
    #     fun = "mean",
    #     geom = "crossbar",
    #     width = 0.5,
    #     colour = "black"
    # ) +
    scale_fill_manual(values = c("1" = "#8c2d04", "2" = "#cc4c02", "3" = "#ec7014", "4" = "#fe9929", "5" = "#fec44f", "6" = "#fee391")) +
    labs(x = "Average 6mA density in gene+-2kb\n(chip/input)", y = "log10 RNA expression (TPM)") +
    # scale_x_discrete(labels = c("Low", "High"), breaks = c(3, 1)) +
    coord_cartesian(ylim = c(1, 500)) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(
            size = 30, color = "black",
        ),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 36, color = "black"),
        axis.title.x = element_text(size = 36, color = "black"),
        legend.title = element_blank(),
        legend.position = "None",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal//gene_bam_RNA_expression_Tstage_violin.pdf", width = 10, height = 10)
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal//gene_bam_RNA_expression_Rstage.pdf", width = 10, height = 10)
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal//gene_bam_RNA_expression_Sstage_violin.pdf", width = 10, height = 10)



summary_info <- data1 %>%
    group_by(gene_flank2kb_density_signal_type) %>%
    summarise(
        data = median(AV)
    )

summary_info$gene_flank2kb_density_signal_type <- as.factor(summary_info$gene_flank2kb_density_signal_type)

ggplot(summary_info, aes(x = gene_flank2kb_density_signal_type, y = data, color = gene_flank2kb_density_signal_type)) +
    # geom_line(aes(group = 1), size = 2) +
    geom_point(size = 7, show.legend = F) +
    stat_smooth(aes(group = 1), method = "lm", col = "black", se = FALSE, show.legend = F) +
    stat_cor(aes(group = 1), method = "spearman", label.x.npc = 0.05, label.y.npc = 0.9, size = 12) +
    scale_color_manual(values = c("1" = "#8c2d04", "2" = "#cc4c02", "3" = "#ec7014", "4" = "#fe9929", "5" = "#fec44f", "6" = "#fee391")) +
    labs(x = "Average 6mA density\nin gene+-2kb(chip/input)", y = "Median of RNA\nexpression (TPM)") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        plot.title = element_text(
            colour = "black",
            size = 40, vjust = 0.5, hjust = 0.5
        ),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(
            size = 30, color = "black",
        ),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 36, color = "black"),
        axis.title.x = element_text(size = 36, color = "black"),
        legend.title = element_blank(),
        legend.position = "None",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )

# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/gene_bam_RNA_expression_Tstage_trendline_tpm.pdf", width = 8, height = 8)
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/gene_bam_RNA_expression_Rstage_trendline_tpm.pdf", width = 8, height = 8)
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/gene_bam_RNA_expression_Sstage_trendline_tpm.pdf", width = 8, height = 8)
