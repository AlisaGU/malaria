#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggpubr)
library(patchwork)
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

get_smrt_gene <- function(smrt_data = NULL) {
    smrt_data1 <- smrt_data %>%
        group_by(gene_name) %>%
        summarise(smrt_count = sum(depth))
    data <- merge(smrt_data1, all_genes_info, by = "gene_name")
    data$smrt_pro <- data$smrt_count / (data$end - data$start_minus_1)
    return(data[, c("gene_name", "smrt_pro")])
}

get_smrt_tss <- function(smrt_data = NULL) {
    smrt_data1 <- smrt_data %>%
        group_by(gene_name) %>%
        summarise(smrt_count = sum(depth))
    tss_info <- unique(smrt_data[, 1:4])
    data <- merge(smrt_data1, tss_info, by = "gene_name")
    data$smrt_pro <- data$smrt_count / (data$end - data$start_minus_1)
    return(data[, c("gene_name", "smrt_pro")])
}

get_smrt1 <- function(smrt_data = NULL) {
    smrt_data1 <- smrt_data %>%
        group_by(gene_name) %>%
        summarise(smrt_count = sum(depth))
    data <- merge(smrt_data1, all_genes_info, by = "gene_name")
    return(data[, c("gene_name", "smrt_count")])
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
    medium_level_index <- (ceiling(nrow(data_sort) * 0.4) + 1):(ceiling(nrow(data_sort) * 0.7))
    lower_level_index <- (ceiling(nrow(data_sort) * 0.7) + 1):nrow(data_sort)
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
# 3. input ---------------------------------------------------------------- TODO:
all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/all_genes_bed"

gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/R_stage/chip.inter"
gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/R_stage/input.inter"


chip_depth <- 54719200
input_depth <- 50938409
RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Rstage_new.txt"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
gene_chip_info <- fread(gene_chip_info_filename, stringsAsFactors = F, header = F)
gene_input_info <- fread(gene_input_info_filename, stringsAsFactors = F, header = F)


RNA <- fread(RNA_filename, stringsAsFactors = F, header = T) %>% na.omit()


##
colnames(all_genes_info) <- c("chrom", "start_minus_1", "end", "gene_name")
colnames(gene_chip_info) <- colnames(gene_input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(RNA)[1] <- "gene"
## density 和 smrt
gene_density_signal_data <- get_6mA(
    chip_info = gene_chip_info, input_info = gene_input_info,
    chip_depth = chip_depth, input_depth = input_depth
)
colnames(gene_density_signal_data)[2] <- "gene_density_signal"


data <- merge(merge(gene_density_signal_data, RNA, by = "gene", all = TRUE), all_genes_info, by.x = "gene", by.y = "gene_name")

data <- data %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")
data <- data[, c(1, 2, 5, 6)]

## 三分组
data$RNA_classIV <- sapply(data$AV, function(x) {
    if (x <= 10) {
        return("low")
    } else if (x > 10 & x <= 1000) {
        return("medium")
    } else if (x > 1000) {
        return("high")
    }
})
data$RNA_classIV <- factor(data$RNA_classIV, levels = c("low", "medium", "high"))

ggplot(data, aes(x = RNA_classIV, y = gene_density_signal, color = RNA_classIV)) +
    geom_violin(position = position_dodge(), size = 1) +
    geom_boxplot(position = position_dodge(), width = 0.5) +
    geom_point(position = position_jitterdodge()) +
    scale_y_continuous(trans = "log10") +
    scale_x_discrete(labels = c("Low", "Medium", "High")) +
    scale_color_manual(values = c("#9e9ac8", "#756bb1", "#54278f")) +
    labs(y = "6mA density", x = "RNA expression") +
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
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/R_stage/RNA_6mA_IV.pdf", width = 7, height = 7)

## 多分组
# data1 <- group_gene_by_coluname(
#     data = data %>% filter(!is.infinite(gene_density_signal) & !is.na(gene_density_signal)),
#     colnumname = "gene_density_signal", number = 4
# )
data1 <- group_gene_by_coluname1(
    data = data %>% filter(!is.infinite(gene_density_signal) & !is.na(gene_density_signal)),
    colnumname = "gene_density_signal"
)
ggplot(data1, aes(x = gene_density_signal_type, y = AV, fill = gene_density_signal_type)) +
    geom_boxplot(position = position_dodge(), width = 0.5, outlier.shape = NA, size = 1) +
    coord_cartesian(ylim = c(0, 400)) +
    scale_fill_manual(values = c("1" = "#8c2d04", "2" = "#cc4c02", "3" = "#ec7014", "4" = "#fe9929", "5" = "#fec44f", "6" = "#fee391")) +
    labs(x = "Average 6mA density in gene\n(chip/input)", y = "RNA expression (TPM)") +
    scale_x_discrete(labels = c("Low", "High"), breaks = c(3, 1)) +
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
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/R_stage/gene_bam_RNA_expression_Rstage_new.pdf", width = 10, height = 10)
