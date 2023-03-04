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

subsample_bam_gene_density <- function(expression_type = NULL, data = NULL, color = NULL) {
    observed <- data.frame(value = data$gene_density_signal[data$RNA_expression_type == expression_type])
    observed$group <- "observed"

    gene_count <- length(which(data$RNA_expression_type == expression_type))

    simulation <- lapply(1:100, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data.frame(value = data$gene_density_signal[index])

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
        labs(y = "Density", x = "gene bam average density") +
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
all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/all_genes_bed"

gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"
tss_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/tss.chip.inter"
tss_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/tss.input.inter"

gene_smrt_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/gene.smrt.bed"
tss_smrt_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/tss.smrt.bed"

chip_depth <- 36971023
input_depth <- 21500156
# RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Tstage.txt"
RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Tstage_new.txt"
TSS_peakpro_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/TSS_peakpro.txt"
upgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/up_gene"
downgene_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_WT210vsKD28/down_gene"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
gene_chip_info <- fread(gene_chip_info_filename, stringsAsFactors = F, header = F)
gene_input_info <- fread(gene_input_info_filename, stringsAsFactors = F, header = F)
tss_chip_info <- fread(tss_chip_info_filename, stringsAsFactors = F, header = F, drop = c(5, 6))
tss_input_info <- fread(tss_input_info_filename, stringsAsFactors = F, header = F, drop = c(5, 6))

gene_smrt <- fread(gene_smrt_filename, stringsAsFactors = F, header = F)
tss_smrt <- fread(tss_smrt_filename, stringsAsFactors = F, header = F, drop = c(5, 6))

RNA <- fread(RNA_filename, stringsAsFactors = F, header = T) %>% na.omit()
TSS_peakpro <- fread(TSS_peakpro_filename, stringsAsFactors = FALSE, header = TRUE)

up_gene <- as.character(unlist(read.table(upgene_filename, stringsAsFactors = F, header = F)))
down_gene <- as.character(unlist(read.table(downgene_filename, stringsAsFactors = F, header = F)))

## process TSS_peakpro
TSS_peakpro <- TSS_peakpro %>%
    group_by(chrom, TSS_start, TSS_end, gene_name, strand) %>%
    summarise(overlap_base = sum(overlap_base))

TSS_peakpro$pro <- TSS_peakpro$overlap_base / (TSS_peakpro$TSS_end - TSS_peakpro$TSS_start)
##
colnames(all_genes_info) <- c("chrom", "start_minus_1", "end", "gene_name")
colnames(gene_chip_info) <- colnames(gene_input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(tss_chip_info) <- colnames(tss_input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
colnames(gene_smrt) <- colnames(tss_smrt) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

colnames(RNA)[1] <- "gene"
## density 和 smrt
gene_density_signal_data <- get_6mA(
    chip_info = gene_chip_info, input_info = gene_input_info,
    chip_depth = chip_depth, input_depth = input_depth
)
colnames(gene_density_signal_data)[2] <- "gene_density_signal"


tss_density_signal_data <- get_6mA(
    chip_info = tss_chip_info, input_info = tss_input_info,
    chip_depth = chip_depth, input_depth = input_depth
)
colnames(tss_density_signal_data)[2] <- "tss_density_signal"

gene_smrtpro <- get_smrt_gene(smrt_data = gene_smrt)
colnames(gene_smrtpro)[2] <- "gene_smrt_pro"


tss_smrtpro <- get_smrt_tss(smrt_data = tss_smrt)
colnames(tss_smrtpro)[2] <- "tss_smrt_pro"
tss_smrtpro1 <- get_smrt1(smrt_data = tss_smrt)
colnames(tss_smrtpro1)[2] <- "tss_smrt"

data <- merge(merge(merge(merge(merge(merge(gene_density_signal_data, tss_density_signal_data, by = "gene", all = TRUE),
    gene_smrtpro,
    by.x = "gene", by.y = "gene_name", all = TRUE
),
tss_smrtpro,
by.x = "gene", by.y = "gene_name", all = TRUE
),
tss_smrtpro1,
by.x = "gene", by.y = "gene_name", all = TRUE
),
TSS_peakpro,
by.x = "gene", by.y = "gene_name", all = TRUE
),
RNA,
by = "gene", all = TRUE
)

data <- data %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")
data <- data[, c(1:6, 12, 15)]
colnames(data)[7] <- "peak_in_TSS_pro"

data$RNA_expression_type <- sapply(data$gene, function(x) {
    if (x %in% up_gene) {
        return("Up")
    } else if (x %in% down_gene) {
        return("Down")
    } else {
        return("Nodiff")
    }
})

## classI:根据表达量分类基因：very low (<0.5), LOW (0.5~10), medium (11~1000), high (>1000)
## 类别：    high      low   medium very low
## 基因数：  3455       243     1712      192
data$RNA_classI <- sapply(data$AV, function(x) {
    if (x <= 0.5) {
        return("very low")
    } else if (x > 0.5 & x <= 10) {
        return("low")
    } else if (x > 10 & x <= 1000) {
        return("medium")
    } else if (x > 1000) {
        return("high")
    }
})
data$RNA_classI <- factor(data$RNA_classI, levels = c("very low", "low", "medium", "high"))

data %>%
    group_by(RNA_classI) %>%
    summarise(median(gene_density_signal, na.rm = TRUE))
ggplot(data, aes(x = RNA_classI, y = gene_density_signal, color = RNA_classI)) +
    geom_violin(position = position_dodge(), size = 1) +
    geom_boxplot(position = position_dodge(), width = 0.5) +
    geom_point(position = position_jitterdodge()) +
    scale_y_continuous(trans = "log10")


## classII: 根据表达量分类基因：LOW (0~10), medium (11~200), high (>200)
## 类别：    high      low   medium
## 基因数：  4528       435     639
data$RNA_classII <- sapply(data$AV, function(x) {
    if (x <= 10) {
        return("low")
    } else if (x > 10 & x <= 200) {
        return("medium")
    } else if (x > 200) {
        return("high")
    }
})
data$RNA_classII <- factor(data$RNA_classII, levels = c("low", "medium", "high"))
data %>%
    group_by(RNA_classII) %>%
    summarise(median(tss_density_signal, na.rm = TRUE))
ggplot(data, aes(x = RNA_classII, y = tss_density_signal, color = RNA_classII)) +
    geom_violin(position = position_dodge(), size = 1) +
    geom_boxplot(position = position_dodge(), width = 0.5) +
    geom_point(position = position_jitterdodge()) +
    scale_y_continuous(trans = "log10")


low_value <- data$tss_density_signal[data$RNA_classII == "low"]
low_value <- low_value[!is.infinite(low_value) & !is.na(low_value)]

medium_value <- data$tss_density_signal[data$RNA_classII == "medium"]
medium_value <- medium_value[!is.infinite(medium_value) & !is.na(medium_value)]

high_value <- data$tss_density_signal[data$RNA_classII == "high"]
high_value <- high_value[!is.infinite(high_value) & !is.na(high_value)]

t.test(low_value, medium_value, alternative = "less")
t.test(low_value, high_value, alternative = "less")
t.test(medium_value, high_value, alternative = "less")

## classIII: 根据表达量分类基因：very low (<0.2), LOW (0.2~10), medium (11~100), high (>100)
## 类别：very low      low   medium     high
## 基因数：     162      273     393     4774
data$RNA_classIII <- sapply(data$AV, function(x) {
    if (x <= 0.2) {
        return("very low")
    } else if (x > 0.2 & x <= 10) {
        return("low")
    } else if (x > 10 & x <= 100) {
        return("medium")
    } else if (x > 100) {
        return("high")
    }
})
data$RNA_classIII <- factor(data$RNA_classIII, levels = c("very low", "low", "medium", "high"))
data %>%
    group_by(RNA_classIII) %>%
    summarise(median(tss_density_signal, na.rm = TRUE))
ggplot(data, aes(x = RNA_classIII, y = tss_density_signal, color = RNA_classIII)) +
    geom_violin(position = position_dodge(), size = 1) +
    geom_boxplot(position = position_dodge(), width = 0.5) +
    geom_point(position = position_jitterdodge()) +
    scale_y_continuous(trans = "log10")

verylow_value <- data$tss_density_signal[data$RNA_classIII == "very low"]
verylow_value <- verylow_value[!is.infinite(verylow_value) & !is.na(verylow_value)]

low_value <- data$tss_density_signal[data$RNA_classIII == "low"]
low_value <- low_value[!is.infinite(low_value) & !is.na(low_value)]

medium_value <- data$tss_density_signal[data$RNA_classIII == "medium"]
medium_value <- medium_value[!is.infinite(medium_value) & !is.na(medium_value)]

high_value <- data$tss_density_signal[data$RNA_classIII == "high"]
high_value <- high_value[!is.infinite(high_value) & !is.na(high_value)]

t.test(verylow_value, low_value, alternative = "less")

t.test(low_value, medium_value, alternative = "less")
t.test(low_value, high_value, alternative = "less")
t.test(medium_value, high_value, alternative = "less")
## classIV:根据表达量分类基因： LOW (<10), medium (11~1000), high (>1000)
## 类别：    high      low   medium
## 基因数：  3455       435     1712
## 选这一版！！！
# TODO:
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

data %>%
    group_by(RNA_classIV) %>%
    summarise(median(gene_density_signal, na.rm = TRUE))
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
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/RNA_6mA_IV.pdf", width = 7, height = 7)

low_value <- data$gene_density_signal[data$RNA_classIV == "low"]
low_value <- low_value[!is.infinite(low_value) & !is.na(low_value)]

medium_value <- data$gene_density_signal[data$RNA_classIV == "medium"]
medium_value <- medium_value[!is.infinite(medium_value) & !is.na(medium_value)]

high_value <- data$gene_density_signal[data$RNA_classIV == "high"]
high_value <- high_value[!is.infinite(high_value) & !is.na(high_value)]

wilcox.test(low_value, medium_value, alternative = "less") # 0.02708
wilcox.test(low_value, high_value, alternative = "less") # p-value < 2.2e-16
wilcox.test(medium_value, high_value, alternative = "less") # p-value < 2.2e-16

## 上下调基因
up <- subsample_smrt(expression_type = "Up", data = data, color = rgb(240, 0, 0, maxColorValue = 255))
down <- subsample_smrt(expression_type = "Down", data = data, color = rgb(0, 0, 252, maxColorValue = 255)) + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
)
up | down
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/SMRT_starting.pdf", width = 14, height = 5)


up <- subsample_bam(expression_type = "Up", data = data, color = rgb(240, 0, 0, maxColorValue = 255))
down <- subsample_bam(expression_type = "Down", data = data, color = rgb(0, 0, 252, maxColorValue = 255)) + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
)
up | down
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/bam_starting.pdf", width = 14, height = 5)


up <- subsample_smrt_density(expression_type = "Up", data = data, color = rgb(240, 0, 0, maxColorValue = 255))
down <- subsample_smrt_density(expression_type = "Down", data = data, color = rgb(0, 0, 252, maxColorValue = 255)) + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
)
up | down
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/SMRT_density_starting.pdf", width = 14, height = 5)


up <- subsample_bam_gene_density(expression_type = "Up", data = data, color = rgb(240, 0, 0, maxColorValue = 255))
down <- subsample_bam_gene_density(expression_type = "Down", data = data, color = rgb(0, 0, 252, maxColorValue = 255))
up | down
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/up_down_gene_density_distri.pdf", width = 14, height = 5)
## 基因SMRT和TSS SMRT
# data1 <- group_gene_by_coluname(
#     data = data %>% filter(!is.infinite(gene_density_signal) & !is.na(gene_density_signal)),
#     colnumname = "gene_density_signal", number = 3
# )
data1 <- group_gene_by_coluname1(
    data = data %>% filter(!is.infinite(gene_density_signal) & !is.na(gene_density_signal)),
    colnumname = "gene_density_signal"
)

ggplot(data1, aes(x = gene_density_signal_type, y = AV, fill = gene_density_signal_type)) +
    geom_boxplot(position = position_dodge(), width = 0.5, outlier.shape = NA, size = 1) +
    # geom_boxplot(position = position_dodge(), width = 0.5, size = 1) +
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
# ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/gene_bam_RNA_expression.pdf", width = 10, height = 10)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/gene_bam_RNA_expression_Tstage_new.pdf", width = 10, height = 10)


data2 <- group_gene_by_coluname(
    data = data %>% filter(!is.infinite(gene_smrt_pro) & !is.na(gene_smrt_pro)),
    colnumname = "gene_smrt_pro", number = 6
)
ggplot(data2, aes(x = gene_smrt_pro_type, y = AV, color = gene_smrt_pro_type)) +
    geom_boxplot(position = position_dodge(), width = 0.2, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 15000))


data3 <- group_gene_by_coluname(
    data = data %>% filter(!is.infinite(tss_density_signal) & !is.na(tss_density_signal)),
    colnumname = "tss_density_signal", number = 10
)
ggplot(data3, aes(x = tss_density_signal_type, y = AV, color = tss_density_signal_type)) +
    geom_boxplot(position = position_dodge(), width = 0.2) +
    coord_cartesian(ylim = c(0, 500))

data4 <- group_gene_by_coluname(
    data = data %>% filter(!is.infinite(tss_smrt_pro) & !is.na(tss_smrt_pro)),
    colnumname = "tss_smrt_pro", number = 6
)
ggplot(data4, aes(x = tss_smrt_pro_type, y = AV, color = tss_smrt_pro_type)) +
    geom_boxplot(position = position_dodge(), width = 0.2, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 15000))



data5 <- group_gene_by_coluname(
    data = data %>% filter(!is.infinite(tss_smrt) & !is.na(tss_smrt)),
    colnumname = "tss_smrt", number = 4
)
ggplot(data5, aes(x = tss_smrt_type, y = AV, color = tss_smrt_type)) +
    geom_boxplot(position = position_dodge(), width = 0.2, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 15000))
