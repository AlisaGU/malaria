#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(data.table)
library(purrr)
library(gridExtra)

# 2. functions ------------------------------------------------------------ TODO:
read_peakbase_motif_and_flank_scale_by_propor_mutation <- function(flank_length = NULL, motif = NULL, variants = NULL) {
    ## read scale_data 碱基在全基因组常染色体区域的平均突变速率
    setwd("../wholeGenome_Euchromatin_twoOutgroup_base_mutation_density")
    scale_data <- read_scale_data()
    ## read propor_factor 对应位置上，碱基的比例
    ori_path <- getwd()
    if (motif == "GAWGAW") {
        setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/", motif))
    } else {
        setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_all_motifs_asWhole/all_motifs/", motif))
    }
    propor_data <- read_propor()
    setwd(ori_path)
    # scale_factor <- get_scale_factor_by_propor(propor_data = propor_data, scale_data = scale_data)

    ## read motif inside
    setwd(paste0("../", motif))
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif_scale_by_propor(site = x, scale_data = scale_data, propor_data = propor_data)
    })

    ## read motif flank region
    setwd("./flank_nomatter_strand.remove_close_overlap")
    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank_scale_by_propor(flank = x, scale_data = scale_data, propor_data = propor_data)
    }

    data_for_plot <- do.call(rbind, data_for_plot)
    data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
    data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data_for_plot$source <- factor(gsub("peak_", "", data_for_plot$source), levels = c("motif", "motif_flank"))
    data_for_plot$site_f <- factor(data_for_plot$site, levels = unique(data_for_plot$site))
    return(data_for_plot)
}

read_scale_data <- function() {
    file_names <- dir()[grep("snps", dir())]

    result <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        var_count_column <- seq(1, ncol(a), by = 2)
        region_len_column <- seq(2, ncol(a), by = 2)
        var_count <- as.matrix(a[, var_count_column])
        region_len <- as.matrix(a[, region_len_column])
        # data <- data.frame(apply(var_count, 1, sum) / apply(region_len, 1, sum))
        data <- data.frame(rep(sum(apply(var_count, 1, sum)) / sum(apply(region_len, 1, sum)), 14))
        rownames(data) <- paste("chrom", 1:14, sep = "")
        colnames(data) <- "scale_factor"
        return(data)
    })
    names(result) <- gsub("_mean_variant_count", "", file_names)
    return(result)
}


read_propor <- function() {
    result <- lapply(1:6, function(x) {
        propor <- read.table(paste0("site", x, "/site", x, ".peak.base_count"), header = F, as.is = T)
        colnames(propor) <- c("A", "T", "C", "G", "all")
        rownames(propor) <- paste("chrom", 1:14, sep = "")
        apply(propor[, 1:4], 2, function(x) {
            x / propor[, 5]
        })
    })
    names(result) <- paste("site", 1:6, sep = "")
    for (flank in 1:10) {
        propor <- read.table(paste0("flank_nomatter_strand.remove_close_overlap/flank_", flank, "/flank", flank, ".peak.base_count"), header = F, as.is = T)
        colnames(propor) <- c("A", "T", "C", "G", "all")
        rownames(propor) <- paste("chrom", 1:14, sep = "")
        result[[paste0("flank_", flank)]] <- t(apply(propor[, 3:4], 1, function(x) {
            x / sum(x)
        }))
    }

    return(result)
}

read_each_site_motif_scale_by_propor <- function(site = NULL, scale_data = NULL, propor_data = NULL) {
    file_names <- dir()[grep("snps", dir())]

    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        var <- gsub("s_mean_variant_count", "", unlist(strsplit(x, ".", fixed = T))[2])
        a <- fread(x, header = F, stringsAsFactors = T, select = 1:8)
        b <- cbind(a[, 1] + a[, 3], a[, 2] + a[, 4], a[, 5] + a[, 7], a[, 6] + a[, 8])
        rownames(b) <- paste("chrom", 1:14, sep = "")

        data <- ""
        if (all(b[, 2] == 0)) {
            data <- b[, 3] / b[, 4] * (propor_data[[site]][, 3] + propor_data[[site]][, 4]) / unlist(scale_data[["C_snps"]])
        } else if (all(b[, 4] == 0)) {
            data <- b[, 1] / b[, 2] * (propor_data[[site]][, 1] + propor_data[[site]][, 2]) / unlist(scale_data[["A_snps"]])
        } else {
            data <- (propor_data[[site]][, 1] + propor_data[[site]][, 2]) * b[, 1] / b[, 2] / unlist(scale_data[["A_snps"]]) +
                b[, 3] / b[, 4] * (propor_data[[site]][, 3] + propor_data[[site]][, 4]) / unlist(scale_data[["C_snps"]])
        }

        rownames(data) <- rownames(b)
        colnames(data) <- "peak_motif"
        data_for_plot <- data.frame(
            value = unlist(data),
            chrom = rep(rownames(data), times = ncol(data)),
            source = rep(colnames(data), each = nrow(data))
        )
        a <- gsub("chrom", "", data_for_plot$chrom)
        a[a != "1"] <- ""
        data_for_plot$label <- a
        data_for_plot$site <- site
        data_for_plot$variant <- gsub("_mean_variant_count", "", gsub(paste0(site, "."), "", x))

        return(data_for_plot)
    })
    result <- do.call(rbind, data)
    return(result)
}

get_scale_factor_by_propor <- function(propor_data = NULL, scale_data = NULL) {
    A_T_C_G_peak_snp_base_mutation_density <- data.frame(
        "A" = scale_data[["A_snps"]], "T" = scale_data[["T_snps"]],
        "C" = scale_data[["C_snps"]], "G" = scale_data[["G_snps"]]
    )

    A_T_C_G_peak_indel_base_mutation_density <- data.frame(
        A = scale_data[["A_indels"]], T = scale_data[["T_indels"]],
        C = scale_data[["C_indels"]], G = scale_data[["G_indels"]]
    )

    result_snp <- lapply(propor_data, function(x) {
        result <- as.data.frame(apply(x * A_T_C_G_peak_snp_base_mutation_density, 1, sum))
        colnames(result) <- "scale_factor"
        return(result)
    })

    result_indel <- lapply(propor_data, function(x) {
        result <- as.data.frame(apply(x * A_T_C_G_peak_indel_base_mutation_density, 1, sum))
        colnames(result) <- "scale_factor"
        return(result)
    })

    result <- list(snp = result_snp, indel = result_indel)

    return(result)
}


read_each_site_flank_scale_by_propor <- function(flank = NULL, scale_data = NULL, propor_data = NULL) {
    setwd(flank)
    file_names <- dir()[grep("snps", dir())]

    file_names <- file_names[grep(flank, file_names)]
    data <- lapply(file_names, function(x) {
        var <- gsub("s_mean_variant_count", "", unlist(strsplit(x, ".", fixed = T))[2])

        a <- fread(x, header = F, stringsAsFactors = T, select = 1:8)
        b <- cbind(a[, 1] + a[, 3], a[, 2] + a[, 4], a[, 5] + a[, 7], a[, 6] + a[, 8])
        rownames(b) <- paste("chrom", 1:14, sep = "")

        data <- ""
        data <- b[, 3] / b[, 4] * (propor_data[[flank]][, "C"] + propor_data[[flank]][, "G"]) / unlist(scale_data[["C_snps"]])
        rownames(data) <- rownames(b)
        colnames(data) <- c("peak_motif_flank")

        data_for_plot <- data.frame(
            value = unlist(data),
            chrom = rep(rownames(data), times = ncol(data)),
            source = rep(colnames(data), each = nrow(data))
        )
        a <- gsub("chrom", "", data_for_plot$chrom)
        a[a != "1"] <- ""
        data_for_plot$label <- a
        data_for_plot$site <- flank
        data_for_plot$variant <- gsub("_mean_variant_count", "", gsub(paste0(flank, "."), "", x))

        return(data_for_plot)
    })
    result <- do.call(rbind, data)
    setwd("../")
    return(result)
}


plot_motif_scale_by_propor <- function(data_for_plot = NULL, motif = NULL, title = NULL) {
    color_value_N2N <- c("motif" = "#b2182b", "motif_flank" = "#2166ac")

    p <- ggplot(data_for_plot, aes(x = site_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA, width = 0.5, size = 1) +
        geom_point(size = 6, position = position_jitterdodge()) +
        scale_color_manual(values = color_value_N2N, labels = c("Motif", "Flank")) +
        # stat_compare_means(
        #     ref.group = "base",
        #     hide.ns = T, vjust = 0.5, show.legend = FALSE,
        #     label = "p.signif", paired = T, size = 10
        # ) +
        # labs(y = "Scaled (by propor) mutation density") +
        labs(y = "Relative mutation density") +
        # scale_y_continuous(limits = c(0, 5)) +
        # facet_wrap(~variant, nrow = 2, ncol = 1) +
        scale_x_discrete(labels = c(paste("S", 1:6, sep = ""), paste("F", 1:10, sep = ""))) +
        theme_bw() +
        geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
        ggtitle(title) +
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
                size = 36, color = "black",
                # angle = 45, hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = element_text(size = 36, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 30, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}


compute_significance <- function(data = NULL) {
    data_split <- split(data, f = data$site)
    lapply(data_split, function(x) {
        value <- x$value

        distri <- sapply(1:100000, function(seed) {
            set.seed(seed)
            sample_value <- sample(value, 14, replace = TRUE)
            statisitic <- mean(sample_value)
            return(statisitic)
        })
        hist(distri)
        distri_sort <- sort(distri, decreasing = FALSE)
        lowCI <- distri_sort[10000 * 0.05]
        highCI <- distri_sort[10000 * 0.95]
        meanValue <- mean(distri_sort)
        return(c(lowCI, highCI, meanValue))
    })
}
# 3. input ---------------------------------------------------------------- TODO:
popu_symbol <- "ESEA_WSEA_OCE_SAM_SAS"
global_dir <- paste0(
    "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/",
    popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site_GAWGAW_asWhole_LI_formula"
)

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
variants <- c("G", "A")
titles <- c(
    "G,A,A/C,G/C,A,A", "G,A,A/T/C,G,A,A/T", "G,A,A,C,A,A",
    "G,A,A,G,A,A", "G,A,T,G,A,A", "G,A,T,G,A,T", "G,A,A/T,G,A,A/T"
)
names(titles) <- c(
    "GAMSAA", "GAHGAW", "GAACAA",
    "GAAGAA", "GATGAA", "GATGAT", "GAWGAW"
)

setwd(global_dir)
ps <- list()


pdf(file = "scaled_by_propor_peakbase_motif_and_flank_global_variants.remove_close_overlap_all_snp_submotif.pdf", height = 8, width = 16)
# for (motif in c("GAMSAA", "GAHGAW", "GAACAA", "GAAGAA", "GATGAA", "GATGAT")) {
for (motif in c("GAAGAA", "GATGAA", "GATGAT")) {
    path <- paste0(global_dir, "/../each_site_all_motifs_asWhole_LI_formula/", motif)
    setwd(path)
    data_for_plot <- read_peakbase_motif_and_flank_scale_by_propor_mutation(flank_length = 10, motif = motif, variants = variants)
    # p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot, motif = motif, title = as.character(titles[motif]))
    p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot %>% filter(variant == "SNPs"), motif = motif, title = "")
    ps[[motif]] <- p
    print(p)
}
dev.off()

pdf(file = "scaled_by_propor_peakbase_motif_and_flank_global_variants.remove_close_overlap_all_snp.pdf", height = 8, width = 16)
for (motif in c("GAWGAW")) {
    path <- paste0(global_dir, "/", motif)
    setwd(path)
    data_for_plot <- read_peakbase_motif_and_flank_scale_by_propor_mutation(flank_length = 10, motif = motif, variants = variants)
    # p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot, motif = motif, title = as.character(titles[motif]))
    p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot %>% filter(variant == "SNPs"), motif = motif, title = "")
    ps[[motif]] <- p
    print(p)
}
dev.off()
