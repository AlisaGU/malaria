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
read_peakbase_motif_and_flank_scale_by_propor_mutation <- function(flank_length = NULL) {
    # read scale_data
    setwd("peak_base_mutation_density")
    scale_data <- read_scale_data()
    # read propor_factor
    setwd(global_dir)
    propor_data <- read_propor()

    # read motif inside
    setwd("each_site_mutation_density")
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif_scale_by_propor(site = x, scale_data = scale_data, propor_data = propor_data)
    })

    # read motif flank region
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
        data <- data.frame(apply(var_count, 1, sum) / apply(region_len, 1, sum))
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
        result[[paste0("flank_", flank)]] <- apply(propor[, 1:4], 2, function(x) {
            x / propor[, 5]
        })
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
        if (all(b[, 2] == 0)) {
            data <- b[, 3] / b[, 4] * (propor_data[[flank]][, 3] + propor_data[[flank]][, 4]) / unlist(scale_data[["C_snps"]])
        } else if (all(b[, 4] == 0)) {
            data <- b[, 1] / b[, 2] * (propor_data[[flank]][, 1] + propor_data[[flank]][, 2]) / unlist(scale_data[["A_snps"]])
        } else {
            data <- b[, 1] / b[, 2] * (propor_data[[flank]][, 1] + propor_data[[flank]][, 2]) / unlist(scale_data[["A_snps"]]) +
                b[, 3] / b[, 4] * (propor_data[[flank]][, 3] + propor_data[[flank]][, 4]) / unlist(scale_data[["C_snps"]])
        }
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
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Motif", "Flank")) +
        labs(y = "Scaled mutation density") +
        facet_wrap(~variant, nrow = 2, ncol = 1) +
        theme_bw() +
        geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
        ggtitle(title) +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
titles <- c("G,A,A/T,G,A,A/T")
names(titles) <- c("GAWGAW")

setwd(global_dir)

pdf(file = "scaled_by_propor_peakbase_motif_and_flank_global_variants.remove_close_overlap_all_snp.pdf", height = 10, width = 15)
data_for_plot <- read_peakbase_motif_and_flank_scale_by_propor_mutation(flank_length = 10)
p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot %>% filter(variant == "SNPs"), title = "GAWGAW")
print(p)
dev.off()
