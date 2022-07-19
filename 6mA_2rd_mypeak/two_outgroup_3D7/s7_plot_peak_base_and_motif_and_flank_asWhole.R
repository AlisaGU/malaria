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
read_peakbase_motif_and_flank_mutation <- function(flank_length = NULL) {
    # read motif inside
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif(site = x)
    })


    # read motif flank region
    setwd("./flank_nomatter_strand.remove_close_overlap")

    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank(flank = x)
    }


    # read peak region base mutation density
    setwd("../../peak_base_mutation_density")
    data_for_plot[[length(data_for_plot) + 1]] <- read_base_mutation_density()



    data_for_plot <- do.call(rbind, data_for_plot)
    data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
    data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data_for_plot$source <- factor(gsub("peak_", "", data_for_plot$source), levels = c("base", "motif", "motif_flank"))
    data_for_plot$site_f <- factor(data_for_plot$site, levels = unique(data_for_plot$site))
    return(data_for_plot)
}

read_each_site_motif <- function(site = NULL) {
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 1] / a[, 2])
        rownames(data) <- rownames(a)
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

read_each_site_flank <- function(flank = NULL) {
    setwd(flank)
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(flank, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 1] / a[, 2])
        colnames(data) <- c("peak_motif_flank")
        rownames(data) <- rownames(a)

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

read_base_mutation_density <- function() {
    file_names <- dir()[grep("snps|indels", dir())]
    file_names <- file_names[grep("global", file_names)]

    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 1] / a[, 2])
        colnames(data) <- c("base")
        rownames(data) <- rownames(a)

        data_for_plot <- data.frame(
            value = unlist(data),
            chrom = rep(rownames(data), times = ncol(data)),
            source = rep(colnames(data), each = nrow(data))
        )
        a <- gsub("chrom", "", data_for_plot$chrom)
        a[a != "1"] <- ""
        data_for_plot$label <- a
        data_for_plot$site <- "base"
        data_for_plot$variant <- gsub("global_", "", gsub("_mean_variant_count", "", x))

        return(data_for_plot)
    })
    result <- do.call(rbind, data)
    setwd("../")
    return(result)
}

plot_motif <- function(data_for_plot = NULL, motif = NULL, title = NULL) {
    color_value_N2N <- c("base" = "black", "motif" = "#b2182b", "motif_flank" = "#2166ac")

    p <- ggplot(data_for_plot, aes(x = site_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Base", "Motif", "Flank")) +
        stat_compare_means(
            ref.group = "base",
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 10
        ) +
        labs(y = "Mutation density") +
        facet_wrap(~variant, nrow = 2, ncol = 1) +
        theme_bw() +
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




read_peakbase_motif_and_flank_scale_mutation <- function(flank_length = NULL, motif = NULL, motif_seq = NULL, variants = NULL) {
    # read scale_factor
    setwd("../peak_base_mutation_density")
    scale_data <- read_scale_data()
    # read motif inside
    setwd(paste0("../", motif))
    motif_seq <- unlist(lapply(strsplit(motif_seq, ","), function(x) {
        gsub("/", "_", x)
    }))
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif_scale(site = x, seq = motif_seq[as.numeric(gsub("site", "", x))], scale_data = scale_data)
    })


    # read motif flank region
    setwd("./flank_nomatter_strand.remove_close_overlap")
    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank_scale(flank = x, scale_data = scale_data)
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
    file_names <- dir()[grep("snps|indels", dir())]
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

read_each_site_motif_scale <- function(site = NULL, seq = NULL, scale_data = NULL) {
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        var <- paste0(seq, "_", unlist(strsplit(gsub("_mean_variant_count", "", x), ".", fixed = T))[2])

        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 1] / a[, 2] / scale_data[[var]])
        rownames(data) <- rownames(a)
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


read_each_site_flank_scale <- function(flank = NULL, seq = "global", scale_data = NULL) {
    setwd(flank)
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(flank, file_names)]
    data <- lapply(file_names, function(x) {
        var <- paste0(seq, "_", unlist(strsplit(gsub("_mean_variant_count", "", x), ".", fixed = T))[2])

        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 1] / a[, 2] / scale_data[[var]])
        colnames(data) <- c("peak_motif_flank")
        rownames(data) <- rownames(a)

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

plot_motif_scale <- function(data_for_plot = NULL, motif = NULL, title = NULL) {
    color_value_N2N <- c("motif" = "#b2182b", "motif_flank" = "#2166ac")

    p <- ggplot(data_for_plot, aes(x = site_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Base", "Motif", "Flank")) +
        # stat_compare_means(
        #     ref.group = "base",
        #     hide.ns = T, vjust = 0.5, show.legend = FALSE,
        #     label = "p.signif", paired = T, size = 10
        # ) +
        labs(y = "Scaled mutation density") +
        scale_y_continuous(limits = c(0, 4)) +
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


read_peakbase_motif_and_flank_scale_by_propor_mutation <- function(flank_length = NULL, motif = NULL, variants = NULL) {
    # read scale_data
    setwd("../peak_base_mutation_density")
    scale_data <- read_scale_data()
    # read propor_factor
    ori_path <- getwd()
    # setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/all_motifs/", motif))
    setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/", motif))

    propor_data <- read_propor()
    setwd(ori_path)
    # scale_factor=mutation_density*propor
    scale_factor <- get_scale_factor_by_propor(propor_data = propor_data, scale_data = scale_data)

    # read motif inside
    setwd(paste0("../", motif))
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif_scale_by_propor(site = x, scale_factor = scale_factor)
    })

    # read motif flank region
    setwd("./flank_nomatter_strand.remove_close_overlap")
    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank_scale_by_propor(flank = x, scale_factor = scale_factor)
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

read_each_site_motif_scale_by_propor <- function(site = NULL, scale_factor = NULL) {
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        var <- gsub("s_mean_variant_count", "", unlist(strsplit(x, ".", fixed = T))[2])
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        data <- data.frame(a[, 1] / a[, 2] / scale_factor[[var]][[site]])
        rownames(data) <- rownames(a)
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


read_each_site_flank_scale_by_propor <- function(flank = NULL, scale_factor = NULL) {
    setwd(flank)
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(flank, file_names)]
    data <- lapply(file_names, function(x) {
        var <- gsub("s_mean_variant_count", "", unlist(strsplit(x, ".", fixed = T))[2])

        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        # a <- a[apply(a, 1, function(x) {
        #     all(x > 0)
        # }), ]
        data <- data.frame(a[, 1] / a[, 2] / scale_factor[[var]][[flank]])
        colnames(data) <- c("peak_motif_flank")
        rownames(data) <- rownames(a)

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
        # stat_compare_means(
        #     ref.group = "base",
        #     hide.ns = T, vjust = 0.5, show.legend = FALSE,
        #     label = "p.signif", paired = T, size = 10
        # ) +
        # labs(y = "Scaled (by propor) mutation density") +
        labs(y = "Scaled mutation density") +
        scale_y_continuous(limits = c(0, 5)) +
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
popu_symbol <- "ESEA_WSEA_OCE_SAM_SAS"
global_dir <- paste0(
    "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/",
    popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site"
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
# pdf(file = "scaled_by_propor_peakbase_motif_and_flank_global_variants.remove_close_overlap_all_snp.pdf", height = 10, width = 15)
# for (motif in c("GAMSAA", "GAHGAW", "GAACAA", "GAAGAA", "GATGAA", "GATGAT")) {
#     path <- paste0(global_dir, "/", motif)
#     setwd(path)
#     data_for_plot <- read_peakbase_motif_and_flank_scale_by_propor_mutation(flank_length = 10, motif = motif, variants = variants)
#     # p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot, motif = motif, title = as.character(titles[motif]))
#     p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot %>% filter(variant == "SNPs"), motif = motif, title = as.character(titles[motif]))
#     ps[[motif]] <- p
#     print(p)
# }
# dev.off()

pdf(file = "scaled_by_propor_peakbase_motif_and_flank_global_variants.remove_close_overlap_all_snp.pdf", height = 10, width = 15)
for (motif in c("GAWGAW")) {
    path <- paste0(global_dir, "/", motif)
    setwd(path)
    data_for_plot <- read_peakbase_motif_and_flank_scale_by_propor_mutation(flank_length = 10, motif = motif, variants = variants)
    # p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot, motif = motif, title = as.character(titles[motif]))
    p <- plot_motif_scale_by_propor(data_for_plot = data_for_plot %>% filter(variant == "SNPs"), motif = motif, title = as.character(titles[motif]))
    ps[[motif]] <- p
    print(p)
}
dev.off()
