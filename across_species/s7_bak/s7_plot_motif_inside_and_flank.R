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
read_motif_and_flank_mutation <- function(flank_length = NULL) {
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_motif(site = x, ratio = FALSE)
    })

    for (x in paste("site", 1:6, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_motif(site = x, ratio = TRUE)
    }
    # setwd("./flank_nomatter_strand.closest")
    setwd("./flank_nomatter_strand.remove_close_overlap")
    # setwd("./flank_nomatter_strand.remove_overlap")


    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank(flank = x, ratio = FALSE)
    }

    for (x in paste("flank_", 1:flank_length, sep = "")) {
        data_for_plot[[length(data_for_plot) + 1]] <- read_each_site_flank(flank = x, ratio = TRUE)
    }
    setwd("../")

    data_for_plot <- do.call(rbind, data_for_plot)
    data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
    data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data_for_plot$source <- factor(gsub("peak_vs_control", "pVSc", gsub("_motif", "", data_for_plot$source)), levels = c("control", "peak", "pVSc"))
    data_for_plot$site_f <- factor(data_for_plot$site, levels = unique(data_for_plot$site))
    return(data_for_plot)
}

read_each_site_motif <- function(site = NULL, ratio = NULL) {
    # file_names <- dir()[grep("allvar|snps|indels", dir())]
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 3] / a[, 4], a[, 1] / a[, 2])
        if (ratio) {
            data <- data[, 2] / data[, 1]
            data <- as.data.frame(data)
            colnames(data) <- "peak_vs_control"
        } else {
            colnames(data) <- c("control_motif", "peak_motif")
        }
        rownames(data) <- rownames(a)

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

read_each_site_flank <- function(flank = NULL, ratio = NULL) {
    setwd(flank)
    # file_names <- dir()[grep("allvar|snps|indels", dir())]
    file_names <- dir()[grep("snps|indels", dir())]

    file_names <- file_names[grep(flank, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        rownames(a) <- paste("chrom", 1:14, sep = "")
        a <- a[apply(a, 1, function(x) {
            all(x > 0)
        }), ]
        data <- data.frame(a[, 3] / a[, 4], a[, 1] / a[, 2])
        if (ratio) {
            data <- data[, 2] / data[, 1]
            data <- as.data.frame(data)
            colnames(data) <- "peak_vs_control"
        } else {
            colnames(data) <- c("control_motif", "peak_motif")
        }
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


plot_motif <- function(data_for_plot = NULL, motif = NULL, title = NULL) {
    ratio <- data_for_plot[data_for_plot$source == "pVSc", ]
    density <- data_for_plot[data_for_plot$source != "pVSc", ]
    ratio_p <- plot_each_type(data = ratio, ylab = "Mutation density ratio\n(peak/control)") +
        geom_hline(yintercept = 1, color = "red", linetype = "dashed")
    density_p <- plot_each_type(data = density, ylab = "Mutation density")

    layout <- rbind(c(1, 1, 1, 1, 2, 2, 2), c(1, 1, 1, 1, 2, 2, 2))
    p <- grid.arrange(grobs = list(density_p, ratio_p), layout_matrix = layout, top = textGrob(title, gp = gpar(col = "red", fontsize = 40)))
    return(p)
}

plot_each_type <- function(data = NULL, ylab = NULL) {
    color_value_N2N <- c("control" = "black", "peak" = "#b2182b", "pVSc" = "#dc863b")

    p <- ggplot(data, aes(x = site_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Control", "Peak", "pVSc")) +
        stat_compare_means(aes(group = source),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 10
        ) +
        labs(y = ylab) +
        facet_wrap(~variant, nrow = 3, ncol = 1, scale = "free_y") +
        theme_bw() +
        # ggtitle(motif) +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 30, color = "black"),
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
titles <- c("G,A,A/C,G/C,A,A", "G,A,A/T/C,G,A,A/T")
names(titles) <- c("GAMSAA", "GAHGAW")


setwd(global_dir)
pdf(file = "motif_and_flank_global_variants.remove_close_overlap.pdf", height = 20, width = 60)
for (motif in c("GAHGAW", "GAMSAA")) {
    path <- paste0(global_dir, "/", motif)
    setwd(path)
    data_for_plot <- read_motif_and_flank_mutation(flank_length = 10)
    p <- plot_motif(data_for_plot = data_for_plot, motif = motif, title = titles[motif])
    print(p)
}
dev.off()




snp <- data_for_plot %>%
    filter(variant == "SNPs") %>%
    filter(source == "pVSc")
a <- aov(value ~ site, snp)
summary(a)
pair <- TukeyHSD(a)
pair_data <- as.data.frame(pair$site)
pair_data[pair_data[, 4] < 0.05, ]


indel <- data_for_plot %>%
    filter(variant == "Indels") %>%
    filter(source == "pVSc")
a <- aov(value ~ site, indel)
summary(a)
pair <- TukeyHSD(a)
pair_data <- as.data.frame(pair$site)
pair_data[pair_data[, 4] < 0.05, ]
