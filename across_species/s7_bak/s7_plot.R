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
# 2. functions ------------------------------------------------------------ TODO:
plot_motifs <- function(popu_symbol = NULL, motif = NULL) {
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/", motif)
    setwd(path)
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site(site = x)
    })
    data_for_plot <- do.call(rbind, data_for_plot)
    data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
    data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data_for_plot$source <- factor(data_for_plot$source, levels = c("control_motif", "peak_motif"))
    data_for_plot$site_f <- factor(data_for_plot$site, levels = unique(data_for_plot$site))

    p <- plot_motif(data_for_plot = data_for_plot, motif = motif)
    return(p)
}

plot_motifs_with_count <- function(popu_symbol = NULL, motif = NULL) {
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/", motif)
    setwd(path)
    data_for_plot <- lapply(paste("site", 1:6, sep = ""), function(x) {
        read_each_site_count(site = x)
    })
    data_for_plot <- do.call(rbind, data_for_plot)
    data_for_plot$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data_for_plot$variant)))))
    data_for_plot$variant_f <- factor(data_for_plot$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data_for_plot$source <- factor(data_for_plot$source, levels = c("control_motif_count", "peak_motif_count", "control_motif_value", "peak_motif_value"))
    data_for_plot$site_f <- factor(data_for_plot$site, levels = unique(data_for_plot$site))

    p <- plot_motif_with_mutation_count(data_for_plot = data_for_plot, motif = motif)
    return(p)
}

read_each_site <- function(site = NULL) {
    # file_names <- dir()[grep(site, dir())]
    file_names <- dir()[grep("allvar|snps|indels", dir())]
    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        data <- data.frame(a[, 3] / a[, 4], a[, 1] / a[, 2])
        colnames(data) <- c("control_motif", "peak_motif")

        rownames(data) <- paste("chrom", 1:14, sep = "")

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

read_each_site_count <- function(site = NULL) {
    # file_names <- dir()[grep(site, dir())]
    file_names <- dir()[grep("allvar|snps|indels", dir())]
    file_names <- file_names[grep(site, file_names)]
    data <- lapply(file_names, function(x) {
        a <- read.table(x, header = F, as.is = T)
        data <- data.frame(a[, 3], a[, 1], a[, 3] / a[, 4], a[, 1] / a[, 2])
        colnames(data) <- c("control_motif_count", "peak_motif_count", "control_motif_value", "peak_motif_value")

        rownames(data) <- paste("chrom", 1:14, sep = "")

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
plot_motif <- function(data_for_plot = NULL, motif = NULL) {
    color_value_N2N <- c("control_motif" = "black", "peak_motif" = "#b2182b")
    p <- ggplot(data_for_plot, aes(x = variant_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Control motif", "Peak motif")) +
        stat_compare_means(aes(group = source),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 10
        ) +
        labs(y = "Mutation density") +
        facet_wrap(~site, nrow = 2, ncol = 3, scale = "free_x") +
        theme_bw() +
        ggtitle(motif) +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
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

plot_motif_with_mutation_count <- function(data_for_plot = NULL, motif = NULL) {
    color_value_N2N <- c("control_motif_value" = "black", "peak_motif_value" = "#b2182b")
    main_plot <- ggplot(data_for_plot %>% filter(source == "control_motif_value" | source == "peak_motif_value"), aes(x = variant_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Control motif", "Peak motif")) +
        stat_compare_means(aes(group = source),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 10
        ) +
        facet_wrap(~site, nrow = 2, ncol = 3, scale = "free_x") +
        theme_bw() +
        ggtitle(motif) +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.text.y = element_text(size = 30, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )

    for (index in paste0("site", 1:6)) {
        inset_plot <- get_inset(data_for_plot = data_for_plot %>% filter((source == "control_motif_count" | source == "peak_motif_count") & site == index))
        main_plot <- main_plot +
            annotation_custom2(
                grob = ggplotGrob(inset_plot),
                data = data.frame(site = index),
                ymin = 0.075, ymax = Inf, xmin = 1.5, xmax = 3.5
            )
    }



    return(main_plot)
}

get_inset <- function(data_for_plot = NULL) {
    p <- ggplot(data = data_for_plot %>% filter(source == "control_motif_count" | source == "peak_motif_count"), aes(x = value, group = variant_f, color = variant_f)) +
        geom_density() +
        theme_bw() +
        theme(
            # panel.border = element_blank(),
            # axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 13, color = "black"),
            axis.text.y = element_blank(),
            panel.background = element_rect(
                fill = "transparent",
                colour = NA_character_
            ), # necessary to avoid drawing panel outline
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            plot.background = element_rect(
                fill = "transparent",
                colour = NA_character_
            ),
            # plot.title = element_text(
            #     colour = "black",
            #     size = 40, vjust = 0.5, hjust = 0.5
            # ),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            # legend.title = element_blank(),
            # legend.text = element_text(size = 24, color = "black"),
            legend.position = "none",
            # plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            # panel.spacing = unit(3, "lines")
        )
    return(p)
}

annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
    layer(
        data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(
            grob = grob,
            xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax
        )
    )
}


plot_flank_site_patterns <- function(popu_symbol = NULL, motif = NULL) {
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/", motif)
    setwd(path)

    data_for_plot <- list()
    # for (x in paste("flank_", c("upstream", "downstream", "nomatter_strand"), sep = "")) {
    for (x in paste("flank_", c("nomatter_strand"), sep = "")) {
        setwd(x)
        a <- list()
        for (y in paste("flank_", 1:18, sep = "")) {
            setwd(y)
            a[[y]] <- read_each_site(site = y)
            setwd("../")
        }
        a <- do.call(rbind, a)
        data_for_plot[[x]] <- a
        setwd("../")
    }

    # p_upstream <- plot_each_flank_type(data = data_for_plot[["flank_upstream"]])
    # p_downstream <- plot_each_flank_type(data = data_for_plot[["flank_downstream"]])
    # p_nomatter_strand <- plot_each_flank_type(data = data_for_plot[["flank_nomatter_strand"]])
    # p <- p_upstream / p_downstream / p_nomatter_strand +
    #     plot_annotation(
    #         title = motif,
    #         theme = theme(plot.title = element_text(size = 36, hjust = 0.5, color = "red"))
    #     )
    # return(p)

    p_nomatter_strand <- plot_each_flank_type(data = data_for_plot[["flank_nomatter_strand"]])
    p_nomatter_strand +
        plot_annotation(
            title = motif,
            theme = theme(plot.title = element_text(size = 36, hjust = 0.5, color = "red"))
        )
    return(p_nomatter_strand)
}

plot_flank_site_patterns_boudary <- function(popu_symbol = NULL, motif = NULL, nrow = NULL, ncol = NULL) {
    path <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/", motif)
    setwd(path)

    data_for_plot <- list()
    for (x in paste("flank_", c("nomatter_strand"), sep = "")) {
        setwd(x)
        a <- list()
        for (y in paste("flank_", 1:24, sep = "")) {
            setwd(y)
            a[[y]] <- read_each_site(site = y)
            setwd("../")
        }
        a <- do.call(rbind, a)
        data_for_plot[[x]] <- a
        setwd("../")
    }

    p_nomatter_strand <- plot_each_flank_type(data = data_for_plot[["flank_nomatter_strand"]], nrow = nrow, ncol = ncol)
    p_nomatter_strand +
        plot_annotation(
            title = motif,
            theme = theme(plot.title = element_text(size = 36, hjust = 0.5, color = "red"))
        )
    return(p_nomatter_strand)
}

plot_each_flank_type <- function(data = NULL, nrow = NULL, ncol = NULL) {
    data$variant <- gsub("snps", "SNPs", gsub("indels", "Indels", gsub("allvar", "All variants", gsub("_single_del", " to sDel", gsub("2", " to ", data$variant)))))
    data$variant_f <- factor(data$variant, levels = c(
        "All variants", "SNPs", "Indels",
        "A to G", "A to T", "A to C", "A to sDel",
        "G to A", "G to T", "G to C", "G to sDel",
        "C to T", "C to A", "C to G", "C to sDel",
        "T to C", "T to A", "T to G", "T to sDel"
    ))
    data$source <- factor(data$source, levels = c("control_motif", "peak_motif"))
    data$site_f <- factor(data$site, levels = unique(data$site))

    p <- plot_motif_flank(data_for_plot = data, nrow = nrow, ncol = ncol)
    return(p)
}

plot_motif_flank <- function(data_for_plot = NULL, nrow = NULL, ncol = NULL) {
    color_value_N2N <- c("control_motif" = "black", "peak_motif" = "#b2182b")
    p <- ggplot(data_for_plot, aes(x = variant_f, y = value, color = source)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(), size = 4) +
        scale_color_manual(values = color_value_N2N, labels = c("Control motif", "Peak motif")) +
        stat_compare_means(aes(group = source),
            hide.ns = T, vjust = 0.5, show.legend = FALSE,
            label = "p.signif", paired = T, size = 10
        ) +
        labs(y = "Mutation density") +
        facet_wrap(~site_f, nrow = nrow, ncol = ncol, scale = "free_x") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
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
Args <- commandArgs()
popu_symbol <- Args[6]

motifs <- c("GAMSAA", "GAHGAW", "GAACAA", "GAAGAA", "GAAGAT", "GACCAA", "GACGAA", "GACGAT", "GATGAA", "GATGAT")
collapsed_motifs <- c("GAMSAA", "GAHGAW", "motifs_summary")
# 4. variable setting of test module--------------------------------------- TODO:

# popu_symbol <- "ESEA_WSEA_OCE_SAM_SAS"
popu_symbol <- "WSEA"

# 5. process -------------------------------------------------------------- TODO:
# output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/motifs_each_site_mutation_pattern.pdf")
# pdf(output, width = 24, height = 16)
# for (motif in motifs) {
#     p <- plot_motifs(popu_symbol = popu_symbol, motif = motif)
#     plot(p)
# }
# dev.off()


# output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/motifs_each_site_global_variants.pdf")
# pdf(output, width = 24, height = 16)
# for (motif in motifs) {
#     p <- plot_motifs(popu_symbol = popu_symbol, motif = motif)
#     plot(p)
# }
# dev.off()


# output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/motifs_each_site_global_variants_with_count.pdf")
# pdf(output, width = 24, height = 16)
# for (motif in motifs) {
#     p <- plot_motifs_with_count(popu_symbol = popu_symbol, motif = motif)
#     plot(p)
# }
# dev.off()


# output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/collapsed_motifs_each_site_global_variants.pdf")
# pdf(output, width = 24, height = 16)
# for (motif in collapsed_motifs) {
#     p <- plot_motifs(popu_symbol = popu_symbol, motif = motif)
#     plot(p)
# }
# dev.off()


output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/", popu_symbol, "/motif/mutation_dentisy_two_outgroup_consistent/each_site/collapsed_motifs_flank_sites_global_variants.pdf")
rowth <- 4
colth <- 6
height <- 10 * rowth
width <- 10 * colth
pdf(output, width = width, height = height)
for (motif in "motifs_summary") {
    p <- plot_flank_site_patterns_boudary(popu_symbol = popu_symbol, motif = motif, nrow = rowth, ncol = colth)
    plot(p)
}
dev.off()
