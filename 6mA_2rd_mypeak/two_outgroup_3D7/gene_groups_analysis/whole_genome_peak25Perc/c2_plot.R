#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(grid)
# 2. functions ------------------------------------------------------------ TODO:
read_peak_proportion_and_mutation_density <- function(gene_class_criterion = NULL, flank_length = NULL) {
    snp_mutation_density <- read.table(paste0(mutation_density_global_dir, "/snps_mean_variant_count_flank", flank_length), as.is = T, stringsAsFactors = F)
    colnames(snp_mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    indel_mutation_density <- read.table(paste0(mutation_density_global_dir, "/indels_mean_variant_count_flank", flank_length), as.is = T, stringsAsFactors = F)
    colnames(indel_mutation_density) <- c("gene_name", "mutation_count", "cds_length")
    mutation_density <- rbind(data.frame(snp_mutation_density, type = "SNP"), data.frame(indel_mutation_density, type = "Indel"))

    mutation_density_CDS <- data.frame(
        gene_name = mutation_density$gene_name,
        mutation_density = mutation_density$mutation_count / mutation_density$cds_length,
        type = mutation_density$type
    )


    peak_length <- read.table(paste0(peak_proportion_global_dir, "/peak_control_length_in_flank", flank_length, "_peak25Perc"), as.is = T, stringsAsFactors = F)
    colnames(peak_length) <- c("gene_name", "peak_length", "control_length", "peak_count")
    peak_proportion <- data.frame(
        gene_name = peak_length$gene_name,
        peak_proportion = apply(peak_length[, 2:3], 1, function(x) {
            x[1] / sum(x)
        }),
        peak_count = peak_length$peak_count
    )

    gene_class <- ""
    if (gene_class_criterion == "h3k9me3") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class_h3k9me3.bed"), as.is = T, stringsAsFactors = F)
    } else if (gene_class_criterion == "core") {
        gene_class <- read.table(paste0(peak_proportion_global_dir, "/genes_class.bed"), as.is = T, stringsAsFactors = F)
    }
    gene_class$length <- gene_class$V3 - gene_class$V2

    colnames(gene_class)[c(4, 6)] <- c("gene_name", "gene_class")
    gene_class <- gene_class[, c(4, 6, 7)]
    result <- merge(merge(mutation_density_CDS, peak_proportion, by = "gene_name"), gene_class, by = "gene_name")
    return(result)
}

read_nonPseudo_gene_list <- function() {
    result <- lapply(c("VAR", "STEVOR", "RIF"), function(gene_group) {
        data <- read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12/", gene_group, "/", gene_group, "_genes_list"))
        result <- data.frame(gene_name = unlist(data), gene_group = gene_group)
        return(result)
    })
    result <- do.call(rbind, result)
    return(result)
}


read_includingPseudo_gene_list <- function() {
    result <- lapply(c("VAR", "STEVOR", "RIF"), function(gene_group) {
        data <- read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_include_pseudo/", gene_group, "/", gene_group, "_genes_list"))
        result <- data.frame(gene_name = unlist(data), gene_group = gene_group)
        return(result)
    })
    result <- do.call(rbind, result)
    return(result)
}

plot_scatter <- function(data_for_plot = NULL, flank_length = NULL, xtype = NULL) {
    x_name <- ""
    if (xtype == "peak_proportion") {
        x_name <- ifelse(flank_length == "0kb",
            "Proportion of peak in gene",
            paste0("Proportion of peak\nin gene and flank ", flank_length)
        )
    } else {
        x_name <- ifelse(flank_length == "0kb",
            "Peak count in gene",
            paste0("Peak count\nin gene and flank ", flank_length)
        )
    }

    p <- ""
    if (xtype == "peak_proportion") {
        p <- ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density, group = type, color = type))
    } else {
        p <- ggplot(data_for_plot, aes(x = peak_count, y = mutation_density, group = type, color = type))
    }


    p <- p +
        geom_point(size = 4, alpha = 0.3) +
        stat_smooth(aes(color = type), method = "lm", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, show.legend = FALSE) +
        labs(x = x_name, y = "Mutation density") +
        scale_color_manual(values = c(
            "SNP" = rgb(190, 66, 65, maxColorValue = 255),
            "Indel" = rgb(24, 109, 167, maxColorValue = 255)
        )) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_different_flank_length <- function(peak_proportion_threshold = NULL, gene_length_threshold = NULL) {
    all_p <- list()
    for (gene_class_criterion in c("core")) {
        for (flank_length in c("0kb", "1kb", "2kb")) {
            data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)

            data_core <- data %>% filter(gene_class == "core")
            data_noncore <- data %>% filter(gene_class == "noncore")


            for (xtype in c("peak_proportion", "peak_count")) {
                core_p <- plot_scatter(
                    data_for_plot = data_core %>% filter(peak_proportion > peak_proportion_threshold) %>% filter(length > gene_length_threshold),
                    flank_length = flank_length, xtype = xtype
                ) +
                    ggtitle(ifelse(gene_class_criterion == "core", "Core genome", "Non-H3K9me3 region"))
                noncore_p <- plot_scatter(
                    data_for_plot = data_noncore %>% filter(peak_proportion > peak_proportion_threshold) %>% filter(length > gene_length_threshold),
                    flank_length = flank_length, xtype = xtype
                ) +
                    ggtitle(ifelse(gene_class_criterion == "core", "Non-core genome", "H3K9me3 region"))

                all_p[[gene_class_criterion]][[flank_length]][[xtype]][["core"]] <- core_p
                all_p[[gene_class_criterion]][[flank_length]][[xtype]][["noncore"]] <- noncore_p
            }
        }
    }
    return(all_p)
}

plot_var_rifin_stevor_as_one_group <- function(data_for_plot, x_name = NULL) {
    p <- ""

    if (x_name == "peak_proportion") {
        p <- ggplot(data_for_plot %>% filter(peak_proportion > 0), aes(x = peak_proportion, y = mutation_density, color = type))
    } else {
        p <- ggplot(data_for_plot %>% filter(peak_count > 0), aes(x = peak_count, y = mutation_density, color = type))
    }

    p <- p +
        geom_point(aes(shape = gene_group), size = 6) +
        stat_smooth(method = "lm", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, show.legend = FALSE) +
        labs(x = ifelse(x_name == "peak_proportion", "Peak proportion", "Peak count"), y = "Mutation density") +
        scale_color_manual(values = c(
            "SNP" = rgb(190, 66, 65, maxColorValue = 255),
            "Indel" = rgb(24, 109, 167, maxColorValue = 255)
        )) +
        scale_shape_manual(values = c("VAR" = 3, "RIF" = 17, "STEVOR" = 16)) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_var_rifin_stevor_as_many_groups <- function(data_for_plot, x_name = NULL) {
    p <- ""

    if (x_name == "peak_proportion") {
        p <- ggplot(data_for_plot %>% filter(peak_proportion > 0), aes(x = peak_proportion, y = mutation_density))
    } else {
        p <- ggplot(data_for_plot %>% filter(peak_count > 0), aes(x = peak_count, y = mutation_density))
    }

    p <- p +
        geom_point(aes(color = gene_group), size = 6, alpha = 0.6) +
        stat_smooth(aes(color = gene_group, group = gene_group), method = "lm", se = FALSE) +
        stat_cor(aes(color = gene_group, group = gene_group), method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, show.legend = FALSE) +
        labs(x = ifelse(x_name == "peak_proportion", "Peak proportion", "Peak count"), y = "Mutation density") +
        facet_wrap(~type) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 44, color = "black", face = "bold"),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

plot_var_rifin_stevor_seperately <- function(data_for_plot, x_name = NULL) {
    ps <- list()
    for (var_type in c("SNP", "Indel")) {
        p <- ""

        if (x_name == "peak_proportion") {
            p <- ggplot(data_for_plot %>% filter(type == var_type) %>% filter(peak_proportion > 0), aes(x = peak_proportion, y = mutation_density))
        } else {
            p <- ggplot(data_for_plot %>% filter(type == var_type) %>% filter(peak_count > 0), aes(x = peak_count, y = mutation_density))
        }

        p <- p +
            geom_point(color = rgb(101, 149, 192, maxColorValue = 255), size = 6, alpha = 0.6) +
            # geom_point(aes(color = gene_class), size = 6, alpha = 0.6) +

            stat_smooth(method = "lm", se = FALSE) +
            stat_cor(method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, show.legend = FALSE) +
            labs(x = ifelse(x_name == "peak_proportion", "Peak proportion", "Peak count"), y = "Mutation density") +
            facet_wrap(~gene_group) +
            ggtitle(var_type) +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                panel.border = element_blank(),
                # strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
                strip.text.x = element_text(size = 44, color = "black", face = "bold"),
                axis.line = element_line(colour = "black"),
                # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
                axis.text.x = element_text(size = 36, color = "black"),
                axis.text.y = element_text(size = 36, color = "black"),
                plot.title = element_text(
                    colour = "black",
                    size = 40, vjust = 0.5, hjust = 0.5
                ),
                axis.title.x = element_text(size = 40, color = "black"),
                axis.title.y = element_text(size = 40, color = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 24, color = "black"),
                legend.position = "bottom",
                plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
                panel.spacing = unit(3, "lines")
            )

        # if(x_name == "peak_proportion"){
        #     p<-p+scale_x_log10()
        # }
        ps[[var_type]] <- p
    }
    result <- ps$SNP / ps$Indel
    return(result)
}

plot_boxplot_youwu6mA <- function(data_for_plot = NULL, title = NULL) {
    data_for_plot$has_peak <- data_for_plot$peak_proportion == 0

    pvalue_snp <- t.test(
        data_for_plot$mutation_density[data_for_plot$has_peak == TRUE & data_for_plot$type == "SNP"],
        data_for_plot$mutation_density[data_for_plot$has_peak == FALSE & data_for_plot$type == "SNP"]
    )$p.value
    pvalue_snp <- format(pvalue_snp, digits = 2)

    pvalue_indel <- t.test(
        data_for_plot$mutation_density[data_for_plot$has_peak == TRUE & data_for_plot$type == "Indel"],
        data_for_plot$mutation_density[data_for_plot$has_peak == FALSE & data_for_plot$type == "Indel"]
    )$p.value
    pvalue_indel <- format(pvalue_indel, digits = 2)

    ypos <- max(data_for_plot$mutation_density) * 0.8

    data_for_plot$type <- factor(data_for_plot$type, levels = c("SNP", "Indel"))
    p <- ggplot(data_for_plot, aes(x = type, y = mutation_density)) +
        geom_boxplot(aes(fill = has_peak)) +
        ggtitle(title) +
        labs(y = "Mutation density") +
        annotate("text", x = 1, y = ypos, label = pvalue_snp, size = 12) +
        annotate("text", x = 2, y = ypos, label = pvalue_indel, size = 12) +
        # stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method = "t.test", size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        scale_fill_manual(values = c(rgb(137, 157, 159, maxColorValue = 255), rgb(236, 174, 34, maxColorValue = 255)), labels = c("Without peak", "With peak")) +
        scale_y_log10() +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

get_outlier_genes <- function(peak_proportion_threshold = NULL, gene_length_threshold = NULL,
                              flank_length = NULL, gene_class_criterion = NULL) {
    data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
    data_for_plot <- data %>%
        filter(gene_class == "core") %>%
        filter(peak_proportion > peak_proportion_threshold) %>%
        filter(length > gene_length_threshold)
    data_for_plot$gene_color <- apply(data_for_plot, 1, function(x) {
        if (x["mutation_density"] > 0.07 & x["peak_proportion"] <= 0.2) {
            return("I")
        } else if (x["mutation_density"] > 0.07 & x["peak_proportion"] > 0.2) {
            return("II")
        } else if (x["mutation_density"] <= 0.07 & x["peak_proportion"] > 0.2) {
            return("III")
        } else {
            return("IV")
        }
    })

    ggplot(data_for_plot, aes(x = peak_proportion, y = mutation_density, group = type, color = type)) +
        geom_point(aes(shape = gene_color), size = 4, alpha = 0.3) +
        stat_smooth(aes(color = type), method = "lm", se = FALSE) +
        stat_cor(method = "spearman", label.x.npc = 0.04, label.y.npc = 0.9, size = 12, show.legend = FALSE) +
        labs(x = "Peak proportion", y = "Mutation density") +
        scale_color_manual(values = c(
            "SNP" = rgb(190, 66, 65, maxColorValue = 255),
            "Indel" = rgb(24, 109, 167, maxColorValue = 255)
        )) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        ) +
        ggtitle("Euchromatin") +
        geom_hline(yintercept = 0.07, color = "red", size = 2) +
        geom_vline(xintercept = 0.2, color = "green", size = 2) +
        annotate("text", x = 0.1, y = 0.12, label = "I", size = 30) +
        annotate("text", x = 0.4, y = 0.12, label = "II", size = 30) +
        annotate("text", x = 0.4, y = 0.05, label = "III", size = 30) +
        annotate("text", x = 0.1, y = 0.05, label = "IV", size = 30)


    gene_I <- unique(data_for_plot$gene_name[data_for_plot$gene_color == "I"])
    gene_II <- unique(data_for_plot$gene_name[data_for_plot$gene_color == "II"])
    gene_III <- unique(data_for_plot$gene_name[data_for_plot$gene_color == "III"])
    gene_IV <- unique(data_for_plot$gene_name[data_for_plot$gene_color == "IV"])
    write.table(gene_I, paste0(mutation_density_global_dir, "/", "I_class_gene_list.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(gene_II, paste0(mutation_density_global_dir, "/", "II_class_gene_list.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(gene_III, paste0(mutation_density_global_dir, "/", "III_class_gene_list.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(gene_IV, paste0(mutation_density_global_dir, "/", "IV_class_gene_list.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

    for (type_class in c("SNP", "Indel")) {
        for (gene_color_class in c("I", "II", "III", "IV")) {
            data_for_write <- data_for_plot %>%
                filter(type == type_class) %>%
                filter(gene_color == gene_color_class)
            # cat(type_class,gene_color_class,"的基因个数:",nrow(data_for_write),"\n")
            data_for_write_sorted <- data_for_write[order(data_for_write$peak_proportion, decreasing = T), ]
            write.table(
                data_for_write_sorted$gene_name,
                paste0(mutation_density_global_dir, "/", type_class, "_", gene_color_class, "_class_gene_list.txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
        }
    }
}
# 3. input ---------------------------------------------------------------- TODO:
mutation_density_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/wholeGenome_peak25Perc/3D7"
peak_proportion_global_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
## peak proportion：全基因组
for (gene_length_threshold in c(300, 1000)) {
    for (peak_proportion_threshold in c(0, 0.01, 0.1)) {
        all_p <- plot_different_flank_length(peak_proportion_threshold = peak_proportion_threshold, gene_length_threshold = gene_length_threshold)
        kb0 <- (all_p$core$"0kb"$peak_proportion$core | all_p$core$"0kb"$peak_proportion$noncore | all_p$core$"0kb"$peak_count$core | all_p$core$"0kb"$peak_count$noncore)

        kb1 <- (all_p$core$"1kb"$peak_proportion$core | all_p$core$"1kb"$peak_proportion$noncore | all_p$core$"1kb"$peak_count$core | all_p$core$"1kb"$peak_count$noncore)

        kb2 <- (all_p$core$"2kb"$peak_proportion$core | all_p$core$"2kb"$peak_proportion$noncore | all_p$core$"2kb"$peak_count$core | all_p$core$"2kb"$peak_count$noncore)

        pdf(paste0(mutation_density_global_dir, "/filter_PeakProGreat", peak_proportion_threshold * 100, "Perc_geneLenGreat", gene_length_threshold, "bp.pdf"), width = 30, height = 24)
        p <- (kb0 / kb1 / kb2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
        # grid::grid.draw(grid::textGrob("0kb", x = 0.02, y = 0.97, gp = grid::gpar(fontsize = 50, col = "red")))
        # grid::grid.draw(grid::textGrob("1kb", x = 0.02, y = 0.67, gp = grid::gpar(fontsize = 50, col = "red")))
        # grid::grid.draw(grid::textGrob("2kb", x = 0.02, y = 0.35, gp = grid::gpar(fontsize = 50, col = "red")))
        print(p)
        dev.off()
    }
}
## peak proportion：最好的结果
all_p <- plot_different_flank_length(peak_proportion_threshold = 0, gene_length_threshold = 1000)
(all_p$core$"0kb"$peak_proportion$core | all_p$core$"0kb"$peak_proportion$noncore) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(mutation_density_global_dir, "/filter_PeakProGreat0Perc_geneLenGreat1000bp_peak_proportion.pdf"), width = 18, height = 8)

## outlier genes


## comparing_gene_with_and_without_peak：全基因组
all_p <- list()
for (flank_length in c("0kb", "1kb", "2kb")) {
    data <- read_peak_proportion_and_mutation_density(gene_class_criterion = "core", flank_length = flank_length)
    data_core <- data %>% filter(gene_class == "core")
    data_noncore <- data %>% filter(gene_class == "noncore")
    # whole_p_youwu6mA <- plot_boxplot_youwu6mA(data_for_plot = data, title = "Whole genome")
    core_p_youwu6mA <- plot_boxplot_youwu6mA(
        data_for_plot = data_core,
        title = ifelse(gene_class_criterion == "core", "Core genome", "Non-H3K9me3 region")
    )
    noncore_p_youwu6mA <- plot_boxplot_youwu6mA(
        data_for_plot = data_noncore,
        title = ifelse(gene_class_criterion == "core", "Non-core genome", "H3K9me3 region")
    )
    # core_p_youwu6mA$labels$y <- ""
    noncore_p_youwu6mA$labels$y <- ""

    # whole_p_youwu6mA | core_p_youwu6mA | noncore_p_youwu6mA
    p <- core_p_youwu6mA | noncore_p_youwu6mA
    all_p[[flank_length]] <- p
}

pdf(paste0(mutation_density_global_dir, "/comparing_gene_with_and_without_peak.pdf"), width = 18, height = 24)
(all_p$"0kb" / all_p$"1kb" / all_p$"2kb") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
grid::grid.draw(grid::textGrob("0kb", x = 0.05, y = 0.97, gp = grid::gpar(fontsize = 50, col = "red")))
grid::grid.draw(grid::textGrob("1kb", x = 0.05, y = 0.67, gp = grid::gpar(fontsize = 50, col = "red")))
grid::grid.draw(grid::textGrob("2kb", x = 0.05, y = 0.35, gp = grid::gpar(fontsize = 50, col = "red")))
dev.off()

## var, stevor and rifin
for (pseudoType in c("nonPseudo", "includingPseudo")) {
    gene_list <- ""
    if (pseudoType == "nonPseudo") {
        gene_list <- read_nonPseudo_gene_list()
    } else {
        gene_list <- read_includingPseudo_gene_list()
    }

    all_p <- list()
    for (gene_class_criterion in c("core")) {
        for (flank_length in c("0kb", "1kb", "2kb")) {
            data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
            data_for_plot <- merge(data, gene_list, by = "gene_name")
            peak_propor_p <- plot_var_rifin_stevor_as_one_group(data_for_plot = data_for_plot, x_name = "peak_proportion")
            peak_count_p <- plot_var_rifin_stevor_as_one_group(data_for_plot = data_for_plot, x_name = "peak_count")

            all_p[[gene_class_criterion]][[flank_length]][["peak_propor"]] <- peak_propor_p
            all_p[[gene_class_criterion]][[flank_length]][["peak_count"]] <- peak_count_p
        }
    }

    pdf(paste0(mutation_density_global_dir, "/var_stevor_rifin_", pseudoType, ".pdf"), width = 30, height = 24)

    ((all_p$core$"0kb"$peak_propor | all_p$core$"0kb"$peak_count) /
        (all_p$core$"1kb"$peak_propor | all_p$core$"1kb"$peak_count) /
        (all_p$core$"2kb"$peak_propor | all_p$core$"2kb"$peak_count)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

    grid::grid.draw(grid::textGrob("0kb", x = 0.03, y = 0.97, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("1kb", x = 0.03, y = 0.67, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("2kb", x = 0.03, y = 0.35, gp = grid::gpar(fontsize = 50, col = "red")))
    dev.off()
}

for (pseudoType in c("nonPseudo", "includingPseudo")) {
    gene_list <- ""
    if (pseudoType == "nonPseudo") {
        gene_list <- read_nonPseudo_gene_list()
    } else {
        gene_list <- read_includingPseudo_gene_list()
    }

    all_p <- list()
    for (gene_class_criterion in c("core")) {
        for (flank_length in c("0kb", "1kb", "2kb")) {
            data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
            data_for_plot <- merge(data, gene_list, by = "gene_name")
            peak_propor_p <- plot_var_rifin_stevor_as_many_groups(data_for_plot = data_for_plot, x_name = "peak_proportion")
            peak_count_p <- plot_var_rifin_stevor_as_many_groups(data_for_plot = data_for_plot, x_name = "peak_count")

            all_p[[gene_class_criterion]][[flank_length]][["peak_propor"]] <- peak_propor_p
            all_p[[gene_class_criterion]][[flank_length]][["peak_count"]] <- peak_count_p
        }
    }

    pdf(paste0(mutation_density_global_dir, "/var_stevor_rifin_many_groups", pseudoType, ".pdf"), width = 30, height = 24)

    ((all_p$core$"0kb"$peak_propor | all_p$core$"0kb"$peak_count) /
        (all_p$core$"1kb"$peak_propor | all_p$core$"1kb"$peak_count) /
        (all_p$core$"2kb"$peak_propor | all_p$core$"2kb"$peak_count)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

    grid::grid.draw(grid::textGrob("0kb", x = 0.03, y = 0.97, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("1kb", x = 0.03, y = 0.67, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("2kb", x = 0.03, y = 0.35, gp = grid::gpar(fontsize = 50, col = "red")))
    dev.off()
}


for (pseudoType in c("nonPseudo", "includingPseudo")) {
    gene_list <- ""
    if (pseudoType == "nonPseudo") {
        gene_list <- read_nonPseudo_gene_list()
    } else {
        gene_list <- read_includingPseudo_gene_list()
    }

    all_p <- list()
    for (gene_class_criterion in c("core")) {
        for (flank_length in c("0kb", "1kb", "2kb")) {
            data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
            data_for_plot <- merge(data, gene_list, by = "gene_name")
            peak_propor_p <- plot_var_rifin_stevor_seperately(data_for_plot = data_for_plot, x_name = "peak_proportion")
            peak_count_p <- plot_var_rifin_stevor_seperately(data_for_plot = data_for_plot, x_name = "peak_count")

            all_p[[gene_class_criterion]][[flank_length]][["peak_propor"]] <- peak_propor_p
            all_p[[gene_class_criterion]][[flank_length]][["peak_count"]] <- peak_count_p
        }
    }

    pdf(paste0(mutation_density_global_dir, "/var_stevor_rifin_separately", pseudoType, ".pdf"), width = 36, height = 16)

    for (flank_length in c("0kb", "1kb", "2kb")) {
        p <- (all_p$core[[flank_length]]$peak_propor | all_p$core[[flank_length]]$peak_count) + plot_annotation(title = flank_length)
        print(p)
    }
    dev.off()
}


## var
for (pseudoType in c("nonPseudo", "includingPseudo")) {
    gene_list <- ""
    if (pseudoType == "nonPseudo") {
        gene_list <- read_nonPseudo_gene_list()
    } else {
        gene_list <- read_includingPseudo_gene_list()
    }

    all_p <- list()
    for (gene_class_criterion in c("core")) {
        for (flank_length in c("0kb", "1kb", "2kb")) {
            data <- read_peak_proportion_and_mutation_density(gene_class_criterion = gene_class_criterion, flank_length = flank_length)
            data_for_plot <- merge(data, gene_list, by = "gene_name")
            peak_propor_p <- plot_var_rifin_stevor_as_one_group(data_for_plot = data_for_plot %>% filter(peak_proportion > 0) %>% filter(gene_group == "VAR"), x_name = "peak_proportion")
            peak_count_p <- plot_var_rifin_stevor_as_one_group(data_for_plot = data_for_plot %>% filter(peak_proportion > 0) %>% filter(gene_group == "VAR"), x_name = "peak_count")

            all_p[[gene_class_criterion]][[flank_length]][["peak_propor"]] <- peak_propor_p
            all_p[[gene_class_criterion]][[flank_length]][["peak_count"]] <- peak_count_p
        }
    }

    pdf(paste0(mutation_density_global_dir, "/var_", pseudoType, ".pdf"), width = 16, height = 24)

    ((all_p$core$"0kb"$peak_propor | all_p$core$"0kb"$peak_count) /
        (all_p$core$"1kb"$peak_propor | all_p$core$"1kb"$peak_count) /
        (all_p$core$"2kb"$peak_propor | all_p$core$"2kb"$peak_count)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

    grid::grid.draw(grid::textGrob("0kb", x = 0.03, y = 0.97, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("1kb", x = 0.03, y = 0.67, gp = grid::gpar(fontsize = 50, col = "red")))
    grid::grid.draw(grid::textGrob("2kb", x = 0.03, y = 0.35, gp = grid::gpar(fontsize = 50, col = "red")))
    dev.off()
}


## var thesis
for (pseudoType in c("nonPseudo", "includingPseudo")) {
    gene_list <- ""
    if (pseudoType == "nonPseudo") {
        gene_list <- read_nonPseudo_gene_list()
    } else {
        gene_list <- read_includingPseudo_gene_list()
    }

    data <- read_peak_proportion_and_mutation_density(gene_class_criterion = "core", flank_length = "0kb")
    data_for_plot <- merge(data, gene_list, by = "gene_name")
    data_for_plot_var <- data_for_plot %>% filter(gene_group == "VAR")
    snp <- ggplot(data = data_for_plot_var %>%
        filter(type == "SNP") %>%
        filter(peak_proportion > 0), aes(x = peak_proportion, y = mutation_density)) +
        geom_point(size = 6, color = rgb(189, 63, 61, maxColorValue = 255)) +
        stat_smooth(method = "lm", se = FALSE, color = "black") +
        stat_cor(method = "spearman", label.x.npc = 0.2, label.y.npc = 0.05, size = 12, show.legend = FALSE) +
        labs(x = "Peak proportion", y = "Mutation density") +
        ggtitle("SNP") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )


    indel <- ggplot(data = data_for_plot_var %>%
        filter(type == "Indel") %>%
        filter(peak_proportion > 0), aes(x = peak_proportion, y = mutation_density)) +
        geom_point(size = 6, color = rgb(29, 113, 169, maxColorValue = 255)) +
        stat_smooth(method = "lm", se = FALSE, color = "black") +
        stat_cor(method = "spearman", label.x.npc = 0.2, label.y.npc = 0.05, size = 12, show.legend = FALSE) +
        labs(x = "Peak proportion", y = "Mutation density") +
        ggtitle("Indel") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_text(size = 40, color = "black"),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "bottom",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )

    snp | indel
    ggsave(paste0(mutation_density_global_dir, "/var_", pseudoType, "_thesis.pdf"), width = 16, height = 8)
}
