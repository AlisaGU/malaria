#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)

        data <- data.frame(
            rbind(
                cbind(data[, 1], data[, 2] / data[, 3], "peak"),
                cbind(data[, 1], data[, 4] / data[, 5], "control")
            ), gene_set
        )
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "mutation_density_value", "peak_type", "gene_set_type")
    result <- as.data.frame(result)
    result$mutation_density_value <- as.numeric(as.character(result$mutation_density_value))
    return(result)
}

read_data_with_fold_enrichment <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)

        data <- data.frame(
            rbind(
                cbind(data[, 1], data[, 2] / data[, 3], "peak", data[, 6]),
                cbind(data[, 1], data[, 4] / data[, 5], "control", data[, 6])
            ), gene_set
        )
    })
    result <- do.call(rbind, result)
    colnames(result) <- c("gene_name", "mutation_density_value", "peak_type", "average_peak_fold_enrichment", "gene_set_type")
    result <- as.data.frame(result)
    result$mutation_density_value <- as.numeric(as.character(result$mutation_density_value))
    result$average_peak_fold_enrichment <- as.numeric(as.character(result$average_peak_fold_enrichment))
    result$fold_enrichment_class <- sapply(result$average_peak_fold_enrichment, function(i) {
        if (is.na(i)) {
            return(NA)
        } else {
            shang <- i %/% 1
            yu <- i %% 1
            return(as.character(ifelse(yu < 0.5, shang, shang + 0.5)))
        }
    })
    return(result)
}

summary_data <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)

        peak_mutation_density_scale_by_all_gene_sets <- sum(data[, 2], na.rm = TRUE) / sum(data[, 3], na.rm = TRUE)
        control_mutation_density_scale_by_all_gene_sets <- sum(data[, 4], na.rm = TRUE) / sum(data[, 5], na.rm = TRUE)
        return(c(control_mutation_density_scale_by_all_gene_sets, peak_mutation_density_scale_by_all_gene_sets))
    })
    result <- do.call(rbind, result)
    rownames(result) <- gene_sets
    colnames(result) <- c("control", "peak")
    return(result)
}

gene_mutation_density_direc_with_peak_strength <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)

        data <- data.frame(data[, 1], (data[, 2] / data[, 3]) / (data[, 4] / data[, 5]), data[, 6], gene_set)
        colnames(data) <- c("gene_name", "peak_divided_control", "average_peak_fold_enrichment", "gene_set_type")
        # data$fold_enrichment_class <- sapply(data$average_peak_fold_enrichment, function(i) {
        #     if (is.na(i)) {
        #         return(NA)
        #     } else {
        #         shang <- i %/% 1
        #         yu <- i %% 1
        #         return(as.character(shang + 0.5 * (yu %/% 0.5)))
        #     }
        # })
        data$fold_enrichment_class <- sapply(data$average_peak_fold_enrichment, function(i) {
            if (is.na(i)) {
                return(NA)
            } else {
                if (i <= 2) {
                    return("weak_peak")
                } else if (i > 2 & i <= 2.5) {
                    return("medium_peak")
                } else if (i > 2.5) {
                    return("strong_peak")
                }
            }
        })
        return(data)
    })
    result <- do.call(rbind, result)
    result <- result[!is.na(result$fold_enrichment_class) & !is.infinite(result$peak_divided_control), ]
    ggplot(result, aes(x = fold_enrichment_class, y = peak_divided_control)) +
        geom_boxplot() +
        geom_point() +
        stat_compare_means() +
        scale_y_continuous(limits = c(0, 10)) +
        labs(y = "peak/control")


    aov1 <- aov(peak_divided_control ~ fold_enrichment_class, result)
    summary(aov1)



    a <- table(result$gene_set_type, result$fold_enrichment_class)
    a <- a[, c("weak_peak", "medium_peak", "strong_peak")]
    sapply(gene_sets, function(gene_set) {
        sapply(c("weak_peak", "medium_peak", "strong_peak"), function(fold_enrichment_class) {
            b <- a[gene_set, ]
            conti_table <- matrix(c(b[fold_enrichment_class], sum(b[names(b) != fold_enrichment_class]), sum(a[, fold_enrichment_class]) - b[fold_enrichment_class], sum(a[rownames(a) != gene_set, colnames(a) != fold_enrichment_class])), byrow = T, nrow = 2)
            rownames(conti_table) <- c(paste0("in ", gene_set), paste0("out of ", gene_set))
            colnames(conti_table) <- c(paste0("in ", fold_enrichment_class), paste0("out of ", fold_enrichment_class))
            fisher.test(conti_table, alternative = "greater")$p.value
        })
    })


    ggplot(result %>% filter(gene_set_type != "DR")) +
        geom_density(aes(average_peak_fold_enrichment, group = gene_set_type, color = gene_set_type))



    data2 <- rbind(cbind(result %>% filter(gene_set_type == "STEVOR"), "class" = "STEVOR"), cbind(result %>% filter(gene_set_type != "STEVOR" & gene_set_type != "DR"), "class" = "non-STEVOR"))
    data2$class <- factor(data2$class, levels = c("STEVOR", "non-STEVOR"))
    ggplot(data2, aes(x = class, y = average_peak_fold_enrichment)) +
        geom_point(aes(color = class), position = position_jitterdodge(), size = 8, alpha = 0.5) +
        geom_boxplot(aes(color = class), size = 3, show.legend = FALSE, outlier.shape = NA, alpha = 0.1) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = "left", label.y.npc = 0.85) +
        scale_color_manual(values = c("STEVOR" = "#db161b", "non-STEVOR" = "#1c4e6c")) +
        theme_bw() +
        labs(y = "Gene-associated peaks'\n fold enrichment") +
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
            # legend.position = "bottom",
            legend.position = c(0.35,0.95),
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    ggsave("STEVOR_NONSTEVOR_gene_associated_fold_enrichment_comparison.pdf", width = 8, height = 7)
}

gene_mutation_density_ratio_cor_with_peak_strength <- function() {
    result <- lapply(gene_sets, function(gene_set) {
        data <- read.table(paste0(gene_set, "/snps_mean_variant_count"), header = F, as.is = T)
        data <- data.frame(data[, 1], (data[, 2] / data[, 3]) / (data[, 4] / data[, 5]), data[, 6], gene_set)
        colnames(data) <- c("gene_name", "peak_divided_control", "average_peak_fold_enrichment", "gene_set_type")
        return(data)
    })
    result <- do.call(rbind, result)
    result <- result[!is.na(result$average_peak_fold_enrichment) & !is.infinite(result$peak_divided_control), ]
    p1 <- ggplot(result, aes(x = average_peak_fold_enrichment, y = peak_divided_control)) +
        geom_point(color = "#8d8d8d", fill = "#78addc", size = 8, alpha = 0.5, shape = 21) +
        labs(y = expression(frac('Mutation density'[Peak], 'Mutation density'[Control])), x = "Gene-associated peaks'\n fold enrichment") +
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
    p2 <- p1 + scale_y_continuous(limits = c(0, 10)) +
        stat_smooth(method = "lm", col = "black") +
        stat_cor(method = "spearman", label.x.npc = 0.3, label.y.npc = 0.9, size = 12) +
        theme(axis.title.y = element_blank())

    p1 + p2
    ggsave("mutation_density_ratio_cor_with_peak_strength.pdf", width = 18, height = 7)



    cor.test(result$average_peak_fold_enrichment, result$peak_divided_control, method = "spearman")


    fit1 <- lm(peak_divided_control ~ average_peak_fold_enrichment, data = result)
    summary(fit1)
}



plot_paired <- function() {
    p <- ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value, color = peak_type)) +
        geom_boxplot(aes(group = peak_type), outlier.shape = NA) +
        stat_compare_means(paired = T, size = 12, aes(label = paste0("p = ", ..p.format..)), label.x.npc = 0.6, label.y.npc = 0.88) +
        scale_color_manual(values = c("control" = "grey", "peak" = "red")) +
        facet_wrap(vars(gene_set_type), scales = "free_y") +
        theme_bw() +
        labs(y = "Mutation density") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
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
    p + geom_point(size = 4, position = position_jitterdodge())
    ggsave("gene_5_sets_peak_control_mutation_density_comparison.pdf", width = 21, height = 12)

    p + geom_point(size = 1) + geom_line(aes(group = gene_name), color = "black")
    ggsave("gene_5_sets_peak_control_mutation_density_comparison_direc.pdf", width = 21, height = 12)
}

plot_independent <- function() {
    p <- ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value, color = peak_type)) +
        geom_boxplot(aes(group = peak_type), outlier.shape = NA) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = 0.6, label.y.npc = 0.88) +
        scale_color_manual(values = c("control" = "grey", "peak" = "red")) +
        facet_wrap(vars(gene_set_type), scales = "free_y") +
        theme_bw() +
        labs(y = "Mutation density") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
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
    p + geom_point(size = 4, position = position_jitterdodge())
    ggsave("gene_5_sets_peak_control_mutation_density_comparison_independ.pdf", width = 21, height = 12)

    p + geom_point(size = 1) + geom_line(aes(group = gene_name), color = "black")
    ggsave("gene_5_sets_peak_control_mutation_density_comparison_independ_direc.pdf", width = 21, height = 12)
}

plot_independent_with_fold_enrichment <- function() {
    p <- ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value)) +
        geom_boxplot(aes(color = peak_type, group = peak_type), outlier.shape = NA) +
        stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size = 12, label.x.npc = 0.6, label.y.npc = 0.88) +
        scale_color_manual(values = c(
            "1" = "#deebf7", "1.5" = "#c6dbef",
            "2" = "#9ecae1", "2.5" = "#6baed6", "3" = "#4292c6",
            "3.5" = "#2171b5", "4" = "#08519c", "4.5" = "#08306b", "NA" = "grey", "control" = "grey", "peak" = "red"
        )) +
        facet_wrap(vars(gene_set_type), scales = "free_y") +
        theme_bw() +
        labs(y = "Mutation density") +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            # panel.grid = element_blank(),
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
    p + geom_point(aes(color = fold_enrichment_class), size = 4)
    ggsave("gene_5_sets_peak_control_mutation_density_comparison_independ.pdf", width = 21, height = 12)

    p + geom_point(aes(color = fold_enrichment_class), size = 1) + geom_line(aes(group = gene_name), color = "black")
    ggsave("gene_5_sets_peak_control_mutation_density_comparison_independ_direc.pdf", width = 21, height = 12)
}

shift_legend2 <- function(p) {
    # check if p is a valid object
    if (!(inherits(p, "gtable"))) {
        if (inherits(p, "ggplot")) {
            gp <- ggplotGrob(p) # convert to grob
        } else {
            message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
            return(p)
        }
    } else {
        gp <- p
    }

    # check for unfilled facet panels
    facet.panels <- grep("^panel", gp[["layout"]][["name"]])
    empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]),
        USE.NAMES = F
    )
    empty.facet.panels <- facet.panels[empty.facet.panels]

    if (length(empty.facet.panels) == 0) {
        message("There are no unfilled facet panels to shift legend into. Returning original plot.")
        return(p)
    }

    # establish name of empty panels
    empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
    names <- empty.facet.panels$name

    # return repositioned legend
    reposition_legend(p, "center", panel = names)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_5_sets")
gene_sets <- c("DR", "HDR", "RNA_translation", "STEVOR", "VAR")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read_data()
summary_data()
plot_paired()
plot_independent()



data_for_plot <- read_data_with_fold_enrichment()
p <- ggplot(data = data_for_plot, aes(x = peak_type, y = mutation_density_value)) +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    facet_wrap(vars(gene_set_type), scales = "free_y") +
    geom_boxplot(aes(color = peak_type, group = peak_type), outlier.shape = NA, show.legend = FALSE, lwd = 1.5) +
    geom_point(aes(fill = fold_enrichment_class), size = 8, shape = 21, color = "white", stroke = 10^(-10)) +
    stat_compare_means(comparisons = list(c("control", "peak")), aes(label = paste0("p = ", ..p.format..)), size = 12) +
    scale_x_discrete(labels = c("Control", "Peak")) +
    scale_fill_manual(values = c(
        "1" = "#deebf7", "1.5" = "#c6dbef",
        "2" = "#9ecae1", "2.5" = "#6baed6", "3" = "#4292c6",
        "3.5" = "#2171b5", "4" = "#08519c", "4.5" = "#08306b", "NA" = "grey"
    )) +
    scale_color_manual(values = c(
        "control" = "grey", "peak" = "red"
    )) +
    theme_bw() +
    labs(y = "Mutation density") +
    theme(
        strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
        strip.text.x = element_text(size = 30, color = "white"),
        panel.border = element_blank(),
        # panel.grid = element_blank(),
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
        legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))

pdf("gene_5_sets_peak_control_mutation_density_comparison_independ_p_value_on_the_top.pdf", width = 18, height = 12)
shift_legend2(p)
dev.off()
