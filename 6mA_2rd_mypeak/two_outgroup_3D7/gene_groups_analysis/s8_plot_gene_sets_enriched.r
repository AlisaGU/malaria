#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(dplyr)
library(ggrepel)
# 2. functions ------------------------------------------------------------ TODO:
plot_gene_sets <- function(data = NULL) {
    p <- ggplot(data = data, aes(x = gene_region, y = logp, group = gene_sets, color = color, label = gene_sets)) +
        geom_point(
            color = dplyr::case_when(
                data$enrichment_pvalue >= 0.05 ~ "grey",
                data$enrichment_pvalue < 0.05 ~ "blue"
            ), size = 4, alpha = 0.8
        ) +
        geom_hline(yintercept = -log(0.05, base = 10), color = "red", linetype = "dashed", size = 1.5) +
        geom_text_repel(aes(label = ifelse(enrichment_pvalue < 0.05, as.character(gene_sets), "")),
            box.padding = 0.35, size = 7,
            point.padding = 0.5, max.overlaps = 12,
            segment.color = "grey50", direction = "x"
        ) +
        labs(y = "-log(pvalue)") +
        scale_x_discrete(labels = c("all peaks\n(exon)", "top 500 peaks\n(exon)", "top500 peaks\nsubmit500\n(exon)", "all peaks\n(ATG 2kb)", "top 500 peaks\n(ATG 2kb)", "top500 peaks\nsubmit500\n(ATG 2kb)")) +
        facet_wrap(~pseudo, nrow = 2, ncol = 1) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 30, color = "black", angle = 45, hjust = 0.7, vjust = 0.7),
            axis.text.x = element_text(size = 24, color = "black"),
            axis.text.y = element_text(size = 24, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 24, color = "black"),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}

# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
data_for_plot_name <- Args[6]

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read.table(data_for_plot_name, as.is = T, header = F, sep = " ")
colnames(data_for_plot) <- c(
    "pseudo", "gene_region", "gene_sets", "overlap_count", "peak_related_gene_count",
    "whole_genome_gene_count_except_peak_related", "gene_group_count", "enrichment_pvalue"
)
data_for_plot$logp <- -log(data_for_plot$enrichment_pvalue, base = 10)
data_for_plot$logp[88] <- 10 # FIXME:

data_for_plot$color <- data_for_plot$enrichment_pvalue < 0.05
data_for_plot$pseudo <- factor(data_for_plot$pseudo, levels = c("include_pseudo", "exclude_pseudo"))
plot_gene_sets(data = data_for_plot)
ggsave(paste0(dirname(data_for_plot_name), "/significant_gene_sets.pdf"), width = 15, height = 15)
