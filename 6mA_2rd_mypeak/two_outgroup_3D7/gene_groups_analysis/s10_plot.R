#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
# 2. functions ------------------------------------------------------------ TODO:
read_data <- function() {
    result <- list()
    for (pseudo in c("exclude_pseudo", "include_pseudo")) {
        a <- list()
        for (gene_group in gene_group_type) {
            data <- read.table(paste0(pseudo, "/", gene_group, "/snps_mean_variant_count"), header = F, as.is = T)
            mutation_count_index <- seq(1, ncol(data), by = 2)
            window_len_index <- seq(2, ncol(data), by = 2)
            snp_mutation_density <- apply(data, 1, function(x) {
                sum(x[mutation_count_index]) / sum(x[window_len_index])
            })

            data <- read.table(paste0(pseudo, "/", gene_group, "/indels_mean_variant_count"), header = F, as.is = T)
            mutation_count_index <- seq(1, ncol(data), by = 2)
            window_len_index <- seq(2, ncol(data), by = 2)
            indel_mutation_density <- apply(data, 1, function(x) {
                sum(x[mutation_count_index]) / sum(x[window_len_index])
            })

            data <- read.table(paste0(
                "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/",
                pseudo, "/", gene_group, "/", gene_group, "_window50_6mA_FC"
            ))
            # m6mA <- apply(data, 1, mean)
            m6mA <- apply(data, 1, median)

            a[[gene_group]] <- data.frame(
                window = paste("window", 1:50, sep = ""), snp = snp_mutation_density,
                indel = indel_mutation_density, m6mA = m6mA, pseudo = pseudo, gene_groups = gene_group
            )
        }
        a <- do.call(rbind, a)
        result[[pseudo]] <- a
    }
    result <- do.call(rbind, result)
    result$gene_groups <- factor(result$gene_groups, levels = gene_group_type)
    return(result)
}

plot_all_genes_scale_by_fixed_factor_fixed_cordi <- function() {
    ratio <- 5

    pdf("all_genes_scale_by_fixed_factor_fixed_cordi.pdf", width = 30, height = 30)
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        p <- ggplot(data = data_for_plot %>% filter(pseudo == pseudo_type)) +
            # geom_point(aes(x = window, y = snp,color="SNPs"),show.legend = T,size=4) +
            geom_line(aes(x = window, y = snp, color = "SNPs", group = 1), size = 3) +
            # geom_point(aes(x = window, y = indel,color="Indels"),size=4) +
            geom_line(aes(x = window, y = indel, color = "Indels", group = 1), size = 3) +
            # geom_point(aes(x = window, y = m6mA / ratio,color="m6mA"),size=4) +
            geom_line(aes(x = window, y = m6mA / ratio, color = "6mA", group = 1), size = 3) +
            # scale_color_manual(values=c(SNPs="red",Indels="blue",m6mA="green"))+
            scale_y_continuous(
                name = "Mutation density",
                sec.axis = sec_axis(~ . * ratio, name = "6mA density fold change")
            ) +
            ggtitle(pseudo_type) +
            scale_x_discrete(labels = c(rep("", 4), 5, rep("", 4), 10, rep("", 4), 15, rep("", 4), 20, rep("", 4), 25, rep("", 4), 30, rep("", 4), 35, rep("", 4), 40, rep("", 4), 45, rep("", 4), 50)) +
            facet_wrap(~gene_groups, nrow = 5, ncol = 3) +
            theme_bw() +
            labs(y = "Mutation density") +
            scale_color_manual(values = c("red", "blue", "green")) +
            theme(
                strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
                strip.text.x = element_text(size = 32, color = "white"),
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
                legend.text = element_text(size = 50, color = "black"),
                legend.position = c(0.87, 0.06),
                plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
                panel.spacing = unit(3, "lines")
            )
        print(p)
    }
    dev.off()
}

plot_all_genes_scale_by_fixed_factor_free_cordi <- function() {
    ratio <- 5

    pdf("all_genes_scale_by_fixed_factor_free_cordi.pdf", width = 30, height = 30)
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        p <- ggplot(data = data_for_plot %>% filter(pseudo == pseudo_type)) +
            geom_line(aes(x = window, y = snp, color = "SNPs", group = 1), size = 3) +
            geom_line(aes(x = window, y = indel, color = "Indels", group = 1), size = 3) +
            geom_line(aes(x = window, y = m6mA / ratio, color = "6mA", group = 1), size = 3) +
            scale_y_continuous(
                name = "Mutation density",
                sec.axis = sec_axis(~ . * ratio, name = "6mA density fold change")
            ) +
            ggtitle(pseudo_type) +
            scale_x_discrete(labels = c(rep("", 4), 5, rep("", 4), 10, rep("", 4), 15, rep("", 4), 20, rep("", 4), 25, rep("", 4), 30, rep("", 4), 35, rep("", 4), 40, rep("", 4), 45, rep("", 4), 50)) +
            facet_wrap(~gene_groups, nrow = 5, ncol = 3, scales = "free_y") +
            theme_bw() +
            labs(y = "Mutation density") +
            scale_color_manual(values = c("red", "blue", "green")) +
            theme(
                strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
                strip.text.x = element_text(size = 32, color = "white"),
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
                legend.text = element_text(size = 50, color = "black"),
                legend.position = c(0.87, 0.06),
                plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
                panel.spacing = unit(3, "lines")
            )
        print(p)
    }
    dev.off()
}

plot_all_genes_scale_by_free_factor_free_cordi <- function() {
    pdf("all_genes_scale_by_free_factor_free_cordi.pdf", width = 42, height = 32)
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        result_plot <- list()
        for (gene_group in gene_group_type) {
            data <- data_for_plot %>%
                filter(pseudo == pseudo_type) %>%
                filter(gene_groups == gene_group)
            # ratio <- max(c(max(data$snp / data$m6mA), max(data$indel / data$m6mA)))
            data1 <- data.frame(
                mutation_density = c(data$snp, data$indel),
                m6ma = c(data$m6mA, data$m6mA)
            )
            y1 <- min(data1$mutation_density)
            y2 <- max(data1$mutation_density)
            x1 <- min(data1$m6ma)
            x2 <- max(data1$m6ma)
            b <- (y2 - y1) / (x2 - x1)
            a <- y1 - b * x1


            result_plot[[as.character(gene_group)]] <- plot_sub(data = data, a = a, b = b, gene_group = gene_group)
        }
        p1 <- (result_plot$BER / result_plot$HDR / result_plot$MMR / result_plot$NER / result_plot$invasion) | (result_plot$rhoptry / result_plot$RIF / result_plot$STEVOR / result_plot$SURF / result_plot$VAR) | (result_plot$PHIST / result_plot[["MC-2TM"]] / result_plot$DnaJ / result_plot$DR / plot_spacer())
        print(p1)
    }
    dev.off()
}

plot_sub <- function(data = NULL, a = NULL, b = NULL, gene_group = NULL) {
    ppp <- ggplot(data = data) +
        geom_line(aes(x = window, y = snp, color = "SNPs", group = 1), size = 3) +
        geom_line(aes(x = window, y = indel, color = "Indels", group = 1), size = 3) +
        # geom_line(aes(x = window, y = m6mA * ratio, color = "6mA", group = 1), size = 3) +
        geom_line(aes(x = window, y = m6mA * b + a, color = "6mA", group = 1), data = data, size = 3) +
        geom_hline(yintercept = a + b, linetype = "dashed", size = 2.5) +
        scale_y_continuous(
            name = "Mutation density",
            # sec.axis = sec_axis(~ . / ratio, name = "6mA density fold change")
            # sec.axis = sec_axis(~ (. - a) / b, name = "6mA density fold change")
            sec.axis = sec_axis(~ (. - a) / b, name = " ")
        ) +
        ggtitle(gene_group) +
        scale_x_discrete(labels = c(
            rep("", 4), 5, rep("", 4), 10, rep("", 4), 15,
            rep("", 4), 20, rep("", 4), 25, rep("", 4), 30, rep("", 4), 35,
            rep("", 4), 40, rep("", 4), 45, rep("", 4), 50
        )) +
        theme_bw() +
        labs(y = "Mutation density") +
        scale_color_manual(values = c("red", "blue", "green")) +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 32, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 50, color = "black"),
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines"),
            legend.position = "none"
        )
    return(ppp)
}

plot_all_genes_6ma_vs_md <- function() {
    pdf("all_genes_6ma_vs_md.pdf", width = 42, height = 32)
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        result_plot <- list()
        for (gene_group in gene_group_type) {
            data <- data_for_plot %>%
                filter(pseudo == pseudo_type) %>%
                filter(gene_groups == gene_group)
            data1 <- data.frame(
                mutation_density = c(data$snp, data$indel),
                m6ma = c(data$m6mA, data$m6mA),
                variant = c(rep("SNPs", nrow(data)), rep("Indels", nrow(data)))
            )
            ggplot(data = data1) +
                geom_point(aes(x = mutation_density, y = m6ma, color = variant)) +
                stat_smooth(aes(x = mutation_density, y = m6ma, color = variant, group = variant), method = "lm", col = "black")
        }
        p1 <- (result_plot$BER / result_plot$HDR / result_plot$MMR / result_plot$NER / result_plot$invasion) | (result_plot$rhoptry / result_plot$RIF / result_plot$STEVOR / result_plot$SURF / result_plot$VAR) | (result_plot$PHIST / result_plot[["MC-2TM"]] / result_plot$DnaJ / result_plot$DR / plot_spacer())
        print(p1)
    }
    dev.off()
}

plot_each_gene_group <- function() {
    pdf("all_genes_each_window.pdf", width = 10, height = 5)
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        for (gene_group in gene_group_type) {
            data <- read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/", pseudo_type, "/", gene_group, "/", gene_group, "_window50_6mA_FC"), header = F, as.is = T)
            gene_name <- read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/", pseudo_type, "/", gene_group, "/", gene_group, "_genes.bed"), header = F, as.is = T)
            gene_name <- unlist(gene_name[, 4])
            data1 <- data.frame(value = unlist(data), window = rep(1:nrow(data), ncol(data)), gene = rep(gene_name, each = nrow(data)))
            data1$window <- factor(data1$window, levels = 1:nrow(data))
            data1$gene <- factor(data1$gene, levels = gene_name)
            p <- ggplot(data = data1) +
                geom_point(aes(x = window, y = value, color = gene)) +
                geom_line(aes(x = window, y = value, group = gene, color = gene)) +
                geom_boxplot(aes(x = window, y = value)) +
                theme_bw() +
                ggtitle(paste0(pseudo_type, "/", gene_group)) +
                theme(
                    legend.title = element_blank(),
                    legend.text = element_text(size = 4, color = "black"),
                    legend.key.size = unit(0.5, "line")
                )
            print(p)
        }
    }
    dev.off()
}

plot_each_gene_6mA_mutation_density <- function() {
    for (pseudo_type in c("exclude_pseudo", "include_pseudo")) {
        for (gene_group in gene_group_type) {
            gene_list <- unlist(read.table(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/", pseudo_type, "/", gene_group, "/", gene_group, "_genes.bed"), header = F, as.is = T)[, 4])

            data_snp <- read.table(paste0(pseudo_type, "/", gene_group, "/snps_mean_variant_count"), header = F, as.is = T)
            data_indel <- read.table(paste0(pseudo_type, "/", gene_group, "/indels_mean_variant_count"), header = F, as.is = T)
            data_6mA <- read.table(paste0(
                "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/",
                pseudo_type, "/", gene_group, "/", gene_group, "_window50_6mA_FC"
            ))

            pdf(paste0(pseudo_type, "_", gene_group, ".pdf"), width = 10, height = 7)
            for (gene_order in 1:length(gene_list)) {
                data <- data.frame(
                    window = 1:50, m6mA = data_6mA[, gene_order],
                    snp = data_snp[, 2 * gene_order - 1] / data_snp[, 2 * gene_order],
                    indel = data_indel[, 2 * gene_order - 1] / data_snp[, 2 * gene_order]
                )
                data$window <- factor(data$window, levels = 1:50)
                gene <- gene_list[gene_order]
                data1 <- data.frame(
                    mutation_density = c(data$snp, data$indel),
                    m6ma = c(data$m6mA, data$m6mA)
                )
                y1 <- min(data1$mutation_density)
                y2 <- max(data1$mutation_density)
                x1 <- min(data1$m6ma)
                x2 <- max(data1$m6ma)
                b <- (y2 - y1) / (x2 - x1)
                a <- y1 - b * x1
                p <- plot_sub(
                    data = data, a = a, b = b,
                    gene_group = gene
                )
                print(p)
            }
            dev.off()
        }
    }
}

read_peak_window_data <- function() {
    result <- list()
    for (pseudo in c("exclude_pseudo", "include_pseudo")) {
        a <- list()
        for (gene_group in gene_group_type) {
            data <- read.table(paste0(pseudo, "/", gene_group, "/snps_mean_variant_count_peak_window_50bp"), header = F, as.is = T)
            snp_mutation_density <- data.frame(snp_count = data[, 1], snp_length = data[, 2], snp_density = data[, 1] / data[, 2])

            data <- read.table(paste0(pseudo, "/", gene_group, "/indels_mean_variant_count_peak_window_50bp"), header = F, as.is = T)
            indel_mutation_density <- data.frame(indel_count = data[, 1], indel_length = data[, 2], indel_density = data[, 1] / data[, 2])

            data <- read.table(paste0(
                "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/",
                # pseudo, "/", gene_group, "/", gene_group, "_50bp_6mA_FC"
                pseudo, "/", gene_group, "/", gene_group, "_50bp_bam_6mA_FC"
            ))
            colnames(data) <- "m6mA"
            a[[gene_group]] <- data.frame(part = 1:nrow(data), snp_mutation_density, indel_mutation_density, m6mA = data, pseudo = pseudo, gene_groups = gene_group)

            # belong_gene_data <- read.table(paste0(
            #     "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/",
            #     pseudo, "/", gene_group, "/", gene_group, "_genes.peak.window.50bp.anno.bed"
            # ))
            # belong_gene_data$V5 <- as.factor(belong_gene_data$V5)
            # belong_gene <- t(sapply(split(belong_gene_data, f = belong_gene_data$V5), function(x) {
            #     b <- unlist(strsplit(x[1, 4], "_"))
            #     return(c(paste(b[1], b[2], sep = "_"), x[1, 6]))
            # }))
            # a[[gene_group]] <- data.frame(part = 1:nrow(data), snp_mutation_density, indel_mutation_density, m6mA = data, pseudo = pseudo, gene_groups = gene_group, belong_gene = belong_gene[,1],region=belong_gene[,2])
        }
        a <- do.call(rbind, a)
        result[[pseudo]] <- a
    }
    result <- do.call(rbind, result)
    result$gene_groups <- factor(result$gene_groups, levels = gene_group_type)
    return(result)
}

plot_gene_groups_scatter <- function(data = NULL, method = NULL) {
    # pdf("peak_window_snp.pdf",width=30,height=30)
    pdf(paste0("peak_window_snp_", method, ".pdf"), width = 30, height = 40)
    for (snp_count_threshold in 0) {
        p <- plot_each_one(data = data %>%
            filter(pseudo == "exclude_pseudo") %>%
            filter(snp_count >= snp_count_threshold), title = paste0("exclude_pseudo_let", snp_count_threshold), method = method)
        print(p)
    }
    dev.off()
}

plot_gene_groups_scatter_limit_region_length <- function(data = NULL, method = NULL) {
    # pdf("peak_window_snp.pdf",width=30,height=30)
    pdf(paste0("peak_window_snp_region_50_", method, ".pdf"), width = 30, height = 40)
    for (snp_count_threshold in 0) {
        p <- plot_each_one(
            data = data %>%
                filter(pseudo == "exclude_pseudo") %>%
                filter(snp_length == 50) %>%
                filter(snp_count >= snp_count_threshold),
            title = paste0("exclude_pseudo_let", snp_count_threshold), method = method
        )
        print(p)
    }
    dev.off()
}

plot_gene_groups_scatter_limit_region_length_count <- function(data = NULL, method = NULL) {
    # pdf("peak_window_snp.pdf",width=30,height=30)
    pdf(paste0("peak_window_snp_region_50_count_", method, ".pdf"), width = 30, height = 40)
    for (snp_count_threshold in 0:5) {
        p <- plot_each_one_count(
            data = data %>%
                filter(pseudo == "exclude_pseudo") %>%
                filter(snp_length == 50) %>%
                filter(snp_count >= snp_count_threshold),
            title = paste0("exclude_pseudo_let", snp_count_threshold), method = method
        )
        print(p)
    }
    dev.off()
}

plot_each_one <- function(data = NULL, title = NULL, method = NULL) {
    p <- ggplot(data = data) +
        geom_point(aes(x = m6mA, y = snp_density), size = 4) +
        # geom_point(aes(x = m6mA, y = snp_density, color = belong_gene), size = 4) +
        geom_smooth(aes(x = m6mA, y = snp_density),
            method = "lm",
            se = TRUE, fullrange = FALSE, level = 0.95
        ) +
        facet_wrap(~gene_groups, nrow = 5, ncol = 3, scales = "free") +
        stat_cor(aes(x = m6mA, y = snp_density), method = method, size = 20) +
        ggtitle(title) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 32, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 50, color = "black"),
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines"),
            legend.position = "none"
        )
    return(p)
}

plot_each_one_count <- function(data = NULL, title = NULL, method = NULL) {
    p <- ggplot(data = data) +
        geom_point(aes(x = m6mA, y = snp_count), size = 4) +
        geom_smooth(aes(x = m6mA, y = snp_count),
            method = "auto",
            se = TRUE, fullrange = FALSE, level = 0.95
        ) +
        facet_wrap(~gene_groups, nrow = 5, ncol = 3, scales = "free") +
        stat_cor(aes(x = m6mA, y = snp_count), method = method, size = 20) +
        ggtitle(title) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 32, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 36, color = "black"),
            axis.text.y = element_text(size = 36, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 40, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 50, color = "black"),
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines"),
            legend.position = "none"
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/two_outgroup_consistent_in_core_and_noncore/gene_groups")
gene_group_type <- c("BER", "HDR", "MMR", "NER", "invasion", "rhoptry", "RIF", "STEVOR", "SURF", "VAR", "PHIST", "MC-2TM", "DnaJ", "DR")
# gene_group_type <- c("BER", "HDR", "MMR", "NER", "rhoptry", "RIF", "STEVOR", "SURF", "PHIST", "MC-2TM", "DnaJ", "DR")

# gene_group_type <- c("BER", "HDR", "MMR", "NER", "invasion", "rhoptry", "STEVOR", "SURF", "VAR", "PHIST", "MC-2TM", "DnaJ", "DR")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- read_data()
data_for_plot$window <- factor(data_for_plot$window, levels = paste("window", 1:50, sep = ""))
plot_all_genes_scale_by_fixed_factor_fixed_cordi()
# plot_all_genes_scale_by_fixed_factor_free_cordi()
plot_all_genes_scale_by_free_factor_free_cordi()

plot_each_gene_group()
plot_each_gene_6mA_mutation_density()


data_for_plot_peak_window <- read_peak_window_data()
plot_gene_groups_scatter(data = data_for_plot_peak_window, method = "pearson")
plot_gene_groups_scatter(data = data_for_plot_peak_window, method = "spearman")
plot_gene_groups_scatter(data = data_for_plot_peak_window, method = "kendall")
plot_gene_groups_scatter_limit_region_length(data = data_for_plot_peak_window, method = "pearson")
plot_gene_groups_scatter_limit_region_length(data = data_for_plot_peak_window, method = "spearman")
plot_gene_groups_scatter_limit_region_length(data = data_for_plot_peak_window, method = "kendall")

plot_gene_groups_scatter_limit_region_length_count(data = data_for_plot_peak_window, method = "pearson")
plot_gene_groups_scatter_limit_region_length_count(data = data_for_plot_peak_window, method = "spearman")
plot_gene_groups_scatter_limit_region_length_count(data = data_for_plot_peak_window, method = "kendall")

data <- data_for_plot_peak_window %>%
    filter(pseudo == "exclude_pseudo") %>%
    filter(gene_groups == "VAR") %>%
    filter(belong_gene == "PF3D7_0100100") %>%
    filter(snp_count > 0)

ggplot(data = data %>% filter(region == "CDS")) +
    geom_point(aes(x = m6mA, y = snp_density, color = region), size = 8) +
    geom_smooth(aes(x = m6mA, y = snp_density),
        method = "lm",
        se = TRUE, fullrange = FALSE, level = 0.95
    ) +
    facet_wrap(~belong_gene, scales = "free") +
    stat_cor(aes(x = m6mA, y = snp_density), method = "spearman", size = 25)

ggsave("a.pdf", width = 80, height = 80, limitsize = FALSE)
ggsave("b.pdf", width = 80, height = 80, limitsize = FALSE)



ggplot(data = data) +
    geom_point(aes(x = m6mA, y = snp_density), size = 8) +
    geom_smooth(aes(x = m6mA, y = snp_density),
        method = "lm",
        se = TRUE, fullrange = FALSE, level = 0.95
    ) +
    stat_cor(aes(x = m6mA, y = snp_density), method = "spearman", size = 15)

fit <- lm(data$snp_density ~ data$m6mA)
summary(fit)
