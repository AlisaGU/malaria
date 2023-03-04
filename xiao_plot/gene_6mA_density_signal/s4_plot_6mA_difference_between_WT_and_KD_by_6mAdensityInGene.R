#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(ggbreak)))

# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(chip_info_name = NULL, input_info_name = NULL, chip_depth = NULL, input_depth = NULL) {
    chip_info <- fread(chip_info_name, stringsAsFactors = F, header = F)
    input_info <- fread(input_info_name, stringsAsFactors = F, header = F)
    colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

    chip_info <- chip_info %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")
    input_info <- input_info %>% filter(chrom != "Pf_M76611" & chrom != "Pf3D7_API_v3")

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
        colnames(data)[2] <- "gene_density_signal"
        return(data)
    }
}

modified <- function(p = NULL) {
    p1 <- p +
        labs(y = "6mA density in gene") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 24, color = "black"),
            axis.text.y = element_text(size = 24, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            ),
            axis.title.y = element_text(size = 20),
            legend.position = "none"
        )
    return(p1)
}

# 3. input ---------------------------------------------------------------- TODO:
WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"

WT_chip_depth <- 36971023
WT_input_depth <- 21500156

KD_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/chip.inter"
KD_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/input.inter"
KD_chip_depth <- 39189761
KD_input_depth <- 27773718
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
gene_WT_6mA_density <- get_6mA(chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename, chip_depth = WT_chip_depth, input_depth = WT_input_depth)
gene_KD_6mA_density <- get_6mA(chip_info_name = KD_chip_info_filename, input_info_name = KD_input_info_filename, chip_depth = KD_chip_depth, input_depth = KD_input_depth)
colnames(gene_WT_6mA_density)[2] <- "WT"
colnames(gene_KD_6mA_density)[2] <- "KD"

gene_WT_6mA_density <- gene_WT_6mA_density[order(gene_WT_6mA_density$WT, decreasing = FALSE), ]
gene_WT_6mA_density$WT_level <- ceiling(1:nrow(gene_WT_6mA_density) / (nrow(gene_WT_6mA_density) / 3)) # level 1的6mA最低


gene_density <- gene_WT_6mA_density %>%
    left_join(gene_KD_6mA_density, by = "gene") %>%
    na.omit()
infinite_index <- unique(c(which(unlist(gene_density$WT) == Inf), which(unlist(gene_density$KD) == Inf)))
gene_density <- gene_density[-infinite_index, ]
gene_density$WT_level <- as.factor(gene_density$WT_level)

gene_density$WTminusKD <- gene_density$WT - gene_density$KD
gene_density$WTdividedKD <- gene_density$WT / gene_density$KD

wt <- ggplot(data = gene_density, aes(x = WT_level, y = WT)) +
    geom_boxplot(color = "black", size = 1, outlier.shape = 1, outlier.size = 0.5) +
    scale_y_break(c(6, 15), scales = 0.1) +
    scale_y_continuous(
        breaks = c(0, 2, 4, 6, 15, 18),
        labels = c(0, 2, 4, 6, 15, 18)
    )

kd <- ggplot(data = gene_density, aes(x = WT_level, y = KD)) +
    geom_boxplot(color = "red", size = 1, outlier.shape = 1, outlier.size = 0.5) +
    scale_y_break(c(5, 15), scales = 0.1) +
    scale_y_continuous(
        breaks = c(0, 2, 4, 15, 18),
        labels = c(0, 2, 4, 15, 18)
    )

difference <- ggplot(data = gene_density, aes(x = WT_level, y = WTminusKD)) +
    geom_boxplot(color = "#2a562c", size = 1, outlier.shape = 1, outlier.size = 0.5) +
    scale_y_break(c(3, 15), scales = 1) +
    scale_y_break(c(-7, -0.5), scales = 15) +
    scale_y_continuous(
        breaks = c(-7, 0, 2, 4, 15, 16),
        labels = c(-7, 0, 2, 4, 15, 16)
    ) +
    geom_hline(yintercept = 0)

divided <- ggplot(data = gene_density, aes(x = WT_level, y = WTdividedKD)) +
    geom_boxplot(color = "#d48c31", size = 1, outlier.shape = 1, outlier.size = 0.5) +
    scale_y_break(c(3, 15), scales = 0.1) +
    scale_y_continuous(
        breaks = c(0, 2, 4, 15, 20),
        labels = c(0, 2, 4, 15, 20)
    ) +
    geom_hline(yintercept = 1)

wt_modi <- modified(wt)
kd_modi <- modified(kd)
diff_modi <- modified(difference)
divi_modi <- modified(divided)

(wt_modi / kd_modi) | (diff_modi / divi_modi)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/difference between WT and KD_3group.pdf", width = 10, height = 5)
