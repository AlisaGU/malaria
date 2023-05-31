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
    colnames(chip_info) <- colnames(input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "strand", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")

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
        colnames(data)[2] <- "gene_flank2kb_density_signal"
        return(data)
    }
}

modified <- function(p = NULL) {
    p1 <- p +
        labs(y = "6mA density in gene +-2kb") +
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
WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Tstage"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Tstage"

WT_chip_depth <- 36971023
WT_input_depth <- 21500156

KD_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/chip.gene_flank2kb.inter"
KD_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/input.gene_flank2kb.inter"
KD_chip_depth <- 39189761
KD_input_depth <- 27773718
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
gene_WT_6mA_density <- get_6mA(chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename, chip_depth = WT_chip_depth, input_depth = WT_input_depth)
gene_KD_6mA_density <- get_6mA(chip_info_name = KD_chip_info_filename, input_info_name = KD_input_info_filename, chip_depth = KD_chip_depth, input_depth = KD_input_depth)
colnames(gene_WT_6mA_density)[2] <- "WT"
colnames(gene_KD_6mA_density)[2] <- "KD"

gene_WT_6mA_density <- gene_WT_6mA_density[order(gene_WT_6mA_density$WT, decreasing = FALSE), ]
gene_WT_6mA_density$WT_level <- ceiling(1:nrow(gene_WT_6mA_density) / (nrow(gene_WT_6mA_density) / 5)) # level 1的6mA最低


gene_density <- gene_WT_6mA_density %>%
    left_join(gene_KD_6mA_density, by = "gene") %>%
    na.omit()

gene_density$WT_level <- as.factor(gene_density$WT_level)

gene_density$WTminusKD <- gene_density$WT - gene_density$KD
gene_density$WTdividedKD <- gene_density$WT / gene_density$KD

wt <- ggplot(data = gene_density, aes(x = WT_level, y = WT)) +
    geom_boxplot(color = "black", size = 1, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 3))

kd <- ggplot(data = gene_density, aes(x = WT_level, y = KD)) +
    geom_boxplot(color = "red", size = 1, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 3))

difference <- ggplot(data = gene_density, aes(x = WT_level, y = WTminusKD)) +
    geom_boxplot(color = "#2a562c", size = 1, outlier.shape = NA) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-0.3, 0.4))

divided <- ggplot(data = gene_density, aes(x = WT_level, y = WTdividedKD)) +
    geom_boxplot(color = "#d48c31", size = 1, outlier.shape = NA) +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(0.7, 1.5))

wt_modi <- modified(wt)
kd_modi <- modified(kd)
diff_modi <- modified(difference)
divi_modi <- modified(divided)

(wt_modi / kd_modi) | (diff_modi / divi_modi)
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/difference between WT and KD_3group.pdf", width = 10, height = 8)

wt_modi / kd_modi / diff_modi
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/difference between WT and KD_3group_1.pdf", width = 8, height = 15)



summary_WT <- gene_density %>%
    group_by(WT_level) %>%
    summarise(
        data = median(WT)
    )
summary_KD <- gene_density %>%
    group_by(WT_level) %>%
    summarise(
        data = median(KD)
    )
summary_WTminusKD <- gene_density %>%
    group_by(WT_level) %>%
    summarise(
        data = median(WTminusKD)
    )
summary_WTdividedKD <- gene_density %>%
    group_by(WT_level) %>%
    summarise(
        data = median(WTdividedKD)
    )
summary_info <- data.frame(rbind(summary_WT, summary_KD, summary_WTminusKD, summary_WTdividedKD),
    type = c(rep("WT", nrow(summary_WT)), rep("KD", nrow(summary_KD)), rep("WTminusKD", nrow(summary_WTminusKD)), rep("WTdividedKD", nrow(summary_WTdividedKD)))
)
summary_info$type <- factor(summary_info$type)
ggplot(summary_info %>% filter(type == "WT" | type == "KD"), aes(x = WT_level, y = data, group = type, color = type)) +
    geom_line() +
    geom_point()
