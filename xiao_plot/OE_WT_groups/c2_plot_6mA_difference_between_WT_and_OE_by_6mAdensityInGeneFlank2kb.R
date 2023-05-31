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
WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/3D7-T_chip.gene_flank2kb.inter"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/3D7-T_input.gene_flank2kb.inter"

WT_chip_depth <- 35068263
WT_input_depth <- 27989115

KD_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/NAMTOE-T_chip.gene_flank2kb.inter"
KD_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/NAMTOE-T_input.gene_flank2kb.inter"
KD_chip_depth <- 32038632
KD_input_depth <- 8218045
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
gene_WT_6mA_density <- get_6mA(chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename, chip_depth = WT_chip_depth, input_depth = WT_input_depth)
gene_KD_6mA_density <- get_6mA(chip_info_name = KD_chip_info_filename, input_info_name = KD_input_info_filename, chip_depth = KD_chip_depth, input_depth = KD_input_depth)
colnames(gene_WT_6mA_density)[2] <- "WT"
colnames(gene_KD_6mA_density)[2] <- "OE"

gene_WT_6mA_density <- gene_WT_6mA_density[order(gene_WT_6mA_density$WT, decreasing = FALSE), ]

## 5组
gene_WT_6mA_density$WT_level <- ceiling(1:nrow(gene_WT_6mA_density) / (nrow(gene_WT_6mA_density) / 5)) # level 1的6mA最低

write.table(gene_WT_6mA_density$gene[gene_WT_6mA_density$WT_level == 1], "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/class1.geneList.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(gene_WT_6mA_density$gene[gene_WT_6mA_density$WT_level == 5], "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/class5.geneList.txt", row.names = F, col.names = F, sep = "\t", quote = F)

gene_density <- gene_WT_6mA_density %>%
    left_join(gene_KD_6mA_density, by = "gene") %>%
    na.omit()

gene_density$WT_level <- as.factor(gene_density$WT_level)

gene_density$OEminusWT <- gene_density$OE - gene_density$WT
gene_density$OEdividedWT <- gene_density$OE / gene_density$WT

wt <- ggplot(data = gene_density, aes(x = WT_level, y = WT)) +
    geom_boxplot(color = "black", size = 1, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 2.5))

kd <- ggplot(data = gene_density, aes(x = WT_level, y = OE)) +
    geom_boxplot(color = "red", size = 1, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 2.5))

difference <- ggplot(data = gene_density, aes(x = WT_level, y = OEminusWT)) +
    geom_boxplot(color = "#2a562c", size = 1, outlier.shape = NA) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-0.2, 0.3))


wt_modi <- modified(wt)
kd_modi <- modified(kd)
diff_modi <- modified(difference)


wt_modi / kd_modi / diff_modi
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_WT_groups/difference between WT and OE_5_groups.pdf", width = 8, height = 15)

## more groups
gene_WT_6mA_density$WT_level <- ceiling(1:nrow(gene_WT_6mA_density) / (nrow(gene_WT_6mA_density) / 5)) # level 1的6mA最低


gene_density <- gene_WT_6mA_density %>%
    left_join(gene_KD_6mA_density, by = "gene") %>%
    na.omit()

gene_density$WT_level <- as.factor(gene_density$WT_level)

gene_density$OEminusWT <- gene_density$OE - gene_density$WT

data_for_plot <- do.call(rbind, list(
    data.frame(gene = gene_density[, 1], value = gene_density[, 2], WT_level = gene_density[, 3], group = "WT"),
    data.frame(gene = gene_density[, 1], value = gene_density[, 4], WT_level = gene_density[, 3], group = "OE"),
    data.frame(gene = gene_density[, 1], value = gene_density[, 5], WT_level = gene_density[, 3], group = "OEminusWT")
))

data_for_plot1 <- data_for_plot %>%
    filter(!is.infinite(value)) %>%
    group_by(group, WT_level) %>%
    summarise(median = median(value))
ggplot(data_for_plot1 %>% filter(group != "OEminusWT"), aes(x = WT_level, y = median, group = group, color = group)) +
    geom_point() +
    geom_line()
