#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggpubr)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
get_6mA <- function(chip_info = NULL, input_info = NULL, chip_depth = NULL, input_depth = NULL) {
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
        return(data)
    }
}

subsample_bam_gene_density <- function(gene_list = NULL, data = NULL, color = NULL) {
    observed <- data.frame(value = data$gene_density_signal[match(gene_list, data$gene)])
    observed$group <- "observed"

    gene_count <- length(gene_list)

    simulation <- lapply(1:100, function(x) {
        set.seed(x)
        index <- sample(1:nrow(data), gene_count)
        data1 <- data.frame(value = data$gene_density_signal[index])

        data1$group <- paste0("simulation", x)
        return(data1)
    })
    simulation <- do.call(rbind, simulation)
    data_for_plot <- rbind(observed, simulation)



    p <- ggplot(data_for_plot, aes(x = value, group = group)) +
        geom_density(data = function(x) {
            x[x$group %in% paste("simulation", 1:100, sep = ""), ]
        }, color = "grey", alpha = 0.2, size = 2) +
        geom_density(data = function(x) {
            x[x$group %in% "observed", ]
        }, color = color, size = 2) +
        labs(y = "Density", x = "gene bam average density") +
        # scale_y_continuous(limits = c(0, 2.3)) +
        scale_x_continuous(trans = "log10") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = rgb(0, 72, 144, maxColorValue = 255)),
            strip.text.x = element_text(size = 30, color = "white"),
            plot.title = element_text(
                colour = "black",
                size = 40, vjust = 0.5, hjust = 0.5
            ),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(
                size = 30, color = "black",
            ),
            axis.ticks.length = unit(.25, "cm"),
            axis.text.x = element_text(size = 30, color = "black"),
            axis.title.y = element_text(size = 36, color = "black"),
            axis.title.x = element_text(size = 36, color = "black"),
            legend.title = element_blank(),
            legend.position = "None",
            plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
            panel.spacing = unit(3, "lines")
        )
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/all_genes_bed"



## T stage
chip_depth <- 36971023
input_depth <- 21500156
gene_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/chip.inter"
gene_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/input.inter"
drug_resistance_name <- "/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/drugResistanceMarker"
# 4. variable setting of test module--------------------------------------- TODO:

# 5. process -------------------------------------------------------------- TODO:
all_genes_info <- fread(all_genes_info_filename, stringsAsFactors = F, header = F)
gene_chip_info <- fread(gene_chip_info_filename, stringsAsFactors = F, header = F)
gene_input_info <- fread(gene_input_info_filename, stringsAsFactors = F, header = F)

drug_resistance <- as.character(unlist(read.table(drug_resistance_name, stringsAsFactors = F, header = F)))


colnames(all_genes_info) <- c("chrom", "start_minus_1", "end", "gene_name")
colnames(gene_chip_info) <- colnames(gene_input_info) <- c("chrom", "start_minus_1", "end", "gene_name", "chrom_depth", "start_minus_1_depth", "end_depth", "depth")
## 整合数据
gene_density_signal_data <- get_6mA(
    chip_info = gene_chip_info, input_info = gene_input_info,
    chip_depth = chip_depth, input_depth = input_depth
)
colnames(gene_density_signal_data)[2] <- "gene_density_signal"
colnames(all_genes_info)[4] <- "gene"

data <- merge(gene_density_signal_data, all_genes_info, by = "gene", all = TRUE)

data <- data %>% filter(!is.na(chrom))

drug_resistance_p <- subsample_bam_gene_density(gene_list = drug_resistance, data = data, color = rgb(240, 0, 0, maxColorValue = 255))

ggsave("", width = 14, height = 5)

data_drugResistance <- data %>%
    filter(gene %in% drug_resistance) %>%
    select(gene, gene_density_signal)
write.table(data_drugResistance, "/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/drug_resistance_gene_6mADensity.txt", row.names = F, col.names = T, quote = F, sep = "\t")
