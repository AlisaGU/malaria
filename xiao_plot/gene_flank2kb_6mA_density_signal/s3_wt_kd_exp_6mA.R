#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

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

# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
WT_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/chip.gene_flank2kb.inter.Tstage"
WT_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/input.gene_flank2kb.inter.Tstage"
WT_chip_depth <- 36971023
WT_input_depth <- 21500156

KD_chip_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/chip.gene_flank2kb.inter"
KD_input_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/input.gene_flank2kb.inter"
KD_chip_depth <- 39189761
KD_input_depth <- 27773718



all_genes_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_flank2kb_6mA_density_signal/gene_flank2kb.bed"
WT_RNA_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/3D7_RNA_Tstage_new.txt" # tpm

# 5. process -------------------------------------------------------------- TODO:
gene_WT_6mA_density <- get_6mA(chip_info_name = WT_chip_info_filename, input_info_name = WT_input_info_filename, chip_depth = WT_chip_depth, input_depth = WT_input_depth)
gene_KD_6mA_density <- get_6mA(chip_info_name = KD_chip_info_filename, input_info_name = KD_input_info_filename, chip_depth = KD_chip_depth, input_depth = KD_input_depth)
colnames(gene_WT_6mA_density)[2] <- "WT"
colnames(gene_KD_6mA_density)[2] <- "KD"


RNA <- fread(RNA_filename, stringsAsFactors = F, header = T) %>% na.omit()
