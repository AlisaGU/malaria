#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
# library(RCircos)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(RCircos)))

# 2. functions ------------------------------------------------------------ TODO:



# 3. input ---------------------------------------------------------------- TODO:
enzyme_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/histone_modification/protein/enzyme_log2ratio_bin", 10000, ".bdg"), header = F, stringsAsFactors = F)
m6A_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/WT_log2ratio_bin", 10000, ".bdg"), header = F, stringsAsFactors = F)
KD_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/KD_log2ratio_bin", 10000, ".bdg"), header = F, stringsAsFactors = F)
motif_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/GAAGAA_genome_position/GAAGAA_bin", 100, "_count.bdg"), header = F, stringsAsFactors = F)
genome_structure <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length", header = F, stringsAsFactors = F)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(motif_data) <- colnames(enzyme_data) <- colnames(m6A_data) <- colnames(KD_data) <- c("Chromosome", "ChromStart", "ChromEnd", "density")

enzyme_data <- enzyme_data[order(enzyme_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
m6A_data <- m6A_data[order(m6A_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
KD_data <- KD_data[order(KD_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
motif_data <- motif_data[order(motif_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")

genome_structure <- data.frame(
    Chromosome = genome_structure$V1, ChromStart = 0,
    ChromEnd = genome_structure$V2, Band = "1", Stain = "gneg"
) %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")


### circlize
track_height <- 0.2
pdf(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos/circos_binsize_motif100_other10000_trachheight20_7.pdf"),
    width = 5, height = 5
)
RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram, # 这是上面load的基因组文件
    chr.exclude <- NULL, # 这个参数是要排除的染色体，这里我选无，也就是画出所有的染色体，你也可以排除x和Y染色体
    tracks.inside = 4, # 这一个参数是指定在染色体圆圈内部一共要画几个圆圈
    tracks.outside = 0
)
dev.off()
