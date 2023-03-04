#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
# library(RCircos)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(circlize)))

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
binsize <- 10000
enzyme_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/histone_modification/protein/enzyme_log2ratio_bin", binsize, ".bdg"), header = F, stringsAsFactors = F)
m6A_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/WT_log2ratio_bin", binsize, ".bdg"), header = F, stringsAsFactors = F)
KD_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/KD_log2ratio_bin", binsize, ".bdg"), header = F, stringsAsFactors = F)
motif_data <- fread(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/GAAGAA_genome_position/GAAGAA_bin", binsize, "_count.bdg"), header = F, stringsAsFactors = F)
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
pdf(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos/circos_binsize", binsize, "bp_trackheight", track_height * 100, ".pdf"),
    width = 10, height = 10
)

circos.par("track.height" = track_height, start.degree = 90)
circos.initializeWithIdeogram(cytoband = genome_structure, plotType = c("labels", "axis"))
# circos.initializeWithIdeogram(chromosome.index = "Pf3D7_01_v3", cytoband = genome_structure, plotType = c("labels", "axis"))
text(0, 0, "Plasmodium\nfalciparum", cex = 1)

circos.genomicTrack(motif_data,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h", col = "#e600ff")
    }
)

circos.genomicTrack(enzyme_data,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h", col = "blue")
    }
)

m6A <- list(m6A_data, KD_data)
circos.genomicTrack(m6A,
    panel.fun = function(region, value, ...) {
        i <- getI(...)
        circos.genomicLines(region, value, col = i, area = TRUE)
    }
)



dev.off()
circos.clear()
