#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
# library(RCircos)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(circlize)))

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
enzyme_data <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/enzyme/wholeGenome/merged.bed", header = F, stringsAsFactors = F)
m6A_data <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/histone_modification/6mA/merged.bed", header = F, stringsAsFactors = F)
KD_data <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/merged.bed", header = F, stringsAsFactors = F)
motif_data <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/GAAGAA_genome_position/GAAGAA.bed", header = F, stringsAsFactors = F)
genome_structure <- fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length", header = F, stringsAsFactors = F)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(enzyme_data) <- colnames(m6A_data) <- colnames(KD_data) <- c("Chromosome", "ChromStart", "ChromEnd", "density")
colnames(motif_data)[1] <- c("Chromosome")
enzyme_data <- enzyme_data[order(enzyme_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
m6A_data <- m6A_data[order(m6A_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
KD_data <- KD_data[order(KD_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
motif_data <- motif_data[order(motif_data[, 1], decreasing = FALSE), ] %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")

genome_structure <- data.frame(
    Chromosome = genome_structure$V1, ChromStart = 0,
    ChromEnd = genome_structure$V2, Band = "1", Stain = "gneg"
) %>% filter(Chromosome != "Pf3D7_API_v3" & Chromosome != "Pf_M76611")
motif_data$status <- ifelse(motif_data$V6 == "+", 1, -1)
motif_data <- motif_data %>% select(Chromosome, V2, V3, status)
colnames(motif_data) <- c("Chromosome", "ChromStart", "ChromEnd", "status")
## 画基因组
# pdf("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos/circos.pdf", width = 10, height = 10)

# RCircos.Set.Core.Components(
#     cyto.info = genome_structure,
#     chr.exclude = NULL,
#     tracks.inside = 10, tracks.outside = 0
# )
# rcircos.params <- RCircos.Get.Plot.Parameters()
# rcircos.params$base.per.unit <- 60000
# RCircos.Reset.Plot.Parameters(rcircos.params)
# RCircos.Set.Plot.Area()
# RCircos.Chromosome.Ideogram.Plot()
# ##
# RCircos.Scatter.Plot(motif_data, data.col = 4, track.num = 1, side = "in", by.fold = 0.5)
# # RCircos.Line.Plot(enzyme_data, data.col = 4, track.num = 6, side = "in")
# dev.off()




### circlize
pdf("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos/circos_circlize.pdf", width = 10, height = 10)

# colnames(genome_structure)<-c("chr","start","end","band","color")
circos.par("track.height" = 0.15)
circos.initializeWithIdeogram(cytoband = genome_structure, plotType = c("labels", "axis"))
text(0, 0, "Plasmodium\nfalciparum", cex = 1)
# circos.track(ylim = c(0, 1))
# circos.genomicIdeogram(cytoband = genome_structure, track.height = 0.2)
# circos.genomicHeatmap(motif_data %>% filter(Chromosome == "Pf3D7_01_v3"), side = "inside", colorRamp2(c(-1,  1), c("green", "red")),border = "white")
{
    # circos.genomicDensity(motif_data %>% filter(Chromosome == "Pf3D7_01_v3"), col = "red", border = "black")
    circos.genomicDensity(motif_data, col = "blue", border = "black", window.size = 200000)
}

# circos.genomicDensity(motif_data %>% filter(Chromosome == "Pf3D7_01_v3"), count_by = "number", col = "red", border = "black")
# circos.text(1, 0.1, "right", sector.index = "Pf3D7_01_v3", track.index = 3)
{
    enzyme_data_subset <- enzyme_data
    circos.genomicTrack(enzyme_data_subset,
        panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, type = "h", col = "green")
        }
    )
}

{
    m6A_data_subset <- m6A_data
    circos.genomicTrack(m6A_data_subset,
        panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, type = "h", col = "red")
        }
    )
}

{
    KD_data_subset <- KD_data
    circos.genomicTrack(KD_data_subset,
        panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, type = "h", col = "#ff8f6b")
        }
    )
}

dev.off()
circos.clear()
