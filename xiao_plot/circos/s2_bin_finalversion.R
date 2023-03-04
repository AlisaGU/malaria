#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
# library(RCircos)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(circlize)))

# 2. functions ------------------------------------------------------------ TODO:
my_ini <- function(cytoband = system.file(
                       package = "circlize", "extdata",
                       "cytoBand.txt"
                   ), species = NULL, sort.chr = TRUE, chromosome.index = NULL,
                   major.by = NULL, plotType = c("ideogram", "axis", "labels"),
                   track.height = NULL, ideogram.height = convert_height(
                       2,
                       "mm"
                   ), ...) {
    if ("sector.names" %in% names(list(...))) {
        stop_wrap("Found argumnet `sector.names` is used. Please use the argument `chromosome.index` instead.")
    }
    e <- try(cytoband <- read.cytoband(cytoband,
        species = species,
        sort.chr = sort.chr, chromosome.index = chromosome.index
    ),
    silent = TRUE
    )
    if (inherits(e, "try-error") && !is.null(species)) {
        e2 <- try(cytoband <- read.chromInfo(
            species = species,
            sort.chr = sort.chr, chromosome.index = chromosome.index
        ),
        silent = TRUE
        )
        if (inherits(e2, "try-error")) {
            message(e)
            message(e2)
            stop_wrap("Cannot download either cytoband or chromInfo file from UCSC.")
        } else {
            message_wrap("Downloading cytoBand file from UCSC failed. Use chromInfo file instead. Note ideogram track will be removed from the plot.")
            plotType <- setdiff(plotType, "ideogram")
            if (is.null(chromosome.index)) {
                chromInfo <- read.chromInfo(species = species)
                chr_len <- sort(chromInfo$chr.len, decreasing = TRUE)
                i <- which(chr_len[seq_len(length(chr_len) -
                    1)] / chr_len[seq_len(length(chr_len) - 1) +
                    1] > 5)[1]
                if (length(i)) {
                    chromosome <- chromInfo$chromosome[chromInfo$chromosome %in%
                        names(chr_len[chr_len >= chr_len[i]])]
                } else {
                    chromosome <- chromInfo$chromosome
                }
                cytoband <- read.chromInfo(
                    species = species,
                    chromosome.index = chromosome, sort.chr = sort.chr
                )
            }
        }
    } else if (inherits(e, "try-error")) {
        stop(e)
    }
    df <- cytoband$df
    chromosome <- cytoband$chromosome
    if (is.null(chromosome)) {
        if (is.factor(cytoband[, 1])) {
            chromosome <- levels(cytoband$df[, 1])
        } else {
            chromosome <- unique(cytoband$df[, 1])
        }
    }
    if (is.null(chromosome.index)) {
        chromosome.index <- chromosome
    }
    df[[1]] <- factor(as.vector(df[[1]]), levels = chromosome.index)
    sn <- unique(as.vector(df[[1]]))
    sn <- 1:14
    o.cell.padding <- circos.par("cell.padding")
    circos.par(cell.padding = c(
        o.cell.padding[1], 0, o.cell.padding[3],
        0
    ))
    circos.genomicInitialize(df,
        sector.names = sn, major.by = 1000000, tickLabelsStartFromZero = FALSE,
        plotType = plotType, track.height = track.height, ...
    )
    if (any(plotType %in% "ideogram")) {
        circos.genomicIdeogram(df, track.height = ideogram.height)
    }
}


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
pdf(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/circos/circos_binsize_motif100_other10000_trachheight20_6.pdf"),
    width = 5, height = 5
)

circos.par("track.height" = track_height, start.degree = 90)
my_ini(cytoband = genome_structure, plotType = c("labels", "axis"))
# circos.initializeWithIdeogram(chromosome.index = "Pf3D7_01_v3", cytoband = genome_structure, plotType = c("labels", "axis"))
# text(0, 0, "Plasmodium\nfalciparum", cex = 1)

circos.genomicTrack(enzyme_data,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h", col = "#614090")
    }
)

circos.genomicTrack(motif_data,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h", col = "#f0c341")
    }
)

m6A <- list(m6A_data, KD_data)
circos.genomicTrack(m6A,
    panel.fun = function(region, value, ...) {
        i <- getI(...)
        # circos.genomicLines(region, value, col = i, area = TRUE)
        circos.genomicLines(region, value, col = ifelse(i == 1, "black", "#c8352d"), area = FALSE, type = "l")
    }
)

dev.off()
circos.clear()
