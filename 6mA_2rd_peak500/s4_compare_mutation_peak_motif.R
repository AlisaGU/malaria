#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(seqinr)
library(stringr)
# 2. functions ------------------------------------------------------------ TODO:
find_seq_with_motif <- function(motif = NULL, seqs = NULL) {
    motif_index <- sapply(seqs, function(x) {
        grepl(motif, x)
    })
    seq_with_motif <- seqs[motif_index]
    return(seq_with_motif)
}

find_motif_and_surr_index_in_genome <- function(motif = NULL, seqs = NULL) {
    seq_count <- length(seqs)
    peak_name <- as.character(sapply(names(seqs), function(x) {
        unlist(strsplit(x, "::"))[2]
    }))
    result <- lapply(1:seq_count, function(index) {
        peak_region <- unlist(strsplit(peak_name[index], ":|-", perl = T))
        chrom <- peak_region[1]
        start <- as.numeric(peak_region[2])
        end <- as.numeric(peak_region[3])

        relative_index <- str_locate_all(seqs[[index]], motif)[[1]]
        absolute_index <- get_motif_absolute_index(peak_start = start, peak_end = end, motif_relative_index = relative_index)
        result <- list(motif = data.frame(peak = peak_name[index], absolute_index$motif), motif_surr = data.frame(peak = peak_name[index], absolute_index$motif_surr))
    })
    names(result) <- peak_name

    motif_index_genome <- lapply(result, function(x) {
        x$motif
    })
    motif_index_genome <- do.call(rbind, motif_index_genome)
    rownames(motif_index_genome) <- NULL

    motif_surr_index_genome <- lapply(result, function(x) {
        x$motif_surr
    })
    motif_surr_index_genome <- do.call(rbind, motif_surr_index_genome)
    rownames(motif_surr_index_genome) <- NULL
    colnames(motif_surr_index_genome) <- c("peak", "start", "end")
    result <- list(motif = motif_index_genome, motif_surr = motif_surr_index_genome)
    return(result)
}

get_motif_absolute_index <- function(peak_start = NULL, peak_end = NULL, motif_relative_index = NULL) {
    peak_region_length <- peak_end - peak_start + 1
    motif_relative_index_vector <- unlist(Map(`:`, motif_relative_index[, 1], motif_relative_index[, 2]))

    motif_surr_relative_index <- (1:peak_region_length)[-motif_relative_index_vector]
    rundiff <- c(1, diff(motif_surr_relative_index))
    difflist <- split(motif_surr_relative_index, cumsum(rundiff != 1))
    motif_surr_relative_index <- t(sapply(difflist, function(x) {
        c(x[1], x[length(x)])
    }))

    motif_absolute_index <- peak_start + motif_relative_index
    motif_surr_absolute_index <- peak_start + motif_surr_relative_index

    return(list(motif = motif_absolute_index, motif_surr = motif_surr_absolute_index))
}
# 3. input ---------------------------------------------------------------- TODO:
motif_bed_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif/GAAGAA"
top500Peak_seq <- read.fasta("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/top500/3D7-T3_summits-500-Top500.fa",
    seqtype = "DNA", as.string = T, set.attributes = F, forceDNAtolower = F
)

motif <- "GAAGAA"
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peakseq_with_motif <- find_seq_with_motif(motif = motif, seqs = top500Peak_seq)
motif_and_surr_index_genome <- find_motif_and_surr_index_in_genome(motif = motif, seqs = peakseq_with_motif)
motif_index_genome <- motif_and_surr_index_genome$motif
motif_surr_index_genome <- motif_and_surr_index_genome$motif_surr

motif_index_genome