#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(seqinr)
library(stringr)
library(data.table)
# 2. functions ------------------------------------------------------------ TODO:
merge_variants <- function(data = NULL) {
    a <- table(data$pos_minus_1)
    unique_index <- a[a == 1]
    unique_pos <- names(unique_index)
    unique_data <- data[match(unique_pos, data$pos_minus_1), ]
    duplicated_data <- data[-match(unique_pos, data$pos_minus_1), ]
    data_split <- split(duplicated_data, duplicated_data$pos_minus_1, drop = T)
    deduplicated <- as.data.frame(t(sapply(data_split, function(x) {
        result <- c(chrom = unlist(x[1, 1]), pos_minus_1 = unlist(x[1, 2]), pos = unlist(x[1, 3]), colSums(x[, -c(1:3)]))
    })))
    deduplicated[, -1] <- sapply(deduplicated[, -1], as.integer)
    result <- rbind(unique_data, deduplicated)
    result <- result[order(result$pos_minus_1), ]
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf")

Args <- commandArgs()
bed_name <- Args[6]
raw_bed <- fread(bed_name, stringsAsFactors = F, header = F)
align_fasta <- read.fasta(Args[7], seqtype = "DNA", as.string = T)
# 4. variable setting of test module--------------------------------------- TODO:
# raw_bed <- fread("Pf3D7_01_v3.raw.bed", stringsAsFactors = F, header = T)
# align_fasta <- read.fasta("Pf3D7_01_v3.P_falciparum.P_billcollinsi.P_reichenowi.noFilterBlock.fasta", seqtype = "DNA", as.string = T)

# 5. process -------------------------------------------------------------- TODO:
colnames(raw_bed) <- c("chrom", "pos_minus_1", "pos", "ts", "tv", "snp", "ins", "del", "indel")
raw_bed <- merge_variants(data = raw_bed)
ref_fasta <- align_fasta["P_falciparum"]
gap_in_ref <- str_locate_all(ref_fasta, "-")[[1]][, 1]

raw_bed$pos <- sapply(raw_bed$pos, function(x) {
    x - length(which(gap_in_ref <= x))
})

raw_bed$pos_minus_1 <- raw_bed$pos - 1
raw_bed <- unique(raw_bed)

write.table(raw_bed, file = gsub(".nonscaled_genome_pos.indel.bed", ".scaled_genome_pos.bed", bed_name), sep = "\t", quote = F, row.names = F, col.names = F)
