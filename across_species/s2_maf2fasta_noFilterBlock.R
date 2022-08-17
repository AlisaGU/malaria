#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(rphast)
library(seqinr)
library(stringr)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf")
Args <- commandArgs()
chrom <- Args[6]

# 4. variable setting of test module--------------------------------------- TODO:
# chrom <- "Pf3D7_11_v3"

# 5. process -------------------------------------------------------------- TODO:
datanames_count <- length(grep(paste0(chrom, ".file"), dir()))
datanames <- paste(paste0(chrom, ".file"), 1:datanames_count, sep = "")
offset <- read.msa(datanames[[1]], format = "MAF")$offset
for (dataname in datanames) {
    msa_data <- read.msa(dataname, format = "MAF")
}

dataseqs <- lapply(datanames, function(dataname) {
    result <- list()
    msa_data <- read.msa(dataname, format = "MAF")
    only_one_seq <- length(msa_data$names) == 1
    if (only_one_seq) {
        base_count <- nchar(msa_data$seq[[1]])
        result$P_falciparum <- paste(rep("-", base_count), collapse = "")

        if (msa_data$names == "P_billcollinsi") {
            result$P_billcollinsi <- msa_data$seq[[1]]
            result$P_reichenowi <- paste(rep("-", base_count), collapse = "")
        }
        if (msa_data$names == "P_reichenowi") {
            result$P_billcollinsi <- paste(rep("-", base_count), collapse = "")
            result$P_reichenowi <- msa_data$seq[[1]]
        }
    } else {
        result$P_falciparum <- msa_data$seq[[which(msa_data$names == "P_falciparum")]]
        result$P_billcollinsi <- msa_data$seq[[which(msa_data$names == "P_billcollinsi")]]
        result$P_reichenowi <- msa_data$seq[[which(msa_data$names == "P_reichenowi")]]
    }
    return(result)
})

concat_seq <- sapply(c("P_falciparum", "P_billcollinsi", "P_reichenowi"), function(species) {
    species_seq <- sapply(dataseqs, function(x) {
        x[[species]]
    })
    result <- paste(species_seq, collapse = "")
    return(result)
})


seq1_N_index <- str_locate_all(concat_seq[["P_falciparum"]], "N")[[1]][, 1] # block连接处
seq2_N_index <- str_locate_all(concat_seq[["P_billcollinsi"]], fixed("*"))[[1]][, 1] # block连接处+gap处
seq3_N_index <- str_locate_all(concat_seq[["P_reichenowi"]], fixed("*"))[[1]][, 1] # block连接处+gap处

seq1_gap_index <- setdiff(seq1_N_index, seq1_N_index)
seq2_gap_index <- setdiff(seq2_N_index, seq1_N_index)
seq3_gap_index <- setdiff(seq3_N_index, seq1_N_index)

gap_index <- list(seq1_gap_index, seq2_gap_index, seq3_gap_index)

for (i in 1:length(concat_seq)) {
    seqname <- names(concat_seq)[i]
    seq <- concat_seq[[i]]
    if (i != 1) {
        seq <- gsub("*", "N", seq, fixed = T)
        a <- s2c(seq)
        a[gap_index[[i]]] <- "-"
        seq <- paste(a, collapse = "")
    }
    write.fasta(
        sequences = paste0(paste(rep("N", offset), collapse = ""), seq),
        names = seqname, file.out = paste0(chrom, ".", seqname, ".noFilterBlock.fasta")
    )
}
