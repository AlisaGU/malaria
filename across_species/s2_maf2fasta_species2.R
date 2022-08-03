#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(rphast)
library(seqinr)
library(stringr)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
Type <- Args[6]
setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/", Type))

dataname <- Args[7]

# 4. variable setting of test module--------------------------------------- TODO:
# dataname <- "Pf3D7_01_v3.P_falciparum.P_reichenowi.P_praefalciparum.filter.maf"
# 5. process -------------------------------------------------------------- TODO:
chrom <- unlist(strsplit(dataname, ".", fixed = T))[1]
msa_data <- read.msa(dataname, format = "MAF")

seq1_N_index <- str_locate_all(msa_data$seq[1], "N")[[1]][, 1] # block连接处
seq2_N_index <- str_locate_all(msa_data$seq[2], fixed("*"))[[1]][, 1] # block连接处+gap处

seq1_gap_index <- setdiff(seq1_N_index, seq1_N_index)
seq2_gap_index <- setdiff(seq2_N_index, seq1_N_index)

gap_index <- list(seq1_gap_index, seq2_gap_index)

for (i in 1:length(msa_data$names)) {
    seqname <- msa_data$names[i]
    seq <- msa_data$seq[i]
    if (i != 1) {
        seq <- gsub("*", "N", seq, fixed = T)
        a <- s2c(seq)
        a[gap_index[[i]]] <- "-"
        seq <- paste(a, collapse = "")
    }
    write.fasta(
        sequences = paste0(paste(rep("N", msa_data$offset), collapse = ""), seq),
        names = seqname, file.out = paste0(chrom, ".", seqname, ".fasta")
    )
}
