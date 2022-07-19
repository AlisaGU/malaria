#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:
library(Biostrings)
library(seqinr)

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_flank_pattern")
motifs <- c("GAAGAA", "TTCTTC", "GAAGAT", "ATCTTC", "GATGAA", "TTCATC", "GATGAT", "ATCATC")
positive_motifs <- c("GAAGAA", "GAAGAT", "GATGAA", "GATGAT")
negative_motifs <- c("TTCTTC", "ATCTTC", "TTCATC", "ATCATC")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
flank_3_info <- read.table("flank_3_info.not_remove_inconsistent", header = F, as.is = T)

context <- flank_3_info$V8

seqs <- t(sapply(context, function(x) {
    c(substr(x, 1, 12), substr(x, 13, 18), substr(x, 19, 20), substr(x, 21, 21), substr(x, 22, 23), substr(x, 24, 29), substr(x, 30, 41))
}))

bed_3_6 <- t(apply(cbind(flank_3_info, seqs), 1, function(x) {
    x <- unlist(x)
    strand <- x[5]
    posi.motif_index.x <- which(sapply(x, function(y) {
        y %in% positive_motifs
    }))
    nega.motif_index.x <- which(sapply(x, function(y) {
        y %in% negative_motifs
    }))
    index_3 <- 0
    index_6 <- 0

    if (strand == "+" & length(nega.motif_index.x) == 1 && nega.motif_index.x == 15) {
        # 第1行
        index_3 <- as.numeric(x[3]) - 5
        index_6 <- as.numeric(x[3]) - 8
    } else if (strand == "+" & length(nega.motif_index.x) == 1 && nega.motif_index.x == 19) {
        # 第124行
        index_3 <- as.numeric(x[3]) + 6
        index_6 <- as.numeric(x[3]) + 3
    } else if (strand == "+" & length(posi.motif_index.x) == 1 && posi.motif_index.x == 19) {
        # 第4行
        index_3 <- as.numeric(x[3]) + 5
        index_6 <- as.numeric(x[3]) + 8
    } else if (strand == "+" & length(posi.motif_index.x) == 1 && posi.motif_index.x == 15) {
        # 第8行
        index_3 <- as.numeric(x[3]) - 6
        index_6 <- as.numeric(x[3]) - 3
    } else if (strand == "-" & length(posi.motif_index.x) == 1 && posi.motif_index.x == 15) {
        # 第7行
        index_3 <- as.numeric(x[3]) + 6
        index_6 <- as.numeric(x[3]) + 3
    } else if (strand == "-" & length(posi.motif_index.x) == 1 && posi.motif_index.x == 19) {
        # 第3行
        index_3 <- as.numeric(x[3]) - 5
        index_6 <- as.numeric(x[3]) - 8
    } else if (strand == "-" & length(nega.motif_index.x) == 1 && nega.motif_index.x == 15) {
        # 第95行
        index_3 <- as.numeric(x[3]) + 5
        index_6 <- as.numeric(x[3]) + 8
    } else if (strand == "-" & length(nega.motif_index.x) == 1 && nega.motif_index.x == 19) {
        # 第151行
        index_3 <- as.numeric(x[3]) - 6
        index_6 <- as.numeric(x[3]) - 3
    }
    return(c(index_3, index_6))
}))
rownames(bed_3_6) <- NULL

motif_3_6_9 <- cbind(flank_3_info, seqs, bed_3_6)
write.table(motif_3_6_9, "motif_3_6_9.txt", row.names = F, col.names = F, sep = "\t", quote = F)

motif_9bp <- apply(seqs, 1, function(x) {
    x <- unlist(x)
    posi.motif_index.x <- which(sapply(x, function(y) {
        y %in% positive_motifs
    }))
    nega.motif_index.x <- which(sapply(x, function(y) {
        y %in% negative_motifs
    }))
    if (length(posi.motif_index.x) == 1 && posi.motif_index.x == 2) {
        return(paste(x[2], x[3], x[4], sep = ""))
    } else if (length(posi.motif_index.x) == 1 && posi.motif_index.x == 6) {
        return(paste(x[4], x[5], x[6], sep = ""))
    } else if (length(nega.motif_index.x) == 1 && nega.motif_index.x == 2) {
        return(as.character(reverseComplement(DNAString(paste(x[2], x[3], x[4], sep = "")))))
    } else if (length(nega.motif_index.x) == 1 && nega.motif_index.x == 6) {
        return(as.character(reverseComplement(DNAString(paste(x[4], x[5], x[6], sep = "")))))
    }
})
flank_3_info <- data.frame(flank_3_info, seqs)
rownames(flank_3_info) <- NULL
flank_3_info <- flank_3_info[, -8]
colnames(flank_3_info) <- c("Seqid", "Start", "End", "Score", "Strand", "Phase", "Coverage", "IPDRatio", "frac", "fracLow", "fracUp", "identificationQv", "left_seq1", "left_seq2", "left_seq3", "modified_base", "right_seq1", "right_seq2", "right_seq3")
write.table(flank_3_info, "flank_3_info_seq.txt", sep = "\t", quote = F, row.names = F)

flank3_9bp <- data.frame(flank_3_info, motif_9bp)
rownames(flank3_9bp) <- NULL
colnames(flank3_9bp) <- c("Seqid", "Start", "End", "Score", "Strand", "Phase", "Coverage", "seqs", "IPDRatio", "frac", "fracLow", "fracUp", "identificationQv", "9bp")

write.table(flank3_9bp, "flank3_9bp.txt", sep = "\t", quote = F, row.names = F)



motif_27bp <- t(apply(seqs, 1, function(x) {
    x <- unlist(x)
    posi.motif_index.x <- which(sapply(x, function(y) {
        y %in% positive_motifs
    }))
    nega.motif_index.x <- which(sapply(x, function(y) {
        y %in% negative_motifs
    }))
    if (length(posi.motif_index.x) == 1 && posi.motif_index.x == 2) {
        return(c(
            substr(x[1], nchar(x[1]) - 8, nchar(x[1])),
            paste(x[2], x[3], x[4], sep = ""),
            paste0(x[5], x[6], substr(x[7], 1, 1))
        ))
    } else if (length(posi.motif_index.x) == 1 && posi.motif_index.x == 6) {
        return(c(
            paste0(substr(x[1], nchar(x[1]), nchar(x[1])), x[2], x[3]),
            paste(x[4], x[5], x[6], sep = ""),
            substr(x[7], 1, 9)
        ))
    } else if (length(nega.motif_index.x) == 1 && nega.motif_index.x == 2) {
        return(
            c(
                as.character(reverseComplement(DNAString(paste0(x[5], x[6], substr(x[7], 1, 1))))),
                as.character(reverseComplement(DNAString(paste(x[2], x[3], x[4], sep = "")))),
                as.character(reverseComplement(DNAString(substr(x[1], nchar(x[1]) - 8, nchar(x[1])))))
            )
        )
    } else if (length(nega.motif_index.x) == 1 && nega.motif_index.x == 6) {
        return(c(
            as.character(reverseComplement(DNAString(substr(x[7], 1, 9)))),
            as.character(reverseComplement(DNAString(paste(x[4], x[5], x[6], sep = "")))),
            as.character(reverseComplement(DNAString(paste0(substr(x[1], nchar(x[1]), nchar(x[1])), x[2], x[3]))))
        ))
    }
}))
flank3_27bp <- data.frame(flank_3_info, motif_27bp)
rownames(flank3_27bp) <- NULL
colnames(flank3_27bp) <- c("Seqid", "Start", "End", "Score", "Strand", "Phase", "Coverage", "seqs", "IPDRatio", "frac", "fracLow", "fracUp", "identificationQv", "upstream_9bp", "motif_9bp", "downstream_9bp")

write.table(flank3_27bp, "flank3_27bp.txt", sep = "\t", quote = F, row.names = F)
