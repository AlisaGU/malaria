#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(seqinr)
library(stringr)
library(parallel)
# 2. functions ------------------------------------------------------------ TODO:
collapse_region <- function(run = NULL) {
    rundiff <- c(1, diff(run))
    difflist <- split(run, cumsum(rundiff != 1))
    result <- t(sapply(difflist, function(x) {
        c(x[1], x[length(x)])
    }))
    result[, 1] <- result[, 1] - 1
    colnames(result) <- c("Start_minus_1", "End")
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
seq_data <- read.fasta(
    file = Args[6],
    seqtype = "DNA", as.string = T
)
chrom <- Args[7]
output <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/", chrom, ".outgroup_region.bed")

# 4. variable setting of test module--------------------------------------- TODO:
# seq_data<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/Pf3D7_01_v3.include_P_billcollinsi.fasta"
# seq_data <- read.fasta(
#     file = seq_data,
#     seqtype = "DNA", as.string = T
# )
# 5. process -------------------------------------------------------------- TODO:
P_falciparum <- s2c(seq_data[[1]])
# outgroup_A <- s2c(seq_data[[3]])
# interested_B <- s2c(seq_data[[2]])
outgroup_A <- s2c(seq_data[[2]])
interested_B <- s2c(seq_data[[3]])

pos <- which(outgroup_A != "n" & outgroup_A != "-")


gap_in_ref <- which(P_falciparum == "-")
pos <- mclapply(1:length(pos), function(x) {
    pos[x] - length(which(gap_in_ref <= pos[x]))
}, mc.cores = 4)
pos <- do.call(c, pos)
consistent_region_bed <- collapse_region(run = unique(pos))
result <- data.frame(chrom = chrom, consistent_region_bed)
write.table(result, output, quote = F, row.names = F, col.names = F, sep = "\t")
