#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:
get_based_vcf <- function(vcf = NULL) {
    based_vcf <- t(apply(vcf, 1, function(x) {
        allele <- c(x[3], unlist(strsplit(x[4], ",")))
        names(allele) <- 0:(length(allele) - 1)
        x[5:7] <- allele[x[5:7]]
        # x[5:7] <- gsub("A|T|C|G", "N", x[5:7])
        return(x)
    }))
    based_vcf[, 2] <- gsub(" ", "", based_vcf[, 2])
    return(based_vcf)
}

# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
vcf_data <- Args[6]
output <- paste0(dirname(vcf_data), "/", gsub(".vcf", "", basename(vcf_data)), "_polysite.bed")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
vcf <- fread(cmd = paste0("cat ", vcf_data, "|grep -v \"^#\"|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s\",b[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, c(8) := NULL]

colnames(vcf) <- c("chrom", "pos", "ref", "alt", "P_falciparum", "outgroupA", "B_species")
based_vcf <- get_based_vcf(vcf = vcf)
insertion_index <- based_vcf[based_vcf[, 5] == "*" & based_vcf[, 6] == "*", 1:2]
deletion_index <- based_vcf[based_vcf[, 5] != "*" & based_vcf[, 6] != "*", 1:2]
indel_data <- unique(rbind(insertion_index, deletion_index))
result <- data.frame(chrom = as.character(indel_data[1, 1]), pos_minus_1 = as.numeric(indel_data[, 2]) - 1, pos = as.numeric(indel_data[, 2]))
write.table(result, quote = F, file = gsub("intermediate1", "indel_possible.bed", vcf_data), sep = "\t", row.names = F, col.names = F)
