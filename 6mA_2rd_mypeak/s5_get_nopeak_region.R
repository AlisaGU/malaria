#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:
get_other_region_bed_for_specific_chrom <- function(chrom = NULL, chrom.len = NULL, chrom_region_bed = NULL) {
    chrom_vector <- rep(0, chrom.len)

    chrom_region_position <- chrom_region_bed
    chrom_region_position$start <- chrom_region_position$start + 1
    peak_region <- Map(`:`, chrom_region_position$start, chrom_region_position$end)
    peak_region <- unlist(peak_region)
    chrom_vector[peak_region] <- 1

    nopeak_region_index <- which(chrom_vector == 0)
    rundiff <- c(1, diff(nopeak_region_index))
    difflist <- split(nopeak_region_index, cumsum(rundiff != 1))
    nopeak_region_position <- t(sapply(difflist, range))
    nopeak_region_bed <- cbind(nopeak_region_position[, 1] - 1, nopeak_region_position[, 2])

    result <- data.frame(chrom = chrom, start = nopeak_region_bed[, 1], end = nopeak_region_bed[, 2])
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
chrom_length_filename <- Args[6]
region_bed_filename <- Args[7]
other_region_bed_filename <- Args[8]
# 4. variable setting of test module--------------------------------------- TODO:
# chrom_length_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length"
# region_bed_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"
# other_region_bed_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/nopeak.bed"
# 5. process -------------------------------------------------------------- TODO:
chrom_length <- fread(chrom_length_filename, sep = "\t", stringsAsFactors = F)
colnames(chrom_length) <- c("chrom", "length")
region_bed <- fread(region_bed_filename, sep = "\t", stringsAsFactors = F)
colnames(region_bed) <- c("chrom", "start", "end")

chroms <- unique(region_bed$chrom)
result <- lapply(chroms, function(chrom) {
    chrom.len <- chrom_length$length[which(chrom_length$chrom == chrom)]
    index_in_region <- which(region_bed$chrom == chrom)
    chrom_region_bed <- region_bed[index_in_region, ]
    result <- get_other_region_bed_for_specific_chrom(
        chrom = chrom,
        chrom.len = chrom.len,
        chrom_region_bed = chrom_region_bed
    )
})

nopeak_region <- do.call(rbind, result)
write.table(nopeak_region, other_region_bed_filename, quote = F, sep = "\t", col.names = F, row.names = F)