#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:
get_around_genome_region <- function(chrom = NULL, start = NULL, end = NULL) {
    # start and end are both 1-base
    chrom_length.i <- chrom_length$length[which(chrom_length$chrom == chrom)]
    left_index <- start - 500
    right_index <- end + 500
    # if(left_index<=0 || right_index>chrom_length.i){
    #     return(c(chrom, start, end, chrom_length.i))
    # }
    around_index <- c()

    if (start == 1) {
        around_index <- c(chrom, end + 1, end + 1000)
    }

    if (left_index >= 1 && right_index <= chrom_length.i) {
        around_index <- rbind(c(chrom, left_index, start - 1), c(chrom, end + 1, right_index))
    }
    if (left_index <= 0 & start != 1) {
        left_borrow_len <- length(left_index:0)
        around_index <- rbind(c(chrom, 1, start - 1), c(chrom, end + 1, right_index + left_borrow_len))
    }
    if (right_index > chrom_length.i) {
        right_borrow_len <- length((chrom_length.i + 1):right_index)
        around_index <- rbind(c(chrom, left_index - right_borrow_len, start - 1), c(chrom, end + 1, chrom_length.i))
    }

    return(around_index)
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/top500")
top500_header <- fread(cmd = "grep \"^>\" \"3D7-T3_summits-500-Top500.fa\"|awk -F\"::\" '{print $2}'", stringsAsFactors = T, header = F)
chrom_length <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length", stringsAsFactors = F)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
top500_header <- as.character(unlist(top500_header))
colnames(chrom_length) <- c("chrom", "length")
genome_position <- t(sapply(top500_header, function(x) {
    chrom <- unlist(strsplit(x, ":"))[1]
    start_end <- unlist(strsplit(x, ":"))[2]
    start <- as.numeric(unlist(strsplit(start_end, "-"))[1]) - 1
    end <- as.numeric(unlist(strsplit(start_end, "-"))[2])
    return(c(chrom, start, end))
}))

around_index <- c()
for (i in 1:nrow(genome_position)) {
    chrom <- genome_position[i, 1]
    start <- as.numeric(genome_position[i, 2]) + 1
    end <- as.numeric(genome_position[i, 3])

    around_index <- rbind(around_index, get_around_genome_region(chrom = chrom, start = start, end = end))
}
around_index[, 2] <- as.numeric(around_index[, 2]) - 1

top500_and_surr_position <- rbind(genome_position, around_index)
write.table(top500_and_surr_position, "top500_and_surr.bed", quote = F, row.names = F, col.names = F, sep = "\t")