#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(seqinr)))
suppressWarnings(suppressMessages(library(data.table)))
# 2. functions ------------------------------------------------------------ TODO:
read_cds <- function(output_dir = NULL, mrnaid = NULL) {
    fa <- read.fasta(
        file = paste0(output_dir, "/", mrnaid, ".cds.fa"),
        as.string = TRUE, seqtype = "DNA",
        set.attributes = F, forceDNAtolower = FALSE
    )
    return(fa)
}

read_gff <- function(mrnaid = NULL) {
    anno <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/PlasmoDB-36_Pfalciparum3D7.gff"
    data <- fread(cmd = paste0("grep \"", mrnaid, "\" ", anno, " |grep \"CDS\""), stringsAsFactors = F, select = c(1, 3, 4, 5, 7))
    if (unique(data$V7) == "+") {
        data$V3 <- paste(data$V3, 1:NROW(data), sep = "")
    } else {
        data$V3 <- paste(data$V3, NROW(data):1, sep = "")
    }
    colnames(data) <- c("chrom", "feature", "start", "end", "strand")
    return(data)
}

tag_cds <- function(cds = NULL, anno = NULL) {
    # 无论基因strand方向是正还是负，下载到的cds序列已经是翻译方向了。即基本上开头都是ATG
    strand <- anno$strand[1]
    bases <- unlist(strsplit(cds[[1]], ""))
    if (strand == "+") {
        names(bases) <- unlist(Map(`:`, anno$start, anno$end))
    } else if (strand == "-") {
        a <- unlist(Map(`:`, anno$start, anno$end))
        names(bases) <- a[length(a):1]
    }
    base_count <- length(bases)
    first <- bases[seq(1, base_count, by = 3)]
    second <- bases[seq(2, base_count, by = 3)]
    third <- bases[seq(3, base_count, by = 3)]

    result <- list()
    result$first_second <- paste(first, second, sep = "")
    result$third <- third
    return(result)
}

get_ffd_pos <- function(pattern = NULL, cds_tagged = NULL) {
    pattern_index <- lapply(pattern, function(x) {
        which(cds_tagged$first_second == x)
    })
    pattern_index <- do.call(c, pattern_index)
    ffd_pos <- as.integer(names(cds_tagged$third[pattern_index]))
    return(ffd_pos)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
mrnaid <- Args[6]
output_dir <- Args[7]
ffd_bed_name <- paste0(output_dir, "/", mrnaid, ".ffd.bed")
four_fold_degenerate_firsttwo_pattern <- c("GC", "CG", "GG", "CT", "CC", "TC", "AC", "GT") #pattern模式是密码子的，即翻译方向的
# 4. variable setting of test module--------------------------------------- TODO:
# mrnaid<-"PF3D7_0100100.1"
# output_dir<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition"
# 5. process -------------------------------------------------------------- TODO:
cds_seq <- read_cds(output_dir = output_dir, mrnaid = mrnaid)
anno <- read_gff(mrnaid = mrnaid)
cds_tagged <- tag_cds(cds = cds_seq, anno = anno)
ffd_pos <- get_ffd_pos(pattern = four_fold_degenerate_firsttwo_pattern, cds_tagged = cds_tagged)
ffd_bed <- data.frame(chrom = anno$chrom[1], start = ffd_pos - 1, end = ffd_pos)
# write.table(ffd_bed, ffd_bed_name, sep = "\t", col.names = F, row.names = F, quote = F)

write.table(format(ffd_bed, scientific = FALSE), ffd_bed_name, sep = "\t", col.names = F, row.names = F, quote = F)