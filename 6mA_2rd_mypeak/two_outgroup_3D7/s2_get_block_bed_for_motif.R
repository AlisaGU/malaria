#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
# 2. functions ------------------------------------------------------------
remove_file <- function(file = NULL) {
    if (file.exists(file)) {
        invisible(file.remove(file))
    }
}
extract_block_loc <- function(motif_info_in_single_seq = NULL) {
    motif_count <- nrow(motif_info_in_single_seq)

    result <- ""
    if (motif_count == 1) {
        result <- block_bed_for_single_motif(motif_info_in_single_seq = motif_info_in_single_seq)
    } else {
        result <- block_bed_for_multiple_motifs(motif_info_in_single_seq = motif_info_in_single_seq)
    }
    return(result)
}

block_bed_for_single_motif <- function(motif_info_in_single_seq = NULL) {
    seq_info <- unlist(strsplit(motif_info_in_single_seq$seq_id, "_|-|:", perl = T))
    seq_start_in_genome <- as.numeric(seq_info[4])
    seq_end_in_genome <- as.numeric(seq_info[5])
    start_in_seq <- as.numeric(motif_info_in_single_seq$start_in_seq)
    end_in_seq <- as.numeric(motif_info_in_single_seq$end_in_seq)
    seq_length <- seq_end_in_genome - seq_start_in_genome + 1

    region_loc_for_one_motif <- get_region_loc_for_one_motif(
        start_in_seq = start_in_seq,
        end_in_seq = end_in_seq,
        seq_length = seq_length
    )
    motif_region_loc <- list(
        motif_loc = c(start_in_seq, end_in_seq),
        region_1_loc = vector2range(run = region_loc_for_one_motif$region_1_loc),
        region_2_loc = vector2range(run = region_loc_for_one_motif$region_2_loc),
        region_3_loc = vector2range(run = region_loc_for_one_motif$region_3_loc),
        region_4_loc = vector2range(run = region_loc_for_one_motif$region_4_loc)
    )
    result <- lapply(motif_region_loc, function(x) {
        if (is.vector(x)) {
            x[1] <- x[1] - 1
        } else {
            x[, 1] <- x[, 1] - 1
        }
        return(x)
    })
    result1 <- lapply(result, function(x) {
        if (!is.null(x)) {
            x + seq_start_in_genome - 1
        }
    })
    return(result1)
}

get_region_loc_for_one_motif <- function(start_in_seq = NULL, end_in_seq = NULL, seq_length = NULL) {
    region_1_loc <- get_region_loc(start_in_seq = start_in_seq, end_in_seq = end_in_seq, region_number = 1, seq_length = seq_length)
    region_2_loc <- get_region_loc(start_in_seq = start_in_seq, end_in_seq = end_in_seq, region_number = 2, seq_length = seq_length)
    region_3_loc <- get_region_loc(start_in_seq = start_in_seq, end_in_seq = end_in_seq, region_number = 3, seq_length = seq_length)
    region_4_loc <- get_region_loc(start_in_seq = start_in_seq, end_in_seq = end_in_seq, region_number = 4, seq_length = seq_length)

    result <- list(
        region_1_loc = region_1_loc,
        region_2_loc = region_2_loc,
        region_3_loc = region_3_loc,
        region_4_loc = region_4_loc
    )
    return(result)
}

get_region_loc <- function(start_in_seq = NULL, end_in_seq = NULL, region_number = NULL, seq_length = NULL) {
    region_length <- 6
    left <- (start_in_seq - region_length * region_number):(start_in_seq - 1 - region_length * (region_number - 1))
    left <- left[left >= 1]
    if (length(left) > 0) {
        left <- c(left[1]:left[length(left)])
    } else {
        left <- NULL
    }

    right <- (end_in_seq + 1 + region_length * (region_number - 1)):(end_in_seq + region_length * region_number)
    right <- right[right <= seq_length]
    if (length(right) > 0) {
        right <- c(right[1]:right[length(right)])
    } else {
        right <- NULL
    }
    result <- c(left, right)
    return(result)
}

vector2range <- function(run = NULL) {
    if (!is.null(run)) {
        rundiff <- c(1, diff(run))
        difflist <- split(run, cumsum(rundiff != 1))
        result <- t(sapply(difflist, range))
        return(result)
    } else {
        return(NULL)
    }
}

block_bed_for_multiple_motifs <- function(motif_info_in_single_seq = NULL) {
    motif_collapsed <- collapse_overlap_motif(data = motif_info_in_single_seq)

    seq_info <- unlist(strsplit(motif_collapsed$seq_id[1], "_|-|:", perl = T))
    seq_start_in_genome <- as.numeric(seq_info[4])
    seq_end_in_genome <- as.numeric(seq_info[5])
    seq_length <- seq_end_in_genome - seq_start_in_genome + 1

    motif_count <- nrow(motif_collapsed)
    if (motif_count == 1) {
        result <- block_bed_for_single_motif(motif_info_in_single_seq = motif_collapsed)
        return(result)
    } else {
        motif_region_loc <- get_motif_region_loc_for_multiple_motifs(
            motif_collapsed = motif_collapsed,
            seq_length = seq_length
        )
        overlap_base_loc <- get_overlap_base_loc(motif_region_loc = motif_region_loc)
        motif_region_loc_except_overlap_base <- get_motif_region_loc_except_overlap_base(
            motif_region_loc = motif_region_loc,
            overlap_base_loc = overlap_base_loc
        )
        motif_region_loc_range_except_overlap_base <- lapply(motif_region_loc_except_overlap_base, function(x) {
            lapply(x, function(y) {
                if (length(y) > 0) {
                    a <- vector2range(run = y)
                    if (is.vector(a)) {
                        a[1] <- a[1] - 1
                    } else {
                        a[, 1] <- a[, 1] - 1
                    }
                    result <- a + seq_start_in_genome - 1
                    return(result)
                }
            })
        })
        return(motif_region_loc_range_except_overlap_base)
    }
}

get_motif_region_loc_for_multiple_motifs <- function(motif_collapsed = NULL, seq_length = seq_length) {
    motif_count <- nrow(motif_collapsed)
    result <- lapply(1:motif_count, function(counth) {
        one_motif_info <- motif_collapsed[counth, ]
        start_in_seq <- as.numeric(one_motif_info$start_in_seq)
        end_in_seq <- as.numeric(one_motif_info$end_in_seq)
        region_loc_for_one_motif <- get_region_loc_for_one_motif(start_in_seq = start_in_seq, end_in_seq = end_in_seq, seq_length = seq_length)
        result <- list(
            motif_loc = start_in_seq:end_in_seq,
            region_1_loc = region_loc_for_one_motif$region_1_loc,
            region_2_loc = region_loc_for_one_motif$region_2_loc,
            region_3_loc = region_loc_for_one_motif$region_3_loc,
            region_4_loc = region_loc_for_one_motif$region_4_loc
        )
        return(result)
    })
    return(result)
}

collapse_overlap_motif <- function(data = NULL) {
    first_overlap_index <- get_first_overlap_index(data = data)

    while (NROW(first_overlap_index) > 0) {
        remaining_row <- (1:nrow(data))[-unique(as.numeric(first_overlap_index))]
        a <- unique(sort(union(unlist(data[first_overlap_index[1], 2]):unlist(data[first_overlap_index[1], 3]), unlist(data[first_overlap_index[2], 2]):unlist(data[first_overlap_index[2], 3]))))
        collapsed_data <- data.frame(seq_id = data$seq_id[1], start_in_seq = a[1], end_in_seq = a[length(a)])
        data <- rbind(collapsed_data, data[remaining_row, ])
        first_overlap_index <- get_first_overlap_index(data = data)
    }

    return(data)
}

get_first_overlap_index <- function(data = NULL) {
    motif_count <- nrow(data)
    overlap_index <- c()
    if (motif_count > 1) {
        for (rowth1 in 1:(motif_count - 1)) {
            for (rowth2 in (rowth1 + 1):motif_count) {
                overlap_status <- length(intersect(
                    unlist(data[rowth1, 2]):unlist(data[rowth1, 3]),
                    unlist(data[rowth2, 2]):unlist(data[rowth2, 3])
                )) > 0
                if (overlap_status) {
                    overlap_index <- rbind(overlap_index, c(rowth1, rowth2))
                    break
                }
            }
            if (NROW(overlap_index) > 0) {
                break
            }
        }
    }
    return(overlap_index)
}

merge_overlap_index <- function(overlap_index = NULL) {
    overlap_count <- NROW(overlap_index)
    merged <- list()
    merge_element <- intersect(overlap_index[, 1], overlap_index[, 2])[1]
    merged <- list()
}

get_overlap_base_loc <- function(motif_region_loc = NULL) {
    a <- unlist(motif_region_loc)
    b <- table(a)
    overlap_base_loc <- as.numeric(names(b)[b > 1])
    return(overlap_base_loc)
}

get_motif_region_loc_except_overlap_base <- function(motif_region_loc = NULL, overlap_base_loc = NULL) {
    region_count <- length(motif_region_loc[[1]])
    result <- lapply(motif_region_loc, function(x) {
        a <- lapply(1:region_count, function(y) {
            if (y == 1) {
                return(x[[y]])
            } else {
                data <- x[[y]]
                data[!data %in% overlap_base_loc]
            }
        })
        names(a) <- names(motif_region_loc[[1]])
        return(a)
    })
    return(result)
}

write_bed <- function(motif_info_bed = NULL, chrom = NULL, output_dir = NULL) {
    seq_count <- length(motif_info_bed)
    for (i in 1:seq_count) {
        data <- motif_info_bed[[i]]
        only_1_motif <- !is.null(names(data))
        if (only_1_motif) {
            write_bed_for_one_motif(data = data, chrom = chrom, output_dir = output_dir)
        } else {
            motif_count <- length(data)
            for (j in 1:motif_count) {
                write_bed_for_one_motif(data = data[[j]], chrom = chrom, output_dir = output_dir)
            }
        }
    }
}
write_bed_for_one_motif <- function(data = NULL, chrom = NULL, output_dir = NULL) {
    motif_file <- paste0(output_dir, "/", chrom, "_collapsed_motif.bed")
    region1_file <- paste0(output_dir, "/", chrom, "_region1.bed")
    region2_file <- paste0(output_dir, "/", chrom, "_region2.bed")
    region3_file <- paste0(output_dir, "/", chrom, "_region3.bed")
    region4_file <- paste0(output_dir, "/", chrom, "_region4.bed")


    result <- lapply(data, function(x) {
        if (!is.null(x)) {
            if (is.vector(x)) {
                return(c(chrom, x))
            } else {
                return(cbind(chrom, x))
            }
        }
    })
    write_file(data = result$motif_loc, filename = motif_file)
    write_file(data = result$region_1_loc, filename = region1_file)
    write_file(data = result$region_2_loc, filename = region2_file)
    write_file(data = result$region_3_loc, filename = region3_file)
    write_file(data = result$region_4_loc, filename = region4_file)
}

write_file <- function(data = NULL, filename = NULL) {
    if (!is.null(data)) {
        write.table(matrix(data, ncol = 3),
            file = filename,
            quote = F, row.names = F, col.names = F, sep = "\t", append = T
        )
    }
}
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
motif_bed_file <- Args[6]
output_dir <- Args[7]


# 4. variable setting of test module---------------------------------------

# motif_bed_file<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_2/motif_loc/Pf3D7_01_v3.motif.bed"
# output_dir<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/region_loc"
# 5. process --------------------------------------------------------------
chrom <- gsub(".motif.bed", "", basename(motif_bed_file))

motif_file <- paste0(output_dir, "/", chrom, "_collapsed_motif.bed")
region1_file <- paste0(output_dir, "/", chrom, "_region1.bed")
region2_file <- paste0(output_dir, "/", chrom, "_region2.bed")
region3_file <- paste0(output_dir, "/", chrom, "_region3.bed")
region4_file <- paste0(output_dir, "/", chrom, "_region4.bed")
sapply(c(motif_file, region1_file, region2_file, region3_file, region4_file), remove_file)


motif_bed <- fread(motif_bed_file, header = F, stringsAsFactors = F, select = c(1:3))
colnames(motif_bed) <- c("seq_id", "start_in_seq", "end_in_seq")
motif_bed$start_in_seq <- motif_bed$start_in_seq + 1
motif_bed <- unique(motif_bed)

motif_info <- split(motif_bed, f = motif_bed$seq_id, drop = T)
motif_info_bed <- lapply(motif_info, function(x) {
    extract_block_loc(motif_info_in_single_seq = x)
})

# for (x in 1:length(motif_info)) {
#     print(x)
#     a <- extract_block_loc(motif_info_in_single_seq = motif_info[[x]])
# }
write_bed(motif_info_bed = motif_info_bed, chrom = chrom, output_dir = output_dir)
