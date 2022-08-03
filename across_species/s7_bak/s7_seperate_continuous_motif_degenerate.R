#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(plyr)

# 2. functions ------------------------------------------------------------ TODO:
get_nonoverlap_motif <- function(data_overlap = NULL) {
    if (nrow(data_overlap) == 1) {
        return(data_overlap)
    } else if (nrow(data_overlap) == 2) {
        data_overlap_order <- data_overlap[order(data_overlap$sPos_minus_1, decreasing = F), ]
        return(data_overlap_order[1, ])
    } else if (nrow(data_overlap) == 3) {
        data_overlap_order <- data_overlap[order(data_overlap$sPos_minus_1, decreasing = F), ]
        return(data_overlap_order[c(1, 3), ])
    } else if (nrow(data_overlap) >= 4) {
        result <- related_motif_count(data_overlap = data_overlap)
        return(result)
    }
}

related_motif_count <- function(data_overlap = NULL) {
    data_overlap_order <- data_overlap[order(data_overlap$sPos_minus_1, decreasing = F), ]

    retain_data <- data_overlap_order[c(1, nrow(data_overlap_order)), ]
    first_nonoverlap_index <- 1
    last_nonoverlap_index <- nrow(data_overlap_order)

    first_marker_region <- unlist(data_overlap_order[first_nonoverlap_index, ])
    last_marker_region <- unlist(data_overlap_order[last_nonoverlap_index, ])
    data_overlap_order <- data_overlap_order[-c(first_nonoverlap_index, last_nonoverlap_index), ]

    while (nrow(data_overlap_order) > 0 & first_nonoverlap_index < last_nonoverlap_index) {
        first_nonoverlap_index <- which(!apply(data_overlap_order, 1, function(x) {
            length(intersect(
                (as.numeric(first_marker_region[2]) + 1):as.numeric(first_marker_region[3]),
                (as.numeric(x[2]) + 1):as.numeric(x[3])
            )) > 0
        }))[1]
        last_nonoverlap_index <- rev(which(!apply(data_overlap_order, 1, function(x) {
            length(intersect(
                (as.numeric(last_marker_region[2]) + 1):as.numeric(last_marker_region[3]),
                (as.numeric(x[2]) + 1):as.numeric(x[3])
            )) > 0
        })))[1]
        overlap_status <- length(intersect(
            (unlist(data_overlap_order[first_nonoverlap_index, 2]) + 1):unlist(data_overlap_order[first_nonoverlap_index, 3]),
            (unlist(data_overlap_order[last_nonoverlap_index, 2]) + 1):unlist(data_overlap_order[last_nonoverlap_index, 3])
        ))
        if (first_nonoverlap_index < last_nonoverlap_index & overlap_status != 0) {
            retain_data <- rbind(retain_data, data_overlap_order[first_nonoverlap_index, ])
            break
        } else if (first_nonoverlap_index < last_nonoverlap_index & overlap_status == 0) {
            retain_data <- rbind(retain_data, unique(data_overlap_order[c(first_nonoverlap_index, last_nonoverlap_index), ]))
            first_marker_region <- unlist(data_overlap_order[first_nonoverlap_index, ])
            last_marker_region <- unlist(data_overlap_order[last_nonoverlap_index, ])
            data_overlap_order <- data_overlap_order[-c(1:first_nonoverlap_index, last_nonoverlap_index:nrow(data_overlap_order)), ]
            if (nrow(data_overlap_order) == 1) {
                break
            }
        } else if (first_nonoverlap_index == last_nonoverlap_index) {
            retain_data <- rbind(retain_data, unique(data_overlap_order[c(first_nonoverlap_index, last_nonoverlap_index), ]))
        } else {
            break
        }
    }
    result <- retain_data[order(retain_data$sPos_minus_1, decreasing = F), ]
    return(result)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
working_dir <- Args[6]
cmd <- Args[7]
result_name <- Args[8]
# 4. variable setting of test module--------------------------------------- TODO:
# working_dir<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/all_motifs"
#   cmd<-"/picb/evolgen/users/gushanshan/software/bedtools/bedtools merge -i motif_bed/Pf3D7_01_v3.motif.peak.overlap_unique.bed | awk '{print $0\"\\t\"NR}' | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a motif_bed/Pf3D7_01_v3.motif.peak.overlap_unique.bed -b stdin -wb"

# 5. process -------------------------------------------------------------- TODO:
setwd(working_dir)
data <- fread(cmd = cmd, header = F, stringsAsFactors = F, drop = 6:8)
colnames(data) <- c("chrom", "sPos_minus_1", "ePos", "strand", "motif", "overlap_class")
result <- ddply(data, .(overlap_class), get_nonoverlap_motif)
write.table(result, result_name, quote = F, sep = "\t", row.names = F, col.names = F)
