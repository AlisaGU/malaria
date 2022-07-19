#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(plyr)
# 2. functions ------------------------------------------------------------ TODO:
only_retain_nonoverlap_motifs <- function(data = NULL) {
    marker_region <- unlist(data[1, ])
    if (NROW(data) == 1) {
        result <- data
    } else {
        result <- data[1, ]
        data <- data[-1, ]
        while (NROW(data) > 0) {
            first_non_overlapped_index <- which(!apply(data, 1, function(x) {
                length(intersect(
                    (as.numeric(marker_region[2]) + 1):as.numeric(marker_region[3]),
                    (as.numeric(x[2]) + 1):as.numeric(x[3])
                )) > 0
            }))[1]

            if (!is.na(first_non_overlapped_index)) {
                new_region <- data[first_non_overlapped_index, ]
                result <- rbind(result, new_region)
                data <- data[-c(1:first_non_overlapped_index), ]
                marker_region <- unlist(new_region)
            } else {
                break
            }
        }
    }
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:
Args <- commandArgs()
cmd <- Args[6]

# 5. process -------------------------------------------------------------- TODO:
data <- fread(cmd = cmd, stringsAsFactors = F)

result <- ddply(data, .(V1), only_retain_nonoverlap_motifs)
write.table(result, "", quote = F, sep = "\t", row.names = F, col.names = F)
