#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)

# 2. functions ------------------------------------------------------------
get_p_distance_parameter <- function(inputFile = NULL) {
    GT <- fread(inputFile, sep = "\t", stringsAsFactors = F, header = F)
    GT[, ncol(GT) := NULL]
    result <- sapply(GT, function(x) {
        base_summary <- as.data.frame(table(x))
        total_base_count <- sum(base_summary[base_summary$x != ".", 2])
        diff_base_count <- sum(base_summary[base_summary$x != "." & base_summary$x != "0", 2])
        return(c(diff_base_count, total_base_count))
    })
    rownames(result) <- c("diff_base_count", "total_base_count")
    return(result)
}


# 3. input ----------------------------------------------------------------
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/OCE/africa_biallelic_variant")

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
files <- dir()
vcf_files <- files[!grepl("vcf.gz", files)]

p_dist_parameter <- list()
for (i in 1:14) {
    inputFile <- vcf_files[[i]]
    p_dist_parameter[[i]] <- get_p_distance_parameter(inputFile = inputFile)
}
names(p_dist_parameter) <- inputFile

combined_p_dist <- matrix(data = 0, nrow = 2, ncol = 6628)
for (i in 1:14) {
    combined_p_dist <- combined_p_dist + p_dist_parameter[[i]]
}
p_dist <- combined_p_dist[1, ] / combined_p_dist[2, ]

p_dist_max_index <- which(p_dist == max(p_dist)) # 2048
