#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------


# 2. functions ------------------------------------------------------------
get_p_value <- function(ancestor_region_combination = NULL) {
    allvar <- get_p_value_for_one_variant(combin = ancestor_region_combination, var_type = "allvar")
    snps <- get_p_value_for_one_variant(combin = ancestor_region_combination, var_type = "snps")
    indels <- get_p_value_for_one_variant(combin = ancestor_region_combination, var_type = "indels")

    result <- list(allvar = allvar, snps = snps, indels = indels)
    return(result)
}
get_p_value_for_one_variant <- function(combin = NULL, var_type = NULL) {
    data <- read.table(paste0(combin, "/", var_type, "_mean_variant_count"),
        stringsAsFactors = F, header = F
    )
    colnames(data) <- c(
        "nopeak", paste("relative", c(0.01, 0.05, 0.1, 0.2, 0.5), sep = ""),
        paste("absolute", c(10, 30, 50, 100, 200), sep = "")
    )
    rownames(data) <- paste("chrom", 1:14, sep = "")

    relative <- c()
    for (i in 2:5) {
        for (j in (i + 1):6) {
            p <- wilcox.test(data[, i], data[, j], paired = T, alternative = "greater")$p.value
            relative <- rbind(relative, c(colnames(data)[i], colnames(data)[j], p))
        }
    }

    absolute <- c()
    for (i in 7:10) {
        for (j in (i + 1):11) {
            p <- wilcox.test(data[, i], data[, j], paired = T, alternative = "greater")$p.value
            absolute <- rbind(absolute, c(colnames(data)[i], colnames(data)[j], p))
        }
    }
    result <- list(relative = relative, absolute = absolute)
    return(result)
}

# 3. input ----------------------------------------------------------------
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/WSEA")

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
whole_chromosome_3D7 <- get_p_value(ancestor_region_combination = "mutation_density_3D7_whole_chromosome")
two_outgroup_two_outgroup <- get_p_value(ancestor_region_combination = "mutation_density_two_group_consistent_two_group_consistent")
