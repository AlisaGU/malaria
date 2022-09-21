#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:
read_degeneracy <- function(filename = NULL) {
    degeneracy <- read.table(filename, as.is = T, header = F)
    colnames(degeneracy) <- c(
        "6mA_count_in_degeneracy_1", "degeneracy_1_totol_count",
        "6mA_count_in_degeneracy_2", "degeneracy_2_totol_count",
        "6mA_count_in_degeneracy_3", "degeneracy_3_totol_count",
        "6mA_count_in_degeneracy_4", "degeneracy_4_totol_count"
    )
    proportion <- degeneracy[, c(1, 3, 5, 7)] / degeneracy[, c(2, 4, 6, 8)]
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/motif_start_site")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
