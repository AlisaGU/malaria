#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf")
data <- fread(cmd = "grep \"^s P_falciparum.Pf3D7_04_v3\" Pf3D7_04_v3.P_falciparum.P_reichenowi.P_praefalciparum.maf|awk '{print $3,$4}'")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

colnames(data) <- c("Start_minus_1", "block_character_len")
data$Start <- data$Start_minus_1 + 1
data$End <- data$Start_minus_1 + data$block_character_len

gap_len <- c(unlist(data[1, "Start_minus_1"]), sapply(2:nrow(data), function(i) {
    unlist(data[i, "Start"] - data[i - 1, "End"] - 1)
}))
