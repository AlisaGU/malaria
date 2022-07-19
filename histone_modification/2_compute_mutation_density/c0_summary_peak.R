#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
peak_info_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/histone_methy/3D7C8-H3K36me3-40h/BAM_broad_nomodel_dup1/3D7C8-H3K36me3-40h_peaks.bed"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
peak_info <- read.table(peak_info_filename, header = F, as.is = T)
peak_length <- peak_info[, 3] - peak_info[, 2] + 1
summary(peak_length)
hist(peak_length)
