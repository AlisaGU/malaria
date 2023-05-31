#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/offtarget/BAM")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
WT <- read.table("3D7-T4_minus_log2ratio_bin5000.bdg", header = F, as.is = T)
KD <- read.table("6mAKD-T4_minus_log2ratio_bin5000.bdg", header = F, as.is = T)
colnames(WT) <- colnames(KD) <- c("chrom", "start", "end", "height")

data <- inner_join(WT, KD, by = c("chrom", "start", "end"))
data$minus <- data$height.y - data$height.x

write.table(data[, c(1, 2, 3, 6)], "KD-WT.bed", quote = F, row.names = F, col.names = F, sep = "\t")
