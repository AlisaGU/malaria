#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
## 这是探究motif作用机制时，统计的各个区间（control、submit~5%···）的motif平均数目

window_5_percent_peak_length <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif_GAWGAW/motif_count/GAWGAW_count_sliding_window", header = F, as.is = T)
colnames(window_5_percent_peak_length) <- c("control", paste("w", 1:10, sep = ""))
rownames(window_5_percent_peak_length) <- paste0("chrom", 1:14, sep = "")
apply(window_5_percent_peak_length, 2, mean)

## 这是用motif预测信号区时，统计的信号区和非信号区的motif平均长度


## 这是用motif预测信号区时，统计的信号区和非信号区的motif平均数目
