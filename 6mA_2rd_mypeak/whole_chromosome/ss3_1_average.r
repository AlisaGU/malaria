#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------


# 2. functions ------------------------------------------------------------


# 3. input ----------------------------------------------------------------
Args <- commandArgs()
filename <- Args[6]

# 4. variable setting of test module---------------------------------------
# filename<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/sliding_window_noCollapseWindow/relative.0.05/Pf3D7_01_v3/w1.mean.txt"

# 5. process --------------------------------------------------------------
data <- read.table(filename, header = F, stringsAsFactors = F)
cat(mean(unlist(data)))
