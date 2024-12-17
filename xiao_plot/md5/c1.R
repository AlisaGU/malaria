#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/md5")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
file1 <- read.table("file1.txt", stringsAsFactors = F, header = F)
file2 <- read.table("file2.txt", stringsAsFactors = F, header = F)
md5sum <- read.table("md5_total", stringsAsFactors = F, header = F)
colnames(md5sum) <- c("md5", "filename")
file1_md5 <- md5sum[match(unlist(file1), md5sum$filename), ]
file2_md5 <- md5sum[match(unlist(file2), md5sum$filename), ]
write.table(file1_md5, "file1_md5.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(file2_md5, "file2_md5.txt", row.names = F, col.names = F, sep = "\t", quote = F)
