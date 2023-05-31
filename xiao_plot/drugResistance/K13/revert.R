#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
filename <- Args[6]
outfilename <- Args[7]
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read.table(filename, header = F, as.is = T, sep = "\t")
data <- data[nrow(data):1, ]
write.table(data, outfilename, row.names = F, col.names = F, sep = "\t", quote = F)
