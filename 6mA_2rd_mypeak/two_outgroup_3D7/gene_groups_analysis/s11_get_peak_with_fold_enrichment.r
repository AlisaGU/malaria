#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- fread(cmd = paste0("grep -v \"#\" 3D7-T3_peaks.xls|grep -v \"chr\"|grep -v \"^$\"|awk '{print $1\"\t\"$2-1\"\t\"$3\"\t\"$8}'"))
result <- as.data.frame(data %>% group_by(V1, V2, V3) %>% filter(V4 == max(V4)) %>% slice(1))
write.table(result, file = "peak_with_fold_enrichment.bed", sep = "\t", quote = F, row.names = F, col.names = F)
