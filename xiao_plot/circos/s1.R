#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(data.table)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
filename <- Args[6]

# 4. variable setting of test module--------------------------------------- TODO:
# filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult1/mutation_density_distribution/wholeGenome/ESEA_WSEA_OCE_SAM_SAS/variant_3D7/all_chrom_bin5000_tmp.bed"

# 5. process -------------------------------------------------------------- TODO:
data <- fread(filename, header = F, stringsAsFactors = F)
snp_data <- data %>%
    group_by(V1, V2, V3) %>%
    summarise(SNP = sum(V9))

indel_data <- data %>%
    group_by(V1, V2, V3) %>%
    summarise(SNP = sum(V12))
write.table(snp_data, gsub("tmp", "snp", filename), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(indel_data, gsub("tmp", "indel", filename), quote = F, row.names = F, col.names = F, sep = "\t")
