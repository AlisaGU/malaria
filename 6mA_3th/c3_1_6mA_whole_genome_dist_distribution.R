#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(ggplot2)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
bedtools <- "/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
m6mA_bed <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/3D7-6mA-9Cell_FRACgt20_COVgt25.bed"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/m6mA_dist_distri")
# cmd <- paste0(bedtools, " closest -d -io -t first -a ", m6mA_bed, " -b ", m6mA_bed, " | awk '{print $NF}'")
cmd <- paste0(bedtools, " closest -d -io -t first -a ", m6mA_bed, " -b ", m6mA_bed, " | awk '{print $NF}'")

dist <- unlist(fread(cmd = cmd, stringsAsFactors = FALSE))
freq <- as.data.frame(table(dist))
global <- ggplot(data = freq, aes(x = dist, y = Freq)) +
    geom_line(aes(group = 1)) +
    theme(axis.text.x = element_blank())
subset_1_100 <- ggplot(data = freq, aes(x = dist, y = Freq)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    scale_x_discrete(limits = factor(seq(1, 100))) +
    theme(axis.text.x = element_text(angle = 30))

subset_1_20 <- ggplot(data = freq, aes(x = dist, y = Freq)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    scale_x_discrete(limits = factor(seq(1, 20)))
global / subset_1_100 / subset_1_20
ggsave("m6mA_dist_distri.pdf", width = 10, height = 21)
