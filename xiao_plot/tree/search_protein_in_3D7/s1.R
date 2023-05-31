#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/search_protein_in_3D7")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
proteinPattern.pp <- read.table("proteinPattern.pp.bed", header = F, as.is = T)
proteinPattern.gxg <- read.table("proteinPattern.gxg.bed", header = F, as.is = T)

common_seqname <- intersect(unique(proteinPattern.pp[, 1]), unique(proteinPattern.gxg[, 1]))
pp <- proteinPattern.pp[match(common_seqname, proteinPattern.pp[, 1]), ]
gxg <- proteinPattern.gxg[match(common_seqname, proteinPattern.gxg[, 1]), ]
pp_gxg <- cbind(pp, gxg)
pp_gxg <- pp_gxg[, c(1, 2, 3, 4, 8, 9, 10)]
colnames(pp_gxg) <- c("proteinName", "PP_startInProtein", "PP_endInProtein", "PP", "GXG_startInProtein", "GXG_endInProtein", "GXG")
pp_gxg$PP_startInProtein <- pp_gxg$PP_startInProtein + 1
pp_gxg$GXG_startInProtein <- pp_gxg$GXG_startInProtein + 1
write.table(pp_gxg, "pp_gxg.txt", row.names = F, col.names = T, sep = "\t", quote = F)
