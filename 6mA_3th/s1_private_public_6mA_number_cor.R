#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
private_sliding_window <- read.table("private/sliding_window/sliding_window_6mA_distri", header = F, as.is = T)
private_sliding_window_6mA_count <- private_sliding_window[, seq(1, ncol(private_sliding_window), by = 2)]
public_sliding_window <- read.table("public/sliding_window/sliding_window_6mA_distri", header = F, as.is = T)
public_sliding_window_6mA_count <- public_sliding_window[, seq(1, ncol(public_sliding_window), by = 2)]
plot(unlist(private_sliding_window_6mA_count), unlist(public_sliding_window_6mA_count))
cor.test(unlist(private_sliding_window_6mA_count), unlist(public_sliding_window_6mA_count))

private_motif_flank <- read.table("private/motif_flank_pattern/motif_flank_pattern_6mA_distri", header = F, as.is = T)
private_motif_flank_6mA_count <- private_motif_flank[, seq(1, ncol(private_motif_flank), by = 2)]
public_motif_flank <- read.table("public/motif_flank_pattern/motif_flank_pattern_6mA_distri", header = F, as.is = T)
public_motif_flank_6mA_count <- public_motif_flank[, seq(1, ncol(public_motif_flank), by = 2)]
plot(unlist(private_motif_flank_6mA_count), unlist(public_motif_flank_6mA_count))
cor.test(unlist(private_motif_flank_6mA_count), unlist(public_motif_flank_6mA_count))
