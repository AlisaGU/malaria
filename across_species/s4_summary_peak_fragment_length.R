#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_reichenowi/peak_nopeak")
peak_fragment_len <- unlist(fread("awk '{print $3-$2}' *.peak.bed|sort -n"))
nopeak_fragment_len <- unlist(fread("awk '{print $3-$2}' *.nopeak.bed|sort -n"))

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
hist(peak_fragment_len)
hist(nopeak_fragment_len)
