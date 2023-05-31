#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/tree.OutGroup.remove2.add_PVC50508.1/domain_tree/remove_outgroup")

header <- read.table("domain.seq.1.1.removeOutgroup.mafft.mega.fas.header", header = FALSE, sep = ">")[, 2]
header <- gsub(" ", "_", header)

domains <- read.table("../domain_index_in_protein.txt", header = T, as.is = T, sep = "\t")
domains$position <- paste(domains$protein, paste(domains$start_of_domain_in_protein, domains$end_of_domain_in_protein, sep = "-"), sep = "_")
domains1 <- domains[match(header, domains$position), ]
