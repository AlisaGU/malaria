#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(tidyr)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
tree_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/tree/RAxML_bipartitions.tree.OutGroup.remove2.mafft1.raxml.tree"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
tree <- read.newick(tree_filename)
tree <- groupClade(tree, .node = 31)

d <- data.frame(
    label = tree$tip.label,
    species = c(
        "P.berghei", "P. yoelii yoelii", "P. falciparum", "P. knowlesi",
        "P.cynomolgi", "P. vivax", "T. gondii", "N. gaditana",
        "A.castellanii", "C. reinhardtii", "T. thermophila"
    ),
    gene = c(
        "PBANKA_1363500", "YM PYYM_1365600", "PF3D7_1350700 (PfNAMT)", "strain H PKNH_1250200",
        "PCYB_121300", "P01 PVP01_1204100", "GT1 TGGT1_272210", "CCMP526 XM_005855954",
        "ELR14023.1 MTA70 family", "MT-A70-like"
    )
)
ggtree(tree) +
    geom_tiplab(size = 5) +
    geom_nodepoint(aes(color = as.numeric(label)), size = 3) +
    # geom_text2(aes(subset = !isTip, label = node),color="red", hjust = -.3)+
    # geom_text2(aes(subset = isTip, label = node),color="blue", hjust = -.3)
    scale_color_gradient(low = "#c67f7f", high = "#a60026", name = "Bootstrap value") +
    geom_treescale()
# geom_nodelab(color = "red")+


viewClade(ggtree(tree, size = 2) + geom_tiplab(), node = 31)
