#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(tidyr)
library(seqinr)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
tree <- Args[6]
seq <- Args[7]
outfile <- gsub(".fasta", ".orderedByTree.fasta", seq)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

tr <- read.newick(tree)
p <- ggtree(tr, size = 1) +
    geom_tiplab() +
    geom_treescale(color = "red", fontsize = 5)

taxa_order <- get_taxa_name(p)
sequence <- read.fasta(file = seq, as.string = TRUE, seqtype = "AA", set.attributes = F)
result_seq <- sequence[match(taxa_order, names(sequence))]
write.fasta(result_seq, names = names(result_seq), file.out = outfile)
pdf(paste0(tree, ".pdf"), width = 15, height = 10)
p
dev.off()
