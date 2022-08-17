#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(seqinr)
library(stringr)
library(data.table)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf")

Args <- commandArgs()
vcf_name <- Args[6]
raw_vcf <- fread(vcf_name, stringsAsFactors = F, header = T)
align_fasta <- read.fasta(Args[7], seqtype = "DNA", as.string = T)
# 4. variable setting of test module--------------------------------------- TODO:
# raw_vcf <- fread("Pf3D7_01_v3.raw.vcf", stringsAsFactors = F, header = T)
# align_fasta <- read.fasta("Pf3D7_01_v3.P_falciparum.P_reichenowi.P_praefalciparum.fasta", seqtype = "DNA", as.string = T)

# 5. process -------------------------------------------------------------- TODO:
ref_fasta <- align_fasta["P_falciparum"]
gap_in_ref <- str_locate_all(ref_fasta, "-")[[1]][, 1]
a <- raw_vcf$POS
b <- sapply(1:length(a), function(x) {
    a[x] - length(which(gap_in_ref <= a[x]))
})
raw_vcf$POS <- b
write.table(raw_vcf, file = gsub(".raw", "", vcf_name), sep = "\t", quote = F, row.names = F, col.names = F)
