#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(seqinr)
# 2. functions ------------------------------------------------------------ TODO:

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno")
gff <- fread(cmd = "grep -v \"^##\" PlasmoDB-36_Pfalciparum3D7.gff", stringsAsFactors = F, sep = "\t")
cds_fasta <- read.fasta(file = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_AnnotatedCDSs.fasta", as.string = TRUE, seqtype = "DNA", set.attributes = F)
degeneracy <- read.table("degeneracy.txt", sep = "\t", header = T, as.is = T)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(gff) <- c("sequenceID", "source", "featureType", "featureStart", "featureEnd", "score", "strand", "phase", "attributes")
degeneracy$type <- paste(degeneracy$aminoAcid, degeneracy$codon_index, sep = "_")

gene <- gff[gff$featureType == "gene", ]
pseudo_gene <- gene[grepl("pseudogene", gene$attributes), ]
pseudo_gene_ID <- sapply(pseudo_gene$attributes, function(x) {
    a <- unlist(strsplit(x, ";"))[1]
    b <- unlist(strsplit(a, "="))[2]
})

cds <- gff[gff$featureType == "CDS", ]
mrna_csdOrder <- lapply(cds$attributes, function(x) {
    a <- unlist(strsplit(x, ";"))[1]
    b <- unlist(strsplit(a, "="))[2]
    d <- unlist(strsplit(b, "-"))
    return(c(paste0(d[1], "-", d[2]), gsub("CDS", "", d[3])))
})
mrna_csdOrder <- do.call(rbind, mrna_csdOrder)
cds$mrna <- mrna_csdOrder[, 1]
cds$cdsOrder <- as.numeric(mrna_csdOrder[, 2])

mrna_dataset <- split(cds, f = cds$mrna)
mrna_dataset_true_gene <- mrna_dataset[-sapply(pseudo_gene_ID, function(x) {
    which(grepl(x, names(mrna_dataset)))
})]
genomePos_codon_index <- lapply(mrna_dataset_true_gene, function(x) {
    chrom <- as.character(unlist(x[1, "sequenceID"]))
    strand <- as.character(unlist(x[1, "strand"]))
    mRNA <- as.character(unlist(strsplit(unlist(strsplit(
        unlist(x[1, "attributes"]),
        ";"
    ))[1], "=|-", perl = TRUE))[2])
    mrna_seq <- cds_fasta[mRNA]
    pep_seq <- translate(s2c(unlist(mrna_seq)))
    pep_seq_cds_each_site <- rep(pep_seq, each = 3)

    seq_num <- lapply(1:nrow(x), function(i) {
        y <- x[i, ]
        seq(as.numeric(unlist(y[, "featureStart"])), as.numeric(unlist(y[, "featureEnd"])))
    })
    seq_num <- do.call(c, seq_num)
    if (strand == "-") {
        seq_num <- rev(seq_num)
    }
    result <- data.frame(chrom = chrom, genomePos = seq_num, codon_index = rep(1:3, length(seq_num) / 3), pep = pep_seq_cds_each_site)
    result$type <- paste(result$pep, result$codon_index, sep = "_")
    result$degeneracy <- degeneracy$degeneracy[match(result$type, degeneracy$type)]
    return(result)
})
genomePos_codon_index <- do.call(rbind, genomePos_codon_index)
write.table(genomePos_codon_index, "cds_each_site_codon_index.txt", row.names = F, col.names = F, quote = F, sep = "\t")
