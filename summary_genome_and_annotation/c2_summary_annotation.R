#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno")
plasmo_anno <- fread(cmd = "zgrep -v \"^#\" PlasmoDB-36_Pfalciparum3D7.gff.gz", header = F, stringsAsFactors = F, sep = "\t")
ncbi_anno <- fread(cmd = "zgrep -v \"^#\" GCF_000002765.5_GCA_000002765_genomic.gff.gz", header = F, stringsAsFactors = F, sep = "\t")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:

colnames(plasmo_anno) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
colnames(ncbi_anno) <- colnames(plasmo_anno)

table(plasmo_anno$feature)

# 1. plasmodb中一个gene可以有多个mRNA
gene_index <- which(plasmo_anno$feature == "gene")
gene_count <- length(gene_index)
search_index <- c(gene_index, nrow(plasmo_anno) + 1)
feature_count_of_gene <- sapply(1:gene_count, function(i) {
    subfeature <- unlist(plasmo_anno[(search_index[i] + 1):(search_index[i + 1] - 1), "feature"])
    subfeature_count <- sum(sapply(subfeature, function(x) {
        # candi_feature <- c("mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA")
        # candi_feature <- c("ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA")
        candi_feature <- "mRNA"
        ifelse(x %in% candi_feature, 1, 0)
    }))
    return(subfeature_count)
})

# plasmodb 和 ncbi 的id  match情况
plasmo_gene <- plasmo_anno[gene_index, ]
plasmo_id <- sapply(plasmo_gene$attribute, function(x) {
    k <- unlist(strsplit(unlist(x), ";", fixed = T))[1]
    gsub("ID=", "", k)
})
names(plasmo_id) <- NULL

ncbi_gene <- ncbi_anno[ncbi_anno$feature == "gene", ]
ncbi_id <- sapply(ncbi_gene$attribute, function(x) {
    k <- unlist(strsplit(unlist(x), ";", fixed = T))[1]
    gsub("ID=gene-", "", k)
})
names(ncbi_id) <- NULL

length(intersect(plasmo_id, ncbi_id))
length(plasmo_id)
length(ncbi_id)