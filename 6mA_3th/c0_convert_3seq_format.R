#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/private")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read.csv("3D7-6mA-3seq_noheader.csv", as.is = T, header = F)
data <- data[1:(nrow(data) - 1), ]
colnames(data) <- c(
    "chr", "modify_type", "start", "end",
    "strand", "coverage", "context", "frac",
    "gene_id", "gene_name", "gene_note",
    "position1", "position2"
)
result <- data.frame(
    Seqid = data$chr, Start = data$start - 1, End = data$end, Score = ".", Strand = data$strand,
    Phase = ".", coverage = data$coverage, context = data$context, IPDRatio = ".", frac = data$frac,
    fracLow = ".", fracUp = ".", identificationQV = "."
)
write.table(result, "3D7-6mA-3seq.bed", sep = "\t", row.names = F, col.names = F, quote = F)


# compare between two files
data_private <- data
data_public <- data.table::fread("../public/3D7-6mA-9Cell.txt", header = T, stringsAsFactors = F)

site_private <- paste(data_private$chr, data_private$start)
site_public <- paste(data_public$Seqid, data_public$Start)

length(site_private)
length(site_public)
length(intersect(site_private, site_public))
