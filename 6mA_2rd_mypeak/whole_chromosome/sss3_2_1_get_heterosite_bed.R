#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(parallel)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
inputFile <- Args[6]
output <- paste0(dirname(inputFile), "/", gsub(".vcf.gz", "", basename(inputFile)), "_polysite.bed")
# 4. variable setting of test module--------------------------------------- TODO:
# inputFile <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_PASS_01.vcf.gz"

# 5. process -------------------------------------------------------------- TODO:
GT <- fread(cmd = paste0("zcat ", inputFile, " |grep -v \"^#\"|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");printf(\"%s%s\",a[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F)
GT[, ncol(GT) := NULL]
gt <- transpose(GT[, 3:ncol(GT)]) # only the first two columns are about position.
alt_allele_class_count <- mclapply(gt, function(x) {
    allele_class <- unique(unlist(strsplit(unlist(x), "/|\\/")))
    allele_class_demissing <- allele_class[allele_class != "."]
    return(length(unique(allele_class_demissing)) - 1)
}, mc.cores = 4)
alt_allele_class_count <- do.call("c", alt_allele_class_count)

allele_class_count_lt1_index <- which(alt_allele_class_count != 0)
result <- data.frame(
    GT[allele_class_count_lt1_index, 1],
    GT[allele_class_count_lt1_index, 2] - 1,
    GT[allele_class_count_lt1_index, 2],
    alt_allele_class_count[allele_class_count_lt1_index]
)
write.table(result, output, quote = F, sep = "\t", row.names = F, col.names = F)