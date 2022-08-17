#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(parallel)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ----------------------------------------------------------------
Args <- commandArgs()
vcf_data <- Args[6]

output <- paste0(dirname(vcf_data), "/", gsub(".vcf.gz", "", basename(vcf_data)), "_polysite_nomatter_ancestor.bed")

bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/CAF/variant/CAF_PASS_01.vcf.gz"

# 5. process --------------------------------------------------------------
vcf <- fread(cmd = paste0(bcftools, " view -H ", vcf_data, " |awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s%s%s\",b[1],OFS,b[2],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)

vcf[, ncol(vcf) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", paste("popu", 1:(ncol(vcf) - 4), sep = ""))
chrom <- vcf$chrom[1]
data_for_mutation_load <- vcf

gt <- transpose(data_for_mutation_load[, 3:ncol(data_for_mutation_load)])
mutation_load <- mclapply(gt, function(x) {
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)

    allele_class <- unique(x[3:length(x)])
    allele_class_nomiss <- allele_class[allele_class != "."]
    exist_allele <- allele[allele_class_nomiss]
    exist_allele_no_spanning_deletion <- exist_allele[exist_allele != "*"] # 排除*和.之后，群体中还存在的allele类型


    allele_base_count <- sapply(exist_allele_no_spanning_deletion, nchar)
    result <- c()
    if (all(allele_base_count > 1)) {
        result <- c(0, length(exist_allele_no_spanning_deletion) - 1)
    } else if (all(allele_base_count == 1)) {
        result <- c(length(exist_allele_no_spanning_deletion) - 1, 0)
    } else {
        length_equal_1 <- length(which(allele_base_count == 1))
        length_larger_than_1 <- length(which(allele_base_count > 1))
        result <- c(length_equal_1 - 1, length_larger_than_1)
    }
    names(result) <- c("snps", "indels")
    return(result)
}, mc.cores = 4)
mutation_load <- do.call(rbind, mutation_load)
result <- data.frame(chrom = chrom, bed0 = data_for_mutation_load$pos - 1, bed1 = data_for_mutation_load$pos, mutation_load = mutation_load)
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
