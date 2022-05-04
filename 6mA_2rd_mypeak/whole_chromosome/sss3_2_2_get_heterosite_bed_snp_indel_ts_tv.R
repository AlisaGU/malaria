#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)
# 2. functions ------------------------------------------------------------ TODO:
get_types <- function(allele = NULL) {
    ref <- allele[allele != "*"][1]
    ref_index <- which(allele == ref)
    alt_index <- (1:length(allele))[-ref_index]

    alt_types <- ""
    if (nchar(ref) == 1) {
        alt_types <- sapply(allele[alt_index], function(allele.i) {
            get_type_for_one_base(ref = ref, alt = allele.i)
        })
    } else {
        alt_types <- sapply(allele[alt_index], function(allele.i) {
            if (nchar(allele.i) == nchar(ref)) {
                get_type_for_one_base(ref = substr(ref, 1, 1), alt = substr(allele.i, 1, 1))
            } else if (nchar(allele.i) != nchar(ref)) {
                return("indels")
            }
        })
    }

    return(c("ref", alt_types))
}

get_type_for_one_base <- function(ref = NULL, alt = NULL) {
    if (alt == "*" | nchar(alt) > 1) {
        return("indels")
    } else if ((ref == "A" & alt == "G") | (ref == "G" & alt == "A") | (ref == "T" & alt == "C") | (ref == "C" & alt == "T")) {
        return("transition")
    } else if ((ref == "A" & (alt == "T" | alt == "C")) | (ref == "G" & (alt == "T" | alt == "C")) | (ref == "C" & (alt == "A" | alt == "G")) | (ref == "T" & (alt == "A" | alt == "G"))) {
        return("transversion")
    }
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
inputFile <- Args[6]
output <- paste0(dirname(inputFile), "/", gsub(".vcf.gz", "", basename(inputFile)), "_polysite.bed")

# 4. variable setting of test module---------------------------------------
# inputFile<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/CAF/variant/CAF_polymor_01.vcf.gz"

# 5. process --------------------------------------------------------------
GT <- fread(cmd = paste0("zcat ", inputFile, " |grep -v \"^#\"|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");printf(\"%s%s\",a[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F)
GT[, ncol(GT) := NULL]
gt <- transpose(GT[, 3:ncol(GT)])
sample_count <- nrow(gt) - 2

alt_allele_class_count <- mclapply(gt, function(x) {
    all_sample_gt <- unlist(strsplit(unlist(x[3:sample_count]), "/|\\/"))
    allele_class <- unique(all_sample_gt)
    allele_class_demissing <- allele_class[allele_class != "."]

    allele <- c(unlist(x[1]), unlist(strsplit(unlist(x[2]), ",")))
    id <- 0:(length(allele) - 1)
    allele <- allele[match(allele_class_demissing, id)]
    type <- get_types(allele = allele)
    major_type <- c("ref", sapply(type[2:length(type)], function(x) {
        if (x == "deletion" | x == "insertion") {
            return("indels")
        } else if (x == "transversion" | x == "transition") {
            return("snps")
        }
    }))
    data <- data.frame(allele, id, type, major_type)

    return()
}, mc.cores = 4)
