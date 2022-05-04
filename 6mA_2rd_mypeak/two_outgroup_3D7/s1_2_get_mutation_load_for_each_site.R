#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)

# 2. functions ------------------------------------------------------------
get_derived_types <- function(ancestor = NULL, derived_allele = NULL) {
    derived_allele_types <- ""
    if (ancestor == "*") {
        derived_allele_types <- rep("insertion", length(derived_allele))
    } else if (ancestor != "*" & nchar(ancestor) == 1) {
        derived_allele_types <- sapply(derived_allele, function(allele.i) {
            get_type_for_one_base(ancestor = ancestor, derived_allele = allele.i)
        })
    } else if (nchar(ancestor) > 1) {
        derived_allele_types <- sapply(derived_allele, function(allele.i) {
            if (nchar(allele.i) == nchar(ancestor)) {
                ancestor1 <- unlist(strsplit(ancestor, ""))
                derived1 <- unlist(strsplit(allele.i, ""))
                different_index <- which(ancestor1 != derived1)
                get_type_for_one_base(
                    ancestor = substr(ancestor, different_index, different_index),
                    derived_allele = substr(allele.i, different_index, different_index)
                )
            } else if (nchar(allele.i) > nchar(ancestor)) {
                return("insertion")
            } else if (nchar(allele.i) < nchar(ancestor)) {
                return("deletion")
            }
        })
    }
    result <- data.frame(derived_allele = derived_allele, id = names(derived_allele), variant_types = derived_allele_types)
    result$whole_types <- sapply(result$variant_types, function(x) {
        if (x == "deletion" | x == "insertion") {
            return("indels")
        } else if (x %in% c("A2G", "G2A", "T2C", "C2T")) {
            return("transition")
        } else if (x %in% c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")) {
            return("transversion")
        }
    })
    return(result)
}

get_type_for_one_base <- function(ancestor = NULL, derived_allele = NULL) {
    if (derived_allele == "*") {
        return("deletion")
    } else if (nchar(derived_allele) > nchar(ancestor)) {
        return("insertion")
    } else {
        return(paste0(ancestor, "2", derived_allele))
    }
}
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
vcf_data <- Args[6]
output <- paste0(dirname(vcf_data), "/", gsub(".vcf.gz", "", basename(vcf_data)), "_polysite.bed")

bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/OCE/variant/OCE_PASS_01.vcf.gz"

# 5. process --------------------------------------------------------------
vcf <- fread(cmd = paste0(bcftools, " view -H ", vcf_data, " |awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s%s%s\",b[1],OFS,b[2],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, c(ncol(vcf) - 2, ncol(vcf)) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", paste("popu", 1:(ncol(vcf) - 5), sep = ""), "outgroup1")
chrom <- vcf$chrom[1]

# vcf_ancestor_nomiss <- vcf[vcf$outgroup1 != ".", ] # 去除祖先型没有call出来的位点
# ancestor_popu_notallsame_index <- which(apply(vcf_ancestor_nomiss[, 5:ncol(vcf_ancestor_nomiss)], 1, function(x) {
#     ancestor <- x[1]
#     popu <- x[-1]
#     popu_nomiss <- popu[popu != "."]
#     any(popu_nomiss != ancestor)
# }))
# data_for_mutation_load <- vcf_ancestor_nomiss[ancestor_popu_notallsame_index, ] # 去除祖先型和该群体的等位完全一致的位点

data_for_mutation_load <- vcf[vcf$outgroup1 == "0", ] # 选择ref和outgroup1碱基状态一致的位点
gt <- transpose(data_for_mutation_load[, 3:(ncol(data_for_mutation_load) - 1)])
popu_count <- nrow(gt) - 2
mutation_load <- mclapply(gt, function(x) {
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)
    ancestor_id <- "0"
    ancestor <- x[1]
    derived_allele <- allele[-which(allele == ancestor)]
    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele) ## TODO:

    allele_class <- unique(x[3:popu_count])
    allele_class_nomiss <- allele_class[allele_class != "." & allele_class != ancestor_id]
    exist_allele_data <- derived_allele_types[match(allele_class_nomiss, derived_allele_types$id), ]

    transition_load <- length(which(exist_allele_data$whole_types == "transition"))
    transversion_load <- length(which(exist_allele_data$whole_types == "transversion"))
    snp_load <- transition_load + transversion_load
    insertion_load <- length(which(exist_allele_data$variant_types == "insertion"))
    deletion_load <- length(which(exist_allele_data$variant_types == "deletion"))
    indel_load <- insertion_load + deletion_load
    ts_direc <- sapply(c("A2G", "G2A", "T2C", "C2T"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    tv_direc <- sapply(c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc)

    return(result)
}, mc.cores = 4)
mutation_load <- do.call(rbind, mutation_load)
result <- data.frame(chrom = chrom, bed0 = data_for_mutation_load$pos - 1, bed1 = data_for_mutation_load$pos, mutation_load = mutation_load)
colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
write.table(format(result, scientific = FALSE), output, quote = F, sep = "\t", row.names = F, col.names = F)
