#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)
library(Biostrings)
# 2. functions ------------------------------------------------------------
get_derived_types_old <- function(ancestor = NULL, derived_allele = NULL) {
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
                get_type_for_one_base(ancestor = substr(ancestor, different_index, different_index), derived_allele = substr(allele.i, different_index, different_index))
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

get_type_for_one_base_old <- function(ancestor = NULL, derived_allele = NULL) {
    if (derived_allele == "*") {
        return("deletion")
    } else if (nchar(derived_allele) > nchar(ancestor)) {
        return("insertion")
    } else {
        return(paste0(ancestor, "2", derived_allele))
    }
}

get_derived_types <- function(ancestor = NULL, derived_allele = NULL) {
    derived_allele_types <- ""
    if (ancestor == "*") {
        # derived_allele_types <- rep("insertion", length(derived_allele))
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
                # return("insertion")
                align <- pairwiseAlignment(ancestor, allele.i)
                ancestor_first_gap_index <- unlist(gregexpr("-", as.character(align@pattern)))[1]
                if (ancestor_first_gap_index > 0) {
                    base <- substr(as.character(align@subject), ancestor_first_gap_index, ancestor_first_gap_index)
                    result <- paste0(base, "_insert")
                    return(result)
                } else {
                    ancestor_length <- nchar(ancestor)
                    base <- substr(allele.i, ancestor_length + 1, ancestor_length + 1)
                    result <- paste0(base, "_insert")
                    return(result)
                }
            } else if (nchar(allele.i) < nchar(ancestor)) {
                # return("deletion")
                if (allele.i == "*") {
                    base <- substr(ancestor, 1, 1)
                    result <- paste0(base, "_del")
                    return(result)
                } else {
                    align <- pairwiseAlignment(ancestor, allele.i)
                    derived_first_gap_index <- unlist(gregexpr("-", as.character(align@subject)))[1]
                    if (derived_first_gap_index > 0) {
                        base <- substr(as.character(align@pattern), derived_first_gap_index, derived_first_gap_index)
                        result <- paste0(base, "_del")
                        return(result)
                    } else {
                        subject_length <- nchar(allele.i)
                        base <- substr(ancestor, subject_length + 1, subject_length + 1)
                        result <- paste0(base, "_del")
                        return(result)
                    }
                }
            }
        })
    }
    result <- data.frame(derived_allele = derived_allele, id = names(derived_allele), variant_types = derived_allele_types)
    result$whole_types <- sapply(result$variant_types, function(x) {
        if (x == "A_insert" | x == "T_insert" | x == "C_insert" | x == "G_insert") {
            return("insertion")
        } else if (x == "A_del" | x == "T_del" | x == "C_del" | x == "G_del") {
            return("deletion")
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
        return(paste0(ancestor, "_del"))
    } else if (nchar(derived_allele) > nchar(ancestor)) {
        base <- substr(derived_allele, 2, 2)
        return(paste0(base, "_insert"))
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
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/CAF/variant/CAF_PASS_01.vcf.gz"

# 5. process --------------------------------------------------------------
vcf <- fread(cmd = paste0(bcftools, " view -H ", vcf_data, " |awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s%s%s\",b[1],OFS,b[2],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, c(6, ncol(vcf)) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", "ancestor", paste("popu", 1:(ncol(vcf) - 5), sep = ""))
chrom <- vcf$chrom[1]

vcf_ancestor_nomiss <- vcf[vcf$ancestor != ".", ] # 去除祖先型没有call出来的位点
ancestor_popu_notallsame_index <- which(apply(vcf_ancestor_nomiss[, 5:ncol(vcf_ancestor_nomiss)], 1, function(x) {
    ancestor <- x[1]
    popu <- x[-1]
    popu_nomiss <- popu[popu != "."]
    any(popu_nomiss != ancestor)
}))
data_for_mutation_load <- vcf_ancestor_nomiss[ancestor_popu_notallsame_index, ] # 去除祖先型和该群体的等位完全一致的位点
gt <- transpose(data_for_mutation_load[, 3:ncol(data_for_mutation_load)])
popu_count <- nrow(gt) - 3

mutation_load <- mclapply(gt, function(x) {
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)
    ancestor_id <- x[3]
    ancestor <- allele[ancestor_id]
    derived_allele <- allele[-which(allele == ancestor)]
    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele) ## TODO:

    allele_class <- unique(x[4:popu_count])
    allele_class_nomiss <- allele_class[allele_class != "." & allele_class != ancestor_id]
    exist_allele_data <- derived_allele_types[match(allele_class_nomiss, derived_allele_types$id), ]

    transition_load <- length(which(exist_allele_data$whole_types == "transition"))
    transversion_load <- length(which(exist_allele_data$whole_types == "transversion"))
    snp_load <- transition_load + transversion_load
    insertion_load <- length(which(exist_allele_data$whole_types == "insertion"))
    deletion_load <- length(which(exist_allele_data$whole_types == "deletion"))
    indel_load <- insertion_load + deletion_load
    ts_direc <- sapply(c("A2G", "G2A", "T2C", "C2T"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    tv_direc <- sapply(c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    insertion_direc <- sapply(c("A_insert", "T_insert", "C_insert", "G_insert"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    del_direc <- sapply(c("A_del", "T_del", "C_del", "G_del"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    # result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc)
    result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc, insertion_direc, del_direc)


    return(result)
}, mc.cores = 4)
mutation_load <- do.call(rbind, mutation_load)
result <- data.frame(chrom = chrom, bed0 = data_for_mutation_load$pos - 1, bed1 = data_for_mutation_load$pos, mutation_load = mutation_load)
colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
# write.table(format(result, scientific = FALSE), output, quote = F, sep = "\t", row.names = F, col.names = F)
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)



# mutation_load <- mclapply(gt, function(x) {
# allele <- c(x[1], unlist(strsplit(x[2], ",")))
# names(allele) <- 0:(length(allele) - 1)
# ancestor_id <- x[3]
# ancestor <- allele[ancestor_id]
# derived_allele <- allele[-which(allele == ancestor)]
#     derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele) ## TODO:

#     allele_class <- unique(x[4:popu_count])
#     allele_class_nomiss <- allele_class[allele_class != "." & allele_class != ancestor_id]
#     exist_allele_data <- derived_allele_types[match(allele_class_nomiss, derived_allele_types$id), ]

#     transition_load <- length(which(exist_allele_data$whole_types == "transition"))
#     transversion_load <- length(which(exist_allele_data$whole_types == "transversion"))
#     snp_load <- transition_load + transversion_load
#     insertion_load <- length(which(exist_allele_data$variant_types == "insertion"))
#     deletion_load <- length(which(exist_allele_data$variant_types == "deletion"))
#     indel_load <- insertion_load + deletion_load
#     ts_direc <- sapply(c("A2G", "G2A", "T2C", "C2T"), function(x) {
#         length(which(exist_allele_data$variant_types == x))
#     })
#     tv_direc <- sapply(c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G"), function(x) {
#         length(which(exist_allele_data$variant_types == x))
#     })
#     result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc)

#     return(result)
# }, mc.cores = 4)
# mutation_load <- do.call(rbind, mutation_load)
# result <- data.frame(chrom = chrom, bed0 = data_for_mutation_load$pos - 1, bed1 = data_for_mutation_load$pos, mutation_load = mutation_load)
# colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
# write.table(format(result, scientific = FALSE), output, quote = F, sep = "\t", row.names = F, col.names = F)
