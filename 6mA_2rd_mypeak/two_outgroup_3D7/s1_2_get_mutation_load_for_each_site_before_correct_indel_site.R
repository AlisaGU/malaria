#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)
library(Biostrings)
library(dplyr)
# 2. functions ------------------------------------------------------------
get_derived_types <- function(ancestor = NULL, derived_allele = NULL) {
    insertion <- c(paste(c("A", "T", "C", "G"), "_single_insert", sep = ""), paste(c("A", "T", "C", "G"), "_multi_insert", sep = ""))
    deletion <- c(paste(c("A", "T", "C", "G"), "_single_del", sep = ""), paste(c("A", "T", "C", "G"), "_multi_del", sep = ""))
    ts <- c("A2G", "G2A", "T2C", "C2T")
    tv <- c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")
    derived_allele_types <- ""

    if (nchar(ancestor) == 1) {
        derived_allele_types <- lapply(derived_allele, function(allele.i) {
            get_type_for_one_base(ancestor = ancestor, derived_allele = allele.i)
        })
    }

    if (nchar(ancestor) > 1) {
        derived_allele_types <- lapply(derived_allele, function(allele.i) {
            get_type_for_ancestor_gt_1(ancestor = ancestor, allele.i = allele.i)
        })
    }

    derived_allele_types <- do.call(rbind, derived_allele_types)
    result <- data.frame(
        derived_allele = derived_allele, id = names(derived_allele),
        variant_types = derived_allele_types
    )
    result$whole_types <- sapply(result$variant_types, function(x) {
        if (x %in% insertion) {
            return("insertion")
        } else if (x %in% deletion) {
            return("deletion")
        } else if (x %in% ts) {
            return("transition")
        } else if (x %in% tv) {
            return("transversion")
        } else {
            return("not_determine")
        }
    })
    return(result)
}

get_type_for_one_base <- function(ancestor = NULL, derived_allele = NULL) {
    # 适用于祖先型只有一个碱基的情况
    insertion_one_base_status <- nchar(derived_allele) == nchar(ancestor) + 1
    insertion_many_base_status <- nchar(derived_allele) > nchar(ancestor) + 1
    snp_status <- nchar(derived_allele) == nchar(ancestor)
    result <- c()

    if (insertion_one_base_status) {
        base <- substr(derived_allele, 2, 2)
        result <- paste0(base, "_single_insert")
    } else if (insertion_many_base_status) {
        base <- substr(derived_allele, 2, 2)
        result <- paste0(base, "_multi_insert")
    } else if (snp_status) {
        result <- paste0(ancestor, "2", derived_allele)
    }

    return(result)
}

get_type_for_ancestor_gt_1 <- function(ancestor = NULL, allele.i = NULL) {
    a_eq_d <- nchar(allele.i) == nchar(ancestor)
    a_lt_1_d <- nchar(allele.i) == nchar(ancestor) + 1
    a_lt_many_d <- nchar(allele.i) > nchar(ancestor) + 1
    a_gt_1_d <- nchar(allele.i) == nchar(ancestor) - 1
    a_gt_many_d <- nchar(allele.i) <= nchar(ancestor) - 1

    result <- ""
    if (a_eq_d) {
        ancestor1 <- unlist(strsplit(ancestor, ""))
        derived1 <- unlist(strsplit(allele.i, ""))
        different_index <- which(ancestor1 != derived1)
        result <- get_type_for_one_base(
            ancestor = substr(ancestor, different_index, different_index),
            derived_allele = substr(allele.i, different_index, different_index)
        )

        return(result)
    }

    if (a_lt_1_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        ancestor_first_gap_index <- unlist(gregexpr("-", as.character(align@pattern)))[1]
        result <- ""
        if (ancestor_first_gap_index > 0) {
            base <- substr(as.character(align@subject), ancestor_first_gap_index, ancestor_first_gap_index)
            result <- paste0(base, "_single_insert")
        } else {
            ancestor_length <- nchar(ancestor)
            base <- substr(allele.i, ancestor_length + 1, ancestor_length + 1)
            result <- paste0(base, "_single_insert")
        }
        return(result)
    }

    if (a_lt_many_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        ancestor_first_gap_index <- unlist(gregexpr("-", as.character(align@pattern)))[1]
        result <- ""
        if (ancestor_first_gap_index > 0) {
            base <- substr(as.character(align@subject), ancestor_first_gap_index, ancestor_first_gap_index)
            result <- paste0(base, "_multi_insert")
        } else {
            ancestor_length <- nchar(ancestor)
            base <- substr(allele.i, ancestor_length + 1, ancestor_length + 1)
            result <- paste0(base, "_multi_insert")
        }
        return(result)
    }

    if (a_gt_1_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        derived_first_gap_index <- unlist(gregexpr("-", as.character(align@subject)))[1]
        result <- ""
        if (derived_first_gap_index > 0) {
            deletion_index <- derived_first_gap_index
            base <- substr(as.character(align@pattern), derived_first_gap_index, derived_first_gap_index)
            result <- paste0(base, "_single_del")
        } else {
            subject_length <- nchar(allele.i)
            base <- substr(ancestor, subject_length + 1, subject_length + 1)
            result <- paste0(base, "_single_del")
        }
        return(result)
    }

    if (a_gt_many_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        derived_first_gap_index <- unlist(gregexpr("-", as.character(align@subject)))[1]
        result <- ""
        if (derived_first_gap_index > 0) {
            base <- substr(as.character(align@pattern), derived_first_gap_index, derived_first_gap_index)
            result <- paste0(base, "_multi_del")
        } else {
            subject_length <- nchar(allele.i)
            base <- substr(ancestor, subject_length + 1, subject_length + 1)
            result <- paste0(base, "_multi_del")
        }
        return(result)
    }
}

# 3. input ----------------------------------------------------------------
Args <- commandArgs()
vcf_data <- Args[6]
output <- paste0(dirname(vcf_data), "/", gsub(".vcf.gz", "", basename(vcf_data)), "_polysite_before_correct_indel_site.bed")

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

# test_length_eq_a_d_diff_sites_count_lt_1(gt=gt)
# test_eq_a_d_diff_sites_index_not_eq_1(gt = gt)
# test_not_determine_exist(gt = gt)

mutation_load <- mclapply(gt, function(x) {
    # genome_index <- unlist(data_for_mutation_load[colth, 2])
    # x <- unlist(gt[, ..colth])
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)
    ancestor_id <- "0"
    ancestor <- x[1]
    derived_allele <- allele[-which(allele == ancestor)]
    derived_allele_no_spanning_deletion <- derived_allele[derived_allele != "*"]
    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele_no_spanning_deletion) ## TODO:

    allele_class <- unique(x[3:popu_count])
    allele_class_nomiss <- allele_class[allele_class != "." & allele_class != ancestor_id]

    allele_index <- match(allele_class_nomiss, derived_allele_types$id)
    allele_index <- allele_index[!is.na(allele_index)]
    exist_allele_data <- derived_allele_types[allele_index, ]

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
    single_insertion_direc <- sapply(c("A_single_insert", "T_single_insert", "C_single_insert", "G_single_insert"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    single_del_direc <- sapply(c("A_single_del", "T_single_del", "C_single_del", "G_single_del"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    multi_insertion_direc <- sapply(c("A_multi_insert", "T_multi_insert", "C_multi_insert", "G_multi_insert"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    multi_del_direc <- sapply(c("A_multi_del", "T_multi_del", "C_multi_del", "G_multi_del"), function(x) {
        length(which(exist_allele_data$variant_types == x))
    })
    result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc, single_insertion_direc, single_del_direc, multi_insertion_direc, multi_del_direc)
    return(result)
}, mc.cores = 4)
mutation_load <- do.call(rbind, mutation_load)
result <- data.frame(chrom = chrom, bed0 = data_for_mutation_load$pos - 1, bed1 = data_for_mutation_load$pos, mutation_load = mutation_load)
colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
