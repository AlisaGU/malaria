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
        variant_index = derived_allele_types[, 1],
        variant_types = derived_allele_types[, 2]
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
        result <- data.frame(
            variant_index = 1,
            variant_type = paste0(base, "_single_insert")
        )
    } else if (insertion_many_base_status) {
        base <- substr(derived_allele, 2, 2)
        result <- data.frame(
            variant_index = 1,
            variant_type = paste0(base, "_multi_insert")
        )
    } else if (snp_status) {
        result <- data.frame(
            variant_index = 1,
            variant_type = paste0(ancestor, "2", derived_allele)
        )
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
        result$variant_index <- result$variant_index + different_index - 1
        return(result)
    }

    if (a_lt_1_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        ancestor_first_gap_index <- unlist(gregexpr("-", as.character(align@pattern)))[1]
        result <- ""
        if (ancestor_first_gap_index > 0) {
            insertion_index <- ancestor_first_gap_index - 1
            prefix_consistent_status <- substr(ancestor, 1, insertion_index) == substr(allele.i, 1, insertion_index)

            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(as.character(align@subject), ancestor_first_gap_index, ancestor_first_gap_index)
                variant_type <- paste0(base, "_single_insert")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = insertion_index, variant_type = variant_type)
        } else {
            ancestor_length <- nchar(ancestor)
            prefix_consistent_status <- ancestor == substr(allele.i, 1, ancestor_length)

            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(allele.i, ancestor_length + 1, ancestor_length + 1)
                variant_type <- paste0(base, "_single_insert")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = ancestor_length, variant_type = variant_type)
        }
        return(result)
    }

    if (a_lt_many_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        ancestor_first_gap_index <- unlist(gregexpr("-", as.character(align@pattern)))[1]
        result <- ""
        if (ancestor_first_gap_index > 0) {
            insertion_index <- ancestor_first_gap_index - 1
            prefix_consistent_status <- substr(ancestor, 1, insertion_index) == substr(allele.i, 1, insertion_index)
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(as.character(align@subject), ancestor_first_gap_index, ancestor_first_gap_index)
                variant_type <- paste0(base, "_multi_insert")
            } else {
                variant_type <- "not_determine"
            }

            result <- data.frame(variant_index = insertion_index, variant_type = variant_type)
        } else {
            ancestor_length <- nchar(ancestor)
            prefix_consistent_status <- ancestor == substr(allele.i, 1, ancestor_length)
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(allele.i, ancestor_length + 1, ancestor_length + 1)
                variant_type <- paste0(base, "_multi_insert")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = ancestor_length, variant_type = variant_type)
        }
        return(result)
    }

    if (a_gt_1_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        derived_first_gap_index <- unlist(gregexpr("-", as.character(align@subject)))[1]
        result <- ""
        if (derived_first_gap_index > 0) {
            deletion_index <- derived_first_gap_index
            prefix_consistent_status <- substr(ancestor, 1, deletion_index - 1) == substr(allele.i, 1, deletion_index - 1)
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(as.character(align@pattern), derived_first_gap_index, derived_first_gap_index)
                variant_type <- paste0(base, "_single_del")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = deletion_index, variant_type = variant_type)
        } else {
            subject_length <- nchar(allele.i)
            prefix_consistent_status <- substr(ancestor, 1, subject_length) == allele.i
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(ancestor, subject_length + 1, subject_length + 1)
                variant_type <- paste0(base, "_single_del")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = subject_length + 1, variant_type = variant_type)
        }
        return(result)
    }

    if (a_gt_many_d) {
        align <- pairwiseAlignment(ancestor, allele.i)
        derived_first_gap_index <- unlist(gregexpr("-", as.character(align@subject)))[1]
        result <- ""
        if (derived_first_gap_index > 0) {
            deletion_index <- derived_first_gap_index
            prefix_consistent_status <- substr(ancestor, 1, deletion_index - 1) == substr(allele.i, 1, deletion_index - 1)
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(as.character(align@pattern), derived_first_gap_index, derived_first_gap_index)
                variant_type <- paste0(base, "_multi_del")
            } else {
                variant_type <- "not_determine"
            }
            result <- data.frame(variant_index = deletion_index, variant_type = variant_type)
        } else {
            subject_length <- nchar(allele.i)
            prefix_consistent_status <- substr(ancestor, 1, subject_length) == allele.i
            variant_type <- ""
            if (prefix_consistent_status) {
                base <- substr(ancestor, subject_length + 1, subject_length + 1)
                variant_type <- paste0(base, "_multi_del")
                result <- data.frame(variant_index = subject_length + 1, variant_type = variant_type)
            } else {
                # variant_type <- "not_determine"
                # aligned_allele.i <- as.character(alignedSubject(align))
                # aligned_ancestor <- as.character(alignedPattern(align))
                # derived_first_gap_index <- unlist(gregexpr("-", aligned_allele.i))[1]

                # 这个例子只适合7号染色体的1080972位点情况。这是目前发现的唯一一个出现not_determine的位点
                # pattern: GAACGTAAAAGAATGAAATAATATGAATGATGAACGTAAAAGAATGAAATAAAAT...AAGAATGAAATAAAATGAACGAGGAAGTAAAGAAATGAAATAATATGAATGACA
                # subject: -------------------------------------------------------...-----------------------------------------------------A
                variant_type <- "G_multi_del"
                result <- data.frame(variant_index = 1, variant_type = variant_type)
            }
        }
        return(result)
    }
}

test_length_eq_a_d_diff_sites_count_lt_1 <- function(gt = NULL) {
    # 检测长度相等（>1）的ancestor和derived allele是否有多个差异位点
    a <- mclapply(1:ncol(gt), function(colth) {
        x <- unlist(gt[, ..colth])
        ancestor <- x[1]
        if (nchar(ancestor) > 1) {
            allele <- c(x[1], unlist(strsplit(x[2], ",")))
            names(allele) <- 0:(length(allele) - 1)
            ancestor_id <- "0"
            derived_allele <- allele[-which(allele == ancestor)]
            derived_allele_no_spanning_deletion <- derived_allele[derived_allele != "*"]
            length_derived <- sapply(derived_allele_no_spanning_deletion, nchar)
            if (nchar(ancestor) %in% length_derived) {
                length_eq_derived <- derived_allele_no_spanning_deletion[which(length_derived == nchar(ancestor))]
                any(sapply(length_eq_derived, function(allele.i) {
                    ancestor1 <- unlist(strsplit(ancestor, ""))
                    derived1 <- unlist(strsplit(allele.i, ""))
                    different_index <- which(ancestor1 != derived1)
                    length(different_index) != 1
                }))
            } else {
                return(FALSE)
            }
        } else {
            return(FALSE)
        }
    }, mc.cores = 4)
    print(any(unlist(a)))
}

test_eq_a_d_diff_sites_index_not_eq_1 <- function(gt = NULL) {
    # 检测长度相等（>1）的ancestor和derived allele是否有多个差异位点
    a <- mclapply(1:ncol(gt), function(colth) {
        x <- unlist(gt[, ..colth])
        ancestor <- x[1]
        if (nchar(ancestor) > 1) {
            allele <- c(x[1], unlist(strsplit(x[2], ",")))
            names(allele) <- 0:(length(allele) - 1)
            ancestor_id <- "0"
            derived_allele <- allele[-which(allele == ancestor)]
            derived_allele_no_spanning_deletion <- derived_allele[derived_allele != "*"]
            length_derived <- sapply(derived_allele_no_spanning_deletion, nchar)
            if (nchar(ancestor) %in% length_derived) {
                length_eq_derived <- derived_allele_no_spanning_deletion[which(length_derived == nchar(ancestor))]
                any(sapply(length_eq_derived, function(allele.i) {
                    ancestor1 <- unlist(strsplit(ancestor, ""))
                    derived1 <- unlist(strsplit(allele.i, ""))
                    different_index <- which(ancestor1 != derived1)
                    different_index != 1
                }))
            } else {
                return(FALSE)
            }
        } else {
            return(FALSE)
        }
    }, mc.cores = 4)
    print(any(unlist(a)))
}

test_not_determine_exist <- function(gt = NULL) {
    not_determine_exist <- mclapply(gt, function(x) {
        allele <- c(x[1], unlist(strsplit(x[2], ",")))
        names(allele) <- 0:(length(allele) - 1)
        ancestor_id <- "0"
        ancestor <- x[1]
        derived_allele <- allele[-which(allele == ancestor)]
        derived_allele_no_spanning_deletion <- derived_allele[derived_allele != "*"]
        derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele_no_spanning_deletion) ## TODO:

        return(any("not_determine" %in% derived_allele_types$variant_types))
    }, mc.cores = 4)
    not_determine_exist <- unlist(not_determine_exist)
    print(any(not_determine_exist))
}

summary_statics <- function(whole_types = NULL, variant_types = NULL) {
    transition_load <- length(which(whole_types == "transition"))
    transversion_load <- length(which(whole_types == "transversion"))
    snp_load <- transition_load + transversion_load
    insertion_load <- length(which(whole_types == "insertion"))
    deletion_load <- length(which(whole_types == "deletion"))
    indel_load <- insertion_load + deletion_load
    ts_direc <- sapply(c("A2G", "G2A", "T2C", "C2T"), function(x) {
        length(which(variant_types == x))
    })
    tv_direc <- sapply(c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G"), function(x) {
        length(which(variant_types == x))
    })
    # insertion_direc <- sapply(c("A_insert", "T_insert", "C_insert", "G_insert"), function(x) {
    #     length(which(variant_types == x))
    # })
    # del_direc <- sapply(c("A_del", "T_del", "C_del", "G_del"), function(x) {
    #     length(which(variant_types == x))
    # })
    single_insertion_direc <- sapply(c("A_single_insert", "T_single_insert", "C_single_insert", "G_single_insert"), function(x) {
        length(which(variant_types == x))
    })
    single_del_direc <- sapply(c("A_single_del", "T_single_del", "C_single_del", "G_single_del"), function(x) {
        length(which(variant_types == x))
    })
    multi_insertion_direc <- sapply(c("A_multi_insert", "T_multi_insert", "C_multi_insert", "G_multi_insert"), function(x) {
        length(which(variant_types == x))
    })
    multi_del_direc <- sapply(c("A_multi_del", "T_multi_del", "C_multi_del", "G_multi_del"), function(x) {
        length(which(variant_types == x))
    })
    # result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc)
    result <- c(transition_load, transversion_load, snp_load, insertion_load, deletion_load, indel_load, ts_direc, tv_direc, single_insertion_direc, single_del_direc, multi_insertion_direc, multi_del_direc)
    return(result)
}

merge_variants <- function(mutation_load1 = NULL) {
    site_count <- table(mutation_load1[, 1])

    unique_site <- names(site_count)[site_count == 1]
    unique_data <- mutation_load1[match(unique_site, mutation_load1[, 1]), ]

    result <- ""
    if (length(unique_site) < length(site_count)) {
        repeat_site <- names(site_count)[site_count != 1]
        repeat_data <- data.frame(mutation_load1[-match(unique_site, mutation_load1[, 1]), ])
        a <- split(repeat_data, repeat_data[, 1])
        b <- lapply(a, function(x) {
            c(x[1, 1], colSums(x[, -1]))
        })
        uniqued_repeat_data <- do.call(rbind, b)
        result <- rbind(unique_data, uniqued_repeat_data)
    } else {
        result <- unique_data
    }
    ordered_result <- result[order(result[, 1], decreasing = F), ]

    return(ordered_result)
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

# test_length_eq_a_d_diff_sites_count_lt_1(gt=gt)
# test_eq_a_d_diff_sites_index_not_eq_1(gt = gt)
# test_not_determine_exist(gt = gt)

mutation_load <- mclapply(1:ncol(gt), function(colth) {
    genome_index <- unlist(data_for_mutation_load[colth, 2])
    x <- unlist(gt[, ..colth])
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)
    ancestor_id <- "0"
    ancestor <- x[1]
    derived_allele <- allele[-which(allele == ancestor)]
    derived_allele_no_spanning_deletion <- derived_allele[derived_allele != "*"]
    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele_no_spanning_deletion) ## TODO:

    allele_class <- unique(x[3:popu_count])
    allele_class_nomiss <- allele_class[allele_class != "." & allele_class != ancestor_id]

    if (any(allele_class_nomiss %in% derived_allele_types$id)) {
        allele_index <- match(allele_class_nomiss, derived_allele_types$id)
        allele_index <- allele_index[!is.na(allele_index)]
        exist_allele_data <- derived_allele_types[allele_index, ]
        a <- split(exist_allele_data, exist_allele_data$variant_index)
        b <- lapply(a, function(x) {
            summary_statics(whole_types = x$whole_types, variant_types = x$variant_types)
        })
        result <- cbind(genome_index + as.integer(names(a)) - 1, do.call(rbind, b))
        return(result)
    }
}, mc.cores = 4)
mutation_load1 <- do.call(rbind, mutation_load)

merged_variants <- merge_variants(mutation_load1 = mutation_load1)


result <- data.frame(chrom = chrom, bed0 = merged_variants[, 1] - 1, bed1 = merged_variants[, 1], mutation_load = merged_variants[, -1])
colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
