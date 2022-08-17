#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table, warn.conflicts = FALSE, quietly = TRUE)
library(parallel)
library(Biostrings)
library(dplyr)
library(seqinr)
# 2. functions ------------------------------------------------------------
get_derived_types <- function(ancestor = NULL, derived_allele = NULL) {
    deletion <- c(paste(c("A", "T", "C", "G"), "_single_del", sep = ""))

    ts <- c("A2G", "G2A", "T2C", "C2T")
    tv <- c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")
    derived_allele_types <- sapply(derived_allele, function(x) {
        get_type_for_one_base(ancestor = ancestor, derived_allele = x)
    })

    result <- data.frame(
        derived_allele = derived_allele, id = names(derived_allele),
        variant_types = derived_allele_types
    )
    result$whole_types <- sapply(result$variant_types, function(x) {
        if (x %in% deletion) {
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
    result <- c()

    if (derived_allele == "*") {
        result <- paste0(ancestor, "_single_del")
    } else {
        result <- paste0(ancestor, "2", derived_allele)
    }

    return(result)
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



delete_continuous_deletion <- function(data = NULL) {
    deletion_index <- which(data$alt == "*")
    rundiff <- c(1, diff(deletion_index))
    difflist <- split(deletion_index, cumsum(rundiff != 1))
    deletion_region_len <- sapply(difflist, length)
    continuous_region_index <- which(deletion_region_len != 1)
    continuous_region_vcf_index <- as.integer(unlist(difflist[continuous_region_index]))
    result <- data[-continuous_region_vcf_index, ]
    return(result)
}

get_deletion_result_from_vcf <- function(vcf = NULL) {
    a <- vcf %>%
        filter(P_falciparum == 0 & grepl("*", alt, fixed = T))
    b <- t(apply(a, 1, function(x) {
        allele <- c(x[3], unlist(strsplit(x[4], ",")))
        names(allele) <- 0:(length(allele) - 1)
        x[5:7] <- allele[x[5:7]]
        x[5:7] <- gsub("A|T|C|G", "N", x[5:7])
        return(x)
    }))
    b <- as.data.table(b)
    vcf_deletion <- split(b, by = c("P_falciparum", "outgroupA", "B_species"), drop = T)
    result <- lapply(vcf_deletion, function(x) {
        result_deletion <- get_deletion_type(Pos = as.numeric(x$pos))
    })
    result <- do.call(rbind, result)
    return(result)
}

get_deletion_result_from_vcf_more_loose <- function(based_vcf = NULL) {
    outgroupA_deletion_index <- get_deletion_index_based_on_ref(based_vcf = based_vcf, species_type = "outgroupA")
    speciesB_deletion_index <- get_deletion_index_based_on_ref(based_vcf = based_vcf, species_type = "speciesB")
    start_pos <- sort(unique(c(outgroupA_deletion_index[, 1], speciesB_deletion_index[, 1])))
    deletion_type <- "single_del"
    # base_type是错的，但因为之后不会用到，只是为了满足格式需要，加上了
    base_type <- toupper(sapply(start_pos, function(x) {
        substr(chrom_sequence, x, x)
    }))

    result <- data.frame(chrom = chrom, bed0 = start_pos - 1, bed1 = start_pos, variant_types = paste(base_type, deletion_type, sep = "_"), whole_types = "deletion")

    a <- t(sapply(1:nrow(result), function(x) {
        summary_statics(whole_types = result$whole_types[x], variant_types = result$variant_types[x])
    }))
    final_result <- data.frame(result[, c("chrom", "bed0", "bed1")], a[, 1:6])
    colnames(final_result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
    return(final_result)
}

get_deletion_index_based_on_ref <- function(based_vcf = NULL, species_type = NULL) {
    ref_nogap_index <- which(based_vcf[, 5] != "*")
    ref_nogap_seq <- based_vcf[ref_nogap_index, 5]
    names(ref_nogap_seq) <- based_vcf[ref_nogap_index, 2]
    species_index <- ifelse(species_type == "outgroupA", 6, 7)

    outgroupA_seq <- sapply(1:length(ref_nogap_index), function(x) {
        start_index <- ref_nogap_index[x]

        if (x == length(ref_nogap_index)) {
            end_index <- start_index
        } else {
            end_index <- ref_nogap_index[x + 1] - 1
        }
        seq1 <- based_vcf[start_index:end_index, species_index]
        if (length(seq1) > 1 & any(seq1 != "*")) {
            seq1 <- seq1[seq1 != "*"]
        } else if (length(seq1) > 1 & all(seq1 == "*")) {
            seq1 <- "*"
        }
        result <- paste(seq1, collapse = "")
        return(result)
    })

    outgroupA_deletion_index <- get_deletion_index_region(Pos = ref_nogap_index[outgroupA_seq == "*"])
    collapse_index <- c()

    if (nrow(outgroupA_deletion_index) > 1) {
        for (i in 1:c(nrow(outgroupA_deletion_index) - 1)) {
            first_end <- outgroupA_deletion_index[i, 2]
            second_end <- outgroupA_deletion_index[i + 1, 1]
            ref_between_seq <- based_vcf[(first_end + 1):(second_end - 1), 5]
            if (all(ref_between_seq == "*")) {
                collapse_index <- c(collapse_index, i)
            }
        }
    }

    if (!is.null(collapse_index)) {
        collapse_index <- t(sapply(collapse_index, function(x) {
            c(x, x + 1)
        }))
        norelated_data <- outgroupA_deletion_index[-as.integer(collapse_index), ]
        related_data <- t(apply(collapse_index, 1, function(x) {
            c(outgroupA_deletion_index[x[1], 1], outgroupA_deletion_index[x[2], 2])
        }))
        collapsed_data <- rbind(norelated_data, related_data)
        collapsed_data <- collapsed_data[order(collapsed_data[, 1]), ]
    } else {
        collapsed_data <- outgroupA_deletion_index
    }
    return(collapsed_data)
}


get_insertion_result_from_vcf_more_loose <- function(based_vcf = NULL) {


}

get_insertion_index_based_on_ref <- function(based_vcf = NULL, species_type = NULL) {
  ref_nogap_index <- which(based_vcf[, 5] != "*")
  ref_nogap_seq <- based_vcf[ref_nogap_index, 5]
  names(ref_nogap_seq) <- based_vcf[ref_nogap_index, 2]

   outgroupA_seq <- get_species_seq_based_on_ref(species_type = "outgroupA", based_vcf = based_vcf)
   speciesB_seq <- get_species_seq_based_on_ref(species_type = "speciesB", based_vcf = based_vcf)
    insertion_index <- sort(union(which(sapply(outgroupA_seq, nchar) > 1), which(sapply(speciesB_seq, nchar) > 1)))
    insertion_pos <- names(ref_nogap_seq)[insertion_index]
    return(insertion_pos)
}

get_species_seq_based_on_ref<-function(species_type=NULL,based_vcf = NULL){
 ref_nogap_index <- which(based_vcf[, 5] != "*")
    ref_nogap_seq <- based_vcf[ref_nogap_index, 5]
    names(ref_nogap_seq) <- based_vcf[ref_nogap_index, 2]
    species_index <- ifelse(species_type == "outgroupA", 6, 7)

    species_seq <- sapply(1:length(ref_nogap_index), function(x) {
        start_index <- ref_nogap_index[x]

        if (x == length(ref_nogap_index)) {
            end_index <- start_index
        } else {
            end_index <- ref_nogap_index[x + 1] - 1
        }
        seq1 <- based_vcf[start_index:end_index, species_index]
        if (length(seq1) > 1 & any(seq1 != "*")) {
            seq1 <- seq1[seq1 != "*"]
        } else if (length(seq1) > 1 & all(seq1 == "*")) {
            seq1 <- "*"
        }
        result <- paste(seq1, collapse = "")
        return(result)
    })
    return(species_seq)
}

get_deletion_index_region <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    region <- t(sapply(difflist, function(x) {
        c(x[1], x[length(x)])
    }))
    return(region)
}

get_deletion_index <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    start_pos <- sapply(difflist, function(x) {
        x[1]
    })
    return(start_pos)
}

get_insertion_index <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    start_pos <- sapply(difflist, function(x) {
        x[1] - 1
    })
    return(start_pos)
}

get_deletion_type <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    # start_pos <- sapply(difflist, function(x) {
    #     x[1] - 1
    # })
    start_pos <- sapply(difflist, function(x) {
        x[1]
    })
    deletion_type <- sapply(difflist, function(x) {
        if (length(x) == 1) {
            return("single_del")
        } else {
            return("multi_del")
        }
    })

    # base_type是错的，但因为之后不会用到，只是为了满足格式需要，加上了
    base_type <- toupper(sapply(start_pos, function(x) {
        substr(chrom_sequence, x, x)
    }))

    result <- data.frame(chrom = chrom, bed0 = start_pos - 1, bed1 = start_pos, variant_types = paste(base_type, deletion_type, sep = "_"), whole_types = "deletion")

    a <- t(sapply(1:nrow(result), function(x) {
        summary_statics(whole_types = result$whole_types[x], variant_types = result$variant_types[x])
    }))
    final_result <- data.frame(result[, c("chrom", "bed0", "bed1")], a[, 1:6])
    colnames(final_result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
    return(final_result)
}

get_insertion_result_from_vcf <- function(vcf = NULL) {
    a <- vcf %>%
        filter(P_falciparum != 0 & grepl("*", alt, fixed = T) & !(P_falciparum == outgroupA & P_falciparum == B_species))
    b <- t(apply(a, 1, function(x) {
        allele <- c(x[3], unlist(strsplit(x[4], ",")))
        names(allele) <- 0:(length(allele) - 1)
        x[5:7] <- allele[x[5:7]]
        x[5:7] <- gsub("A|T|C|G", "N", x[5:7])
        return(x)
    }))
    b <- as.data.table(b)
    vcf_insertion <- split(b, by = c("P_falciparum", "outgroupA", "B_species"), drop = T)
    result <- lapply(vcf_insertion, function(x) {
        result_insertion <- get_insertion_type(Pos = as.numeric(x$pos))
    })
    result <- do.call(rbind, result)
    return(result)
}

get_insertion_type <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    start_pos <- sapply(difflist, function(x) {
        x[1] - 1
    })
    insertion_type <- sapply(difflist, function(x) {
        if (length(x) == 1) {
            return("single_insert")
        } else {
            return("multi_insert")
        }
    })

    # base_type是错的，但因为之后不会用到，只是为了满足格式需要，加上了
    base_type <- toupper(sapply(start_pos, function(x) {
        substr(chrom_sequence, x, x)
    }))

    result <- data.frame(chrom = chrom, bed0 = start_pos - 1, bed1 = start_pos, variant_types = paste(base_type, insertion_type, sep = "_"), whole_types = "insertion")

    a <- t(sapply(1:nrow(result), function(x) {
        summary_statics(whole_types = result$whole_types[x], variant_types = result$variant_types[x])
    }))
    final_result <- data.frame(result[, c("chrom", "bed0", "bed1")], a[, 1:6])
    colnames(final_result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
    return(final_result)
}

merge_variants <- function(data = NULL) {
    a <- table(data$bed0)
    unique_index <- a[a == 1]
    unique_pos <- names(unique_index)
    unique_data <- data[match(unique_pos, data$bed0), ]
    duplicated_data <- data[-match(unique_pos, data$bed0), ]
    data_split <- split(duplicated_data, duplicated_data$bed0, drop = T)
    deduplicated <- as.data.frame(t(sapply(data_split, function(x) {
        result <- c(chrom = unlist(x[1, 1]), bed0 = unlist(x[1, 2]), bed1 = unlist(x[1, 3]), colSums(x[, -c(1:3)]))
    })))
    deduplicated[, -1] <- sapply(deduplicated[, -1], as.integer)
    result <- rbind(unique_data, deduplicated)
    result <- result[order(result$bed0), ]
    return(result)
}
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
vcf_data <- Args[6]
output <- paste0(dirname(vcf_data), "/", gsub(".vcf", "", basename(vcf_data)), "_polysite.bed")
genome <- read.fasta("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta", seqtype = "DNA", as.string = T)
chrom <- unlist(strsplit(basename(vcf_data), ".", fixed = T))[[1]]
bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/Pf3D7_01_v3.species3.vcf.gz"

# 5. process --------------------------------------------------------------
chrom_sequence <- genome[[chrom]]
vcf <- fread(cmd = paste0("cat ", vcf_data, "|grep -v \"^#\"|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s\",b[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, c(8) := NULL]
# vcf[, ncol(vcf) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", "P_falciparum", "outgroupA", "B_species")

vcf_snp <- vcf %>%
    filter(P_falciparum == 0) %>%
    filter(grepl("A|T|C|G", alt))
vcf_snp$alt <- gsub(",*", "", gsub("*,", "", vcf_snp$alt, fixed = T), fixed = T)
gt_snp <- transpose(vcf_snp[, 3:4])
mutation_load_snp <- mclapply(1:ncol(gt_snp), function(colth) {
    x <- unlist(gt_snp[, ..colth])
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)

    ancestor_id <- "0"
    ancestor <- x[1]
    derived_allele <- allele[-which(allele == ancestor)]

    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele) ## TODO:
    result <- summary_statics(whole_types = derived_allele_types$whole_types, variant_types = derived_allele_types$variant_types)
    return(result)
}, mc.cores = 4)

mutation_load_snp1 <- do.call(rbind, mutation_load_snp)
result_snp <- data.frame(chrom = chrom, bed0 = vcf_snp[, 2] - 1, bed1 = vcf_snp[, 2], mutation_load = mutation_load_snp1[, 1:6])
colnames(result_snp)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")


based_vcf <- t(apply(vcf, 1, function(x) {
    allele <- c(x[3], unlist(strsplit(x[4], ",")))
    names(allele) <- 0:(length(allele) - 1)
    x[5:7] <- allele[x[5:7]]
    # x[5:7] <- gsub("A|T|C|G", "N", x[5:7])
    return(x)
}))
based_vcf[, 2] <- gsub(" ", "", based_vcf[, 2])



result_deletion <- get_deletion_result_from_vcf_more_loose(based_vcf = based_vcf)
result_insertion <- get_insertion_result_from_vcf(based_vcf = based_vcf)

colnames(result_snp) <- colnames(result_deletion)
result <- rbind(rbind(result_snp, result_deletion), result_insertion)
result <- merge_variants(data = result)
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
