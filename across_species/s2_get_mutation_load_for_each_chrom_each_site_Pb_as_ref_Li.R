#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)
library(Biostrings)
library(dplyr)
library(seqinr)
# 2. functions ------------------------------------------------------------

get_based_vcf <- function(vcf = NULL) {
    based_vcf <- t(apply(vcf, 1, function(x) {
        if (x[3] == "N") {
            x[5:7] <- "N"
        } else {
            allele <- c(x[3], unlist(strsplit(x[4], ",")))
            names(allele) <- 0:(length(allele) - 1)
            x[5:7] <- allele[x[5:7]]
        }
        return(x)
    }))
    based_vcf[, 2] <- gsub(" ", "", based_vcf[, 2])

    return(based_vcf)
}

get_indel <- function(based_vcf = NULL, base_class = NULL, outgroup_nogap_index_in_based_vcf = NULL, outgroup_nogap_seq = NULL) {
    standard_strain_seq <- get_species_seq_based_on_ref(
        species_type = "standard_strain", base_class = base_class,
        based_vcf = based_vcf, outgroup_nogap_index_in_based_vcf = outgroup_nogap_index_in_based_vcf
    )
    standard_strain_seq_base_count <- sapply(standard_strain_seq, nchar)
    standard_strain_seq_base_count[standard_strain_seq == "*"] <- 0


    speciesB_seq <- get_species_seq_based_on_ref(
        species_type = "speciesB", base_class = base_class,
        based_vcf = based_vcf, outgroup_nogap_index_in_based_vcf = outgroup_nogap_index_in_based_vcf
    )
    speciesB_seq_base_count <- sapply(speciesB_seq, nchar)
    speciesB_seq_base_count[speciesB_seq == "*"] <- 0

    chrom <- as.character(based_vcf[1, 1])
    # a <- cbind(outgroup_nogap_seq, standard_strain_seq, speciesB_seq)
    deletion_result <- get_deletion_events_in_branch_B(
        outgroup_nogap_seq = outgroup_nogap_seq, speciesB_seq = speciesB_seq, standard_strain_seq = standard_strain_seq
    )
    insertion_result <- get_insertion_events_in_branch_B(
        outgroup_nogap_seq = outgroup_nogap_seq, speciesB_seq_base_count = speciesB_seq_base_count, standard_strain_seq_base_count = standard_strain_seq_base_count
    )

    result <- cbind(chrom, rbind(deletion_result, insertion_result))
    return(result)
}

get_species_seq_based_on_ref <- function(species_type = NULL, base_class = NULL, based_vcf = NULL, outgroup_nogap_index_in_based_vcf = NULL) {
    species_index <- ""
    if (species_type == "outgroupA") {
        species_index <- 6
    } else if (species_type == "standard_strain") {
        species_index <- 5
    } else if (species_type == "speciesB") {
        species_index <- 7
    }
    species_seq <- mclapply(1:length(outgroup_nogap_index_in_based_vcf), function(x) {
        if (base_class[x] == "N_class") {
            return("N")
        } else {
            start_index <- outgroup_nogap_index_in_based_vcf[x]

            if (x == length(outgroup_nogap_index_in_based_vcf)) {
                end_index <- start_index
            } else {
                end_index <- outgroup_nogap_index_in_based_vcf[x + 1] - 1
            }
            seq1 <- based_vcf[start_index:end_index, species_index]
            if (length(seq1) > 1 & any(seq1 != "*")) {
                seq1 <- seq1[seq1 != "*"]
            } else if (length(seq1) > 1 & all(seq1 == "*")) {
                seq1 <- "*"
            }
            result <- paste(seq1, collapse = "")
            return(result)
        }
    }, mc.cores = 4)
    species_seq <- do.call(c, species_seq)
    return(species_seq)
}

get_deletion_events_in_branch_B <- function(outgroup_nogap_seq = NULL, speciesB_seq = NULL, standard_strain_seq = NULL) {
    # outgroupNongap_speciesBGap_index_in_outgroupSeq:要求外群必须得有碱基（A/T/C/G），感兴趣物种（这里是P_reichenowi）必须得缺失
    outgroupNongap_speciesBGap_index_in_outgroupSeq <- which(outgroup_nogap_seq != "*" & speciesB_seq == "*" & outgroup_nogap_seq != "N")
    outgroupNongap_speciesBGap_continuous_index_in_outgroupSeq <- get_continuous_index_region(Pos = outgroupNongap_speciesBGap_index_in_outgroupSeq)
    seq_len <- length(outgroup_nogap_seq)
    deletion_index_in_outgroupSeq <- apply(outgroupNongap_speciesBGap_continuous_index_in_outgroupSeq, 1, function(x) {
        x <- as.numeric(x)
        standard_strain_seq_in_this_region <- standard_strain_seq[x[1]:x[2]]
        standard_strain_exist_base_in_this_region <- length(which(standard_strain_seq_in_this_region != "*")) > 0
        standard_strain_nogap_status_in_del_index <- standard_strain_seq_in_this_region[1] != "*" # 要求标准株在deletion位置必须有碱基，否则无法定位
        if (standard_strain_exist_base_in_this_region & standard_strain_nogap_status_in_del_index) {
            return(x[1])
        } else {
            expand_index <- get_expand_index(x = x, seq_len = seq_len)
            standard_strain_seq_in_expand_index <- standard_strain_seq[expand_index]
            standard_strain_base_count_in_expand_region <- length(which(standard_strain_seq_in_expand_index != "*" & standard_strain_seq_in_expand_index != "N"))
            if (standard_strain_base_count_in_expand_region < length(expand_index) & standard_strain_nogap_status_in_del_index) {
                return(x[1])
            } else {
                return(0)
            }
        }
    })
    deletion_index_in_outgroupSeq <- deletion_index_in_outgroupSeq[deletion_index_in_outgroupSeq != 0]
    deletion_index_in_fasta <- as.numeric(names(outgroup_nogap_seq)[deletion_index_in_outgroupSeq])
    result <- data.frame(bed0 = deletion_index_in_fasta - 1, bed1 = deletion_index_in_fasta, ts = 0, tv = 0, snp = 0, ins = 0, del = 1, indel = 1)
    return(result)
}

get_expand_index <- function(x = NULL, seq_len = NULL) {
    left <- ""
    if (x[1] == 1) {
        left <- NULL
    } else {
        left <- x[1] - 1
    }

    right <- ""
    if (x[2] == seq_len) {
        right <- NULL
    } else {
        right <- x[2] + 1
    }

    return(c(left, right))
}

get_insertion_events_in_branch_B <- function(outgroup_nogap_seq = NULL, speciesB_seq_base_count = NULL, standard_strain_seq_base_count = NULL) {
    # insertion_index_in_outgroupSeq: 要求外群必须得有碱基（A/T/C/G），感兴趣物种（这里是P_reichenowi）对应位点碱基数要大于1，标准株不得缺失（否则无法定位），标准株的碱基数不能等于感兴趣物种（否则变异发生在标准株和感兴趣物种的祖先枝上）
    insertion_index_in_outgroupSeq <- which(outgroup_nogap_seq != "*" & outgroup_nogap_seq != "N" & speciesB_seq_base_count > 1 & standard_strain_seq_base_count != speciesB_seq_base_count & standard_strain_seq_base_count > 0)
    insertion_index_in_fasta <- as.numeric(names(insertion_index_in_outgroupSeq))
    result <- data.frame(bed0 = insertion_index_in_fasta - 1, bed1 = insertion_index_in_fasta, ts = 0, tv = 0, snp = 0, ins = 1, del = 0, indel = 1)
    return(result)
}




get_continuous_index_region <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    region <- t(sapply(difflist, function(x) {
        c(x[1], x[length(x)])
    }))
    rownames(region) <- NULL
    colnames(region) <- c("start", "end")
    return(region)
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
# output <- paste0(dirname(vcf_data), "/", gsub(".vcf", "", basename(vcf_data)), "_only_indel_polysite.bed")
output <- paste0(dirname(vcf_data), "/", gsub(".vcf", "", basename(vcf_data)), "_only_indel_polysite_standardStrainBaseExist.bed")

# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/Pf3D7_01_v3.vcf.include_P_billcollinsi.intermediate"
# 5. process --------------------------------------------------------------
vcf <- fread(vcf_data, sep = "\t", stringsAsFactors = F, header = T, drop = c(3, 6:9))
colnames(vcf) <- c("chrom", "pos", "P_f_allele", "alt", "P_falciparum", "outgroupA", "B_species")

chrom <- as.character(vcf[1, 1])

based_vcf <- get_based_vcf(vcf = vcf)
outgroup_nogap_index_in_based_vcf <- which(based_vcf[, "outgroupA"] != "*") # TODO:ref selection

base_class <- apply(based_vcf[outgroup_nogap_index_in_based_vcf, ], 1, function(x) {
    if (x[6] == "N") {
        return("N_class")
    } else {
        return("nonN_class")
    }
})
outgroup_nogap_seq <- based_vcf[outgroup_nogap_index_in_based_vcf, "outgroupA"] # TODO:ref selection
names(outgroup_nogap_seq) <- based_vcf[outgroup_nogap_index_in_based_vcf, "pos"] # 名字是碱基在fasta文件中的位置

result_indel <- get_indel(
    based_vcf = based_vcf, base_class = base_class,
    outgroup_nogap_index_in_based_vcf = outgroup_nogap_index_in_based_vcf, outgroup_nogap_seq = outgroup_nogap_seq
)

result <- merge_variants(data = result_indel)
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
