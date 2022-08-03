#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)
library(parallel)
library(Biostrings)
library(dplyr)
# 2. functions ------------------------------------------------------------
get_derived_types <- function(ancestor = NULL, derived_allele = NULL) {
    deletion <- c(paste(c("A", "T", "C", "G"), "_single_del", sep = ""))

    ts <- c("A2G", "G2A", "T2C", "C2T")
    tv <- c("A2C", "A2T", "G2C", "G2T", "C2A", "C2G", "T2A", "T2G")
    derived_allele_types <- get_type_for_one_base(ancestor = ancestor, derived_allele = derived_allele)

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
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
vcf_data <- Args[6]
output <- paste0(dirname(vcf_data), "/", gsub(".vcf.gz", "", basename(vcf_data)), "_polysite.bed")

bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/Pf3D7_01_v3.species3.vcf.gz"

# 5. process --------------------------------------------------------------
vcf <- fread(cmd = paste0("zcat ", vcf_data, "|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s\",b[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, ncol(vcf) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", "P_falciparum", "B_species")
chrom <- vcf$chrom[1]

vcf_snp_single_deletion <- delete_continuous_deletion(data = vcf)

gt <- transpose(vcf_snp_single_deletion[, 3:4])
mutation_load <- mclapply(1:ncol(gt), function(colth) {
    genome_index <- unlist(vcf_snp_single_deletion[colth, 2])
    x <- unlist(gt[, ..colth])
    allele <- c(x[1], unlist(strsplit(x[2], ",")))
    names(allele) <- 0:(length(allele) - 1)

    ancestor_id <- "0"
    ancestor <- x[1]
    derived_allele <- allele[-which(allele == ancestor)]

    derived_allele_types <- get_derived_types(ancestor = ancestor, derived_allele = derived_allele) ## TODO:
    result <- summary_statics(whole_types = derived_allele_types$whole_types, variant_types = derived_allele_types$variant_types)
    return(result)
}, mc.cores = 4)

mutation_load1 <- do.call(rbind, mutation_load)
result <- data.frame(chrom = chrom, bed0 = vcf_snp_single_deletion[, 2] - 1, bed1 = vcf_snp_single_deletion[, 2], mutation_load = mutation_load1)
colnames(result)[4:9] <- c("transition_load", "transversion_load", "snp_load", "insertion_load", "deletion_load", "indel_load")
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
