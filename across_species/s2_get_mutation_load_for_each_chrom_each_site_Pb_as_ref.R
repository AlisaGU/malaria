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
        allele <- c(x[3], unlist(strsplit(x[4], ",")))
        names(allele) <- 0:(length(allele) - 1)
        x[5:7] <- allele[x[5:7]]
        return(x)
    }))
    based_vcf[, 2] <- gsub(" ", "", based_vcf[, 2])
    return(based_vcf)
}

get_indel_from_vcf_more_loose <- function(based_vcf = NULL, ref_nogap_index = NULL, ref_nogap_seq = NULL) {
    standard_strain_seq <- get_species_seq_based_on_ref(
        species_type = "standard_strain",
        based_vcf = based_vcf, ref_nogap_index = ref_nogap_index
    )
    standard_strain_seq_base_count <- sapply(standard_strain_seq, nchar)
    standard_strain_seq_base_count[standard_strain_seq == "*"] <- 0
    speciesB_seq <- get_species_seq_based_on_ref(
        species_type = "speciesB",
        based_vcf = based_vcf, ref_nogap_index = ref_nogap_index
    )
    speciesB_seq_base_count <- sapply(speciesB_seq, nchar)
    speciesB_seq_base_count[speciesB_seq == "*"] <- 0

    chrom <- as.character(based_vcf[1, 1])
    # a <- cbind(ref_nogap_seq, standard_strain_seq, speciesB_seq)
    deletion_result <- get_deletion_events_in_branch_B(
        standard_strain_seq_base_count = standard_strain_seq_base_count,
        speciesB_seq_base_count = speciesB_seq_base_count, speciesB_seq = speciesB_seq, ref_nogap_seq = ref_nogap_seq
    )
    insertion_result <- get_insertion_events_in_branch_B(
        standard_strain_seq_base_count = standard_strain_seq_base_count,
        speciesB_seq_base_count = speciesB_seq_base_count, ref_nogap_seq = ref_nogap_seq
    )

    result <- cbind(chrom, rbind(deletion_result, insertion_result))
    return(result)
}

get_species_seq_based_on_ref <- function(species_type = NULL, based_vcf = NULL, ref_nogap_index = NULL) {
    species_index <- ""
    if (species_type == "outgroupA") {
        species_index <- 6
    } else if (species_type == "standard_strain") {
        species_index <- 5
    } else if (species_type == "speciesB") {
        species_index <- 7
    }
    species_seq <- mclapply(1:length(ref_nogap_index), function(x) {
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
    }, mc.cores = 4)
    species_seq <- do.call(c, species_seq)
    return(species_seq)
}

get_deletion_events_in_branch_B <- function(standard_strain_seq_base_count = NULL, speciesB_seq_base_count = NULL, speciesB_seq = NULL, ref_nogap_seq = NULL) {
    speciesB_deletion_index <- as.numeric(get_continuous_index_region(Pos = as.numeric(based_vcf[based_vcf[, "B_species"] == "*", 2]))[, 1])

    deletion_index_raw <- as.numeric(names(ref_nogap_seq)[which(speciesB_seq_base_count == 0 & standard_strain_seq_base_count > 0)])
    final_deletion_index <- intersect(speciesB_deletion_index, deletion_index_raw)

    result <- data.frame(bed0 = final_deletion_index - 1, bed1 = final_deletion_index, ts = 0, tv = 0, snp = 0, ins = 0, del = 1, indel = 1)
    return(result)
}

get_insertion_events_in_branch_B <- function(standard_strain_seq_base_count = NULL, speciesB_seq_base_count = NULL, ref_nogap_seq = NULL) {
    insertion_index <- which(speciesB_seq_base_count > standard_strain_seq_base_count & standard_strain_seq_base_count == 1)
    fasta_pos <- as.numeric(names(ref_nogap_seq)[insertion_index])
    result <- data.frame(bed0 = fasta_pos - 1, bed1 = fasta_pos, ts = 0, tv = 0, snp = 0, ins = 1, del = 0, indel = 1)
    return(result)
}




get_continuous_index_region <- function(Pos = NULL) {
    rundiff <- c(1, diff(Pos))
    difflist <- split(Pos, cumsum(rundiff != 1))
    region <- t(sapply(difflist, function(x) {
        c(x[1], x[length(x)])
    }))
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
output <- paste0(dirname(vcf_data), "/", gsub(".vcf", "", basename(vcf_data)), "_only_indel_polysite.bed")
bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# vcf_data <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf/Pf3D7_01_v3.species3.vcf.gz"

# 5. process --------------------------------------------------------------
vcf <- fread(cmd = paste0("cat ", vcf_data, "|grep -v \"^#\"|awk 'BEGIN{OFS=\"\t\"}{for(i=1;i<=2;i++){printf(\"%s%s\",$i,OFS)}for(i=4;i<=5;i++){printf(\"%s%s\",$i,OFS)}for(i=10;i<=NF;i++){split($i,a,\":\");split(a[1],b,/[|/]/);printf(\"%s%s\",b[1],OFS)};printf(ORS)}'"), sep = "\t", stringsAsFactors = F, header = F)
vcf[, c(8) := NULL]
# vcf[, ncol(vcf) := NULL]
colnames(vcf) <- c("chrom", "pos", "ref", "alt", "P_falciparum", "outgroupA", "B_species")

chrom <- as.character(vcf[1, 1])

based_vcf <- get_based_vcf(vcf = vcf)
ref_nogap_index <- which(based_vcf[, 6] != "*") # TODO:ref selection
ref_nogap_seq <- based_vcf[ref_nogap_index, 6] # TODO:ref selection
names(ref_nogap_seq) <- based_vcf[ref_nogap_index, 2]

result_indel <- get_indel_from_vcf_more_loose(
    based_vcf = based_vcf,
    ref_nogap_index = ref_nogap_index, ref_nogap_seq = ref_nogap_seq
)

result <- merge_variants(data = result_indel)
fwrite(result, output, scipen = 10, sep = "\t", row.names = F, quote = F, col.names = F)
