#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(data.table)
# 2. functions ------------------------------------------------------------ TODO:
get_whole_chrom_level <- function() {

}
select_submit_flank_position <- function(peak_start = NULL, peak_end = NULL, chrom_treatment_height = NULL, chrom_control_height = NULL, chrom_treatment_genome_index = NULL, chrom_control_genome_index = NULL) {
    peak_treat_index <- which(chrom_treatment_genome_index >= peak_start & chrom_treatment_genome_index <= peak_end)
    peak_treat <- chrom_treatment_height[peak_treat_index]
    names(peak_treat) <- chrom_treatment_genome_index[peak_treat_index]

    peak_control_index <- which(chrom_control_genome_index >= peak_start & chrom_control_genome_index <= peak_end)
    peak_control <- chrom_control_height[peak_control_index]
    names(peak_control) <- chrom_control_genome_index[peak_control_index]

    if (all(names(peak_treat) == names(peak_control))) {
        peak_fold_enrichment <- peak_treat / peak_control
        max_index <- which(peak_fold_enrichment == max(peak_fold_enrichment))[1]
        submit <- names(peak_fold_enrichment)[max_index]
        submit_flank_genome_index <- get_flank_region(index = as.numeric(submit), region_start = peak_start, region_end = peak_end)
        return(c(submit, submit_flank_genome_index))
    } else {
        print("error")
    }
}

get_flank_region <- function(index = NULL, region_start = NULL, region_end = NULL) {
    # start and end are both 1-base
    left_index <- index - 10
    right_index <- index + 10
    left_index <- ifelse(left_index > region_start, left_index, region_start)
    right_index <- ifelse(right_index < region_end, right_index, region_end)
    return(c(left_index, right_index))
}

get_peak_mutation_count <- function(chrom = NULL, peak_submit_region_start = NULL, peak_submit_region_end = NULL, variantfile = NULL) {
    peak_genome_index <- paste0(chrom, ":", peak_submit_region_start - 1, "-", peak_submit_region_end)
    command <- paste0(bcftools, " view -r ", peak_genome_index, " ", variantfile, " | grep -v \"^#\"|awk -F\"\\t\" '$2>=", peak_submit_region_start, " && $2<=", peak_submit_region_end, "'|wc -l")
    peak_mutation_count <- as.numeric(system(command, intern = T))
    scaled_peak_mutation_count <- peak_mutation_count / (peak_submit_region_end - peak_submit_region_start + 1)
}


get_mean_nopeak_mutation_count_from_sample20 <- function(chrom = NULL, nopeak_region_start = NULL, nopeak_region_end = NULL, peak_length = NULL, variantfile = NULL) {
    mutation <- sapply(1:20, function(seed) {
        set.seed(seed)
        nopeak_start.seed <- sample((nopeak_region_start + 1):(nopeak_region_end - peak_length + 1), 1)
        nopeak_end.seed <- nopeak_start.seed + peak_length

        nopeak_genome_index <- paste0(chrom, ":", nopeak_start.seed, "-", nopeak_end.seed)
        command <- paste0(bcftools, " view -r ", nopeak_genome_index, " ", variantfile, " | grep -v \"^#\"|awk -F\"\\t\" '$2>=", nopeak_start.seed, " && $2<=", nopeak_end.seed, "'|wc -l")
        nopeak_mutation_count <- as.numeric(system(command, intern = T))
        return(nopeak_mutation_count)
    })
    mean_mutation <- mean(mutation) / peak_length
    return(mean_mutation)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
region_peak_nopeak_filename <- Args[6]
allVariant_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
treadment_bdg_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_treat_pileup.bdg"
control_bdg_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_control_lambda.bdg"
# 4. variable setting of test module--------------------------------------- TODO:
region_peak_nopeak_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/intergenic_nopeak_closest_to_peak.bed"
# 5. process -------------------------------------------------------------- TODO:

region_peak_nopeak <- fread(region_peak_nopeak_filename, stringsAsFactors = F, select = c(1:3, 7:9))
colnames(region_peak_nopeak) <- c("peak_chrom", "peak_start", "peak_end", "nopeak_chrom", "nopeak_start", "nopeak_end")
treatment_bdg <- fread(treadment_bdg_filename, skip = 1, stringsAsFactors = F, sep = "\t")
colnames(treatment_bdg) <- c("chrom", "start", "end", "value")
treatment_bdg <- treatment_bdg[which(treatment_bdg$end > treatment_bdg$start), ]
control_bdg <- fread(control_bdg_filename, , skip = 1, stringsAsFactors = F, sep = "\t")
colnames(control_bdg) <- c("chrom", "start", "end", "value")
control_bdg <- control_bdg[which(control_bdg$end > control_bdg$start), ]

chroms <- unique(region_peak_nopeak$peak_chrom)
for (chrom in chroms) {
    chrom_index_in_treat <- which(treatment_bdg$chrom == chrom)
    treat_bdg_chrom <- treatment_bdg[chrom_index_in_treat, ]
    treat.height.len <- treat_bdg_chrom$end - treat_bdg_chrom$start
    treat.height.value <- treat_bdg_chrom$value
    treat.height <- rep(treat.height.value, treat.height.len)
    treat.height_genome_position <- unlist(Map(`:`, treat_bdg_chrom$start + 1, treat_bdg_chrom$end))

    chrom_index_in_control <- which(control_bdg$chrom == chrom)
    control_bdg_chrom <- control_bdg[chrom_index_in_control, ]
    control.height.len <- control_bdg_chrom$end - control_bdg_chrom$start
    control.height.value <- control_bdg_chrom$value
    control.height <- rep(control.height.value, control.height.len)
    control.height_genome_position <- unlist(Map(`:`, control_bdg_chrom$start + 1, control_bdg_chrom$end))

    peak_chrom_data <- region_peak_nopeak[region_peak_nopeak$peak_chrom == chrom, ]
    submit_flank <- t(apply(peak_chrom_data, 1, function(x) {
        select_submit_flank_position(
            peak_start = as.numeric(x[2]) + 1, peak_end = as.numeric(x[3]),
            chrom_treatment_height = treat.height, chrom_control_height = control.height, chrom_treatment_genome_index = treat.height_genome_position,
            chrom_control_genome_index = control.height_genome_position
        )
    }))
    peak_chrom_data$peak_submit <- as.numeric(submit_flank[, 1])
    peak_chrom_data$peak_submit_lf <- as.numeric(submit_flank[, 2])
    peak_chrom_data$peak_submit_rf <- as.numeric(submit_flank[, 3])
    chrom_number <- unlist(strsplit(chrom, "_"))[2]
    all_variant <- paste0(allVariant_dir, "/WSEA_", chrom_number, ".vcf.gz")
    bisnp <- paste0(allVariant_dir, "/WSEA_bisnp_", chrom_number, ".vcf.gz")

    # mutation_count <- t(apply(peak_chrom_data, 1, function(x) {
    #     peak_submit_region_start <- as.numeric(x[8])
    #     peak_submit_region_end <- as.numeric(x[9])
    #     nopeak_region_start <- as.numeric(x[5])
    #     nopeak_region_end <- as.numeric(x[6])
    #     peak_length <- peak_submit_region_end - peak_submit_region_start + 1

    #     allvar_peak_mutation_count <- get_peak_mutation_count(
    #         chrom = chrom,
    #         peak_submit_region_start = peak_submit_region_start,
    #         peak_submit_region_end = peak_submit_region_end,
    #         variantfile = all_variant
    #     )
    #     bisnp_peak_mutation_count <- get_peak_mutation_count(
    #         chrom = chrom,
    #         peak_submit_region_start = peak_submit_region_start,
    #         peak_submit_region_end = peak_submit_region_end,
    #         variantfile = bisnp
    #     )

    #     allvar_nopeak_mutation_count <- get_mean_nopeak_mutation_count_from_sample20(
    #         chrom = chrom,
    #         nopeak_region_start = nopeak_region_start,
    #         nopeak_region_end = nopeak_region_end,
    #         peak_length = peak_length,
    #         variantfile = all_variant
    #     )
    #     bisnp_nopeak_mutation_count <- get_mean_nopeak_mutation_count_from_sample20(
    #         chrom = chrom,
    #         nopeak_region_start = nopeak_region_start,
    #         nopeak_region_end = nopeak_region_end,
    #         peak_length = peak_length,
    #         variantfile = bisnp
    #     )

    #     result <- c(allvar_peak_mutation_count, allvar_nopeak_mutation_count, bisnp_peak_mutation_count, bisnp_nopeak_mutation_count)
    # }))

    mutation_count <- c()
    for (i in 1:nrow(peak_chrom_data)) {
        print(i)
        x <- unlist(peak_chrom_data[i, ])
        peak_submit_region_start <- as.numeric(x[8])
        peak_submit_region_end <- as.numeric(x[9])
        nopeak_region_start <- as.numeric(x[5])
        nopeak_region_end <- as.numeric(x[6])
        peak_length <- peak_submit_region_end - peak_submit_region_start + 1

        allvar_peak_mutation_count <- get_peak_mutation_count(
            chrom = chrom,
            peak_submit_region_start = peak_submit_region_start,
            peak_submit_region_end = peak_submit_region_end,
            variantfile = all_variant
        )
        bisnp_peak_mutation_count <- get_peak_mutation_count(
            chrom = chrom,
            peak_submit_region_start = peak_submit_region_start,
            peak_submit_region_end = peak_submit_region_end,
            variantfile = bisnp
        )

        allvar_nopeak_mutation_count <- get_mean_nopeak_mutation_count_from_sample20(
            chrom = chrom,
            nopeak_region_start = nopeak_region_start,
            nopeak_region_end = nopeak_region_end,
            peak_length = peak_length,
            variantfile = all_variant
        )
        bisnp_nopeak_mutation_count <- get_mean_nopeak_mutation_count_from_sample20(
            chrom = chrom,
            nopeak_region_start = nopeak_region_start,
            nopeak_region_end = nopeak_region_end,
            peak_length = peak_length,
            variantfile = bisnp
        )

        result <- c(allvar_peak_mutation_count, allvar_nopeak_mutation_count, bisnp_peak_mutation_count, bisnp_nopeak_mutation_count)
        mutation_count <- rbind(mutation_count, result)
    }
    write.table(mutation_count, "mutation_count.txt", sep = "\t", row.names = F, col.names = F, quote = F)
    colnames(mutation_count) <- c("allvar_peak", "allvar_surr", "bisnp_peak", "bisnp_surr")
    peak_count <- nrow(mutation_count)
    data_for_plot <- data.frame(
        mutation_type = c(rep("allvar", peak_count * 2), rep("bisnp", peak_count * 2)),
        genome_region = c(
            rep("peak", peak_count), rep("nopeak", peak_count),
            rep("peak", peak_count), rep("nopeak", peak_count)
        ), value = c(
            mutation_count[, 1],
            mutation_count[, 2],
            mutation_count[, 3],
            mutation_count[, 4]
        )
    )

    ggboxplot(data = data_for_plot, x = "mutation_type", y = "value", fill = "genome_region") +
        stat_compare_means(aes(group = genome_region), label = "p.format", method = "wilcox.test", paired = T) +
        scale_fill_manual(values = c(peak = "#efc000", nopeak = "#868686")) +
        labs(y = "Scale mutation count") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 20, color = "black"),
            legend.title = element_blank()
        )
    ggsave("a.pdf")
}