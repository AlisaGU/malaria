#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
# 2. functions ------------------------------------------------------------ TODO:
get_peak_mutation_site_count_all_variant <- function(chrom = NULL, start = NULL, end = NULL) {
    peak_genome_index <- paste0(chrom, ":", start, "-", end)
    allvar_peak_command <- paste0(bcftools, " view -r ", peak_genome_index, " ", allVariant, " | grep -v \"^#\"|wc -l")
    allvar_peak_mutation_count <- as.numeric(system(allvar_peak_command, intern = T))
    scaled_allvar_peak_mutation_count <- allvar_peak_mutation_count / (end - start + 1)

    surr_genome_region <- matrix(get_around_genome_region(chrom = chrom, start = start, end = end), ncol = 3)
    surr_genome_index <- apply(surr_genome_region, 1, function(x) {
        paste0(x[1], ":", x[2], "-", x[3])
    })
    allvar_surr_mutation_count <- sum(sapply(surr_genome_index, function(x) {
        allvar_surr_command.i <- paste0(bcftools, " view -r ", x, " ", allVariant, " | grep -v \"^#\"|wc -l")
        allvar_surr_mutation_count.i <- as.numeric(system(allvar_surr_command.i, intern = T))
    }))
    scaled_allvar_surr_mutation_count <- allvar_surr_mutation_count / 1000
    return(c(scaled_allvar_peak_mutation_count, scaled_allvar_surr_mutation_count))
}

get_peak_mutation_site_count_bisnp <- function(chrom = NULL, start = NULL, end = NULL) {
    peak_genome_index <- paste0(chrom, ":", start, "-", end)
    bisnp_peak_command <- paste0(bcftools, " view -r ", peak_genome_index, " ", bisnp, " | grep -v \"^#\"|wc -l")
    bisnp_peak_mutation_count <- as.numeric(system(bisnp_peak_command, intern = T))
    scaled_bisnp_peak_mutation_count <- bisnp_peak_mutation_count / (end - start + 1)

    surr_genome_region <- matrix(get_around_genome_region(chrom = chrom, start = start, end = end), ncol = 3)
    surr_genome_index <- apply(surr_genome_region, 1, function(x) {
        paste0(x[1], ":", x[2], "-", x[3])
    })
    bisnp_surr_mutation_count <- sum(sapply(surr_genome_index, function(x) {
        bisnp_surr_command.i <- paste0(bcftools, " view -r ", x, " ", bisnp, " | grep -v \"^#\"|wc -l")
        bisnp_surr_mutation_count.i <- as.numeric(system(bisnp_surr_command.i, intern = T))
    }))
    scaled_bisnp_surr_mutation_count <- bisnp_surr_mutation_count / 1000
    return(c(scaled_bisnp_peak_mutation_count, scaled_bisnp_surr_mutation_count))
}

get_around_genome_region <- function(chrom = NULL, start = NULL, end = NULL) {
    # start and end are both 1-base
    chrom_length.i <- chrom_length$length[which(chrom_length$chrom == chrom)]
    left_index <- start - 500
    right_index <- end + 500
    around_index <- c()

    if (start == 1) {
        around_index <- c(chrom, end + 1, end + 1000)
    }

    if (left_index >= 1 && right_index <= chrom_length.i) {
        around_index <- rbind(c(chrom, left_index, start - 1), c(chrom, end + 1, right_index))
    }
    if (left_index <= 0 & start != 1) {
        left_borrow_len <- length(left_index:0)
        around_index <- rbind(c(chrom, 1, start - 1), c(chrom, end + 1, right_index + left_borrow_len))
    }
    if (right_index > chrom_length.i) {
        right_borrow_len <- length((chrom_length.i + 1):right_index)
        around_index <- rbind(c(chrom, left_index - right_borrow_len, start - 1), c(chrom, end + 1, chrom_length.i))
    }

    return(around_index)
}
# 3. input ---------------------------------------------------------------- TODO:
allVariant <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_peak500_surr.vcf.gz"
bisnp <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_peak500_surr_bisnp.vcf.gz"
bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
top500_peak_index <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/top500/top500.bed", as.is = T)
chrom_length <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length", stringsAsFactors = F)
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
colnames(top500_peak_index) <- c("chrom", "start", "end")
colnames(chrom_length) <- c("chrom", "length")
top500_peak_index$start <- top500_peak_index$start + 1

peak_mutation_count <- t(apply(top500_peak_index, 1, function(x) {
    chrom <- unlist(x[1])
    start <- as.numeric(x[2])
    end <- as.numeric(x[3])
    allvar_count <- get_peak_mutation_site_count_all_variant(chrom = chrom, start = start, end = end)
    bisnp_count <- get_peak_mutation_site_count_bisnp(chrom = chrom, start = start, end = end)
    return(c(allvar_count, bisnp_count))
}))
colnames(peak_mutation_count) <- c("allvar_peak", "allvar_surr", "bisnp_peak", "bisnp_surr")
rownames(peak_mutation_count) <- paste(paste(top500_peak_index[, 1], top500_peak_index[, 2], sep = ":"), top500_peak_index[, 3], sep = "-")
write.table(peak_mutation_count, "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif/peak_mutation.txt", quote = F, sep = "\t")

peak_count <- nrow(peak_mutation_count)
data_for_plot <- data.frame(
    mutation_type = c(rep("allvar", peak_count * 2), rep("bisnp", peak_count * 2)),
    genome_region = c(
        rep("peak", peak_count), rep("surr", peak_count),
        rep("peak", peak_count), rep("surr", peak_count)
    ), value = c(
        peak_mutation_count[, 1],
        peak_mutation_count[, 2],
        peak_mutation_count[, 3],
        peak_mutation_count[, 4]
    )
)


ggboxplot(data = data_for_plot, x = "mutation_type", y = "value", fill = "genome_region") +
    stat_compare_means(aes(group = genome_region), label = "p.format") +
    scale_fill_manual(values = c(peak = "#efc000", surr = "#868686")) +
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
ggsave("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_motif/mutation_count_compar_between_peak_and_surr.pdf", width = 7, height = 7)