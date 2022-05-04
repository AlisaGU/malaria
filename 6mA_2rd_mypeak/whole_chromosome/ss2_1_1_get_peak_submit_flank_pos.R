#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)

# 2. functions ------------------------------------------------------------ TODO:
get_treat_control_level_for_each_site <- function(treatment_level = NULL, control_level = NULL) {
    treat.height.len <- treatment_level$ePos - treatment_level$sPos_minus_1
    treat.height.value <- treatment_level$level
    treat.height <- rep(treat.height.value, treat.height.len)
    treat.height_genome_position <- unlist(Map(`:`, treatment_level$sPos_minus_1 + 1, treatment_level$ePos))

    control.height.len <- control_level$ePos - control_level$sPos_minus_1
    control.height.value <- control_level$level
    control.height <- rep(control.height.value, control.height.len)
    control.height_genome_position <- unlist(Map(`:`, control_level$sPos_minus_1 + 1, control_level$ePos))

    result <- list(
        treat = list(level = treat.height, pos = treat.height_genome_position),
        control = list(level = control.height, pos = control.height_genome_position)
    )
}

select_submit_flank_position <- function(peak_start = NULL, peak_end = NULL,
                                         treatment_level_each_site = NULL, control_level_each_site = NULL,
                                         treatment_pos_each_site = NULL, control_pos_each_site = NULL,
                                         flank_length_type = NULL, flank_length = NULL) {
    # peak_start and peak_end are both 1-based.
    peak_treat_index <- which(treatment_pos_each_site >= peak_start & treatment_pos_each_site <= peak_end)
    peak_treat <- treatment_level_each_site[peak_treat_index]
    names(peak_treat) <- treatment_pos_each_site[peak_treat_index]

    peak_control_index <- which(control_pos_each_site >= peak_start & control_pos_each_site <= peak_end)
    peak_control <- control_level_each_site[peak_control_index]
    names(peak_control) <- control_pos_each_site[peak_control_index]

    flank_specific_length <- get_flank_specific_length(
        peak_start = peak_start, peak_end = peak_end,
        flank_length_type = flank_length_type, flank_length = flank_length
    )

    if (all(names(peak_treat) == names(peak_control))) {
        peak_fold_enrichment <- peak_treat / peak_control
        max_index <- which(peak_fold_enrichment == max(peak_fold_enrichment))
        submit <- names(peak_fold_enrichment)[max_index]
        submit_flank_genome_index <- t(sapply(submit, function(x) {
            get_flank_region(
                index = as.numeric(x), region_start = peak_start,
                region_end = peak_end, flank_specific_length = flank_specific_length
            )
        }))
        return(cbind(submit, submit_flank_genome_index))
    } else {
        print("error")
    }
}

get_flank_specific_length <- function(peak_start = NULL, peak_end = NULL,
                                      flank_length_type = NULL, flank_length = NULL) {
    flank_specific_length <- 0
    if (flank_length_type == "absolute") {
        flank_specific_length <- as.numeric(flank_length)
    } else if (flank_length_type == "relative") {
        peak_length <- peak_end - peak_start + 1
        flank_specific_length <- ceiling(as.numeric(flank_length) * peak_length / 2)
    }
    return(flank_specific_length)
}

get_flank_region <- function(index = NULL, region_start = NULL, region_end = NULL, flank_specific_length = NULL) {
    # start and end are both 1-base
    left_index <- index - flank_specific_length
    right_index <- index + flank_specific_length
    left_index <- ifelse(left_index > region_start, left_index, region_start)
    right_index <- ifelse(right_index < region_end, right_index, region_end)
    return(c(left_index, right_index))
}

get_collapsed_region <- function(start = NULL, end = NULL) {
    a <- unique(unlist(Map(`:`, start, end)))
    rundiff <- c(1, diff(a))
    difflist <- split(a, cumsum(rundiff != 1))
    collapsed_region <- t(sapply(difflist, range))
    collapsed_region[, 1] <- collapsed_region[, 1] - 1
    return(collapsed_region)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
chrom <- Args[6]
peak_pos_filename <- Args[7]
flank_length_type <- Args[8]
flank_length <- Args[9]
treatment_filename <- Args[10]
control_filename <- Args[11]
output <- Args[12]
# 4. variable setting of test module--------------------------------------- TODO:
# chrom <- "Pf3D7_01_v3"
# peak_pos_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/Pf3D7_01_v3.peak.bed"
# flank_length_type <- "absolute"
# flank_length <- "10"
# treatment_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/modifi_level/Pf3D7_01_v3.treat.bed"
# control_filename <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/modifi_level/Pf3D7_01_v3.cont.bed"
# 5. process -------------------------------------------------------------- TODO:
peak_pos <- fread(peak_pos_filename, stringsAsFactors = F, sep = "\t")
colnames(peak_pos) <- c("chrom", "sPos_minus_1", "ePos")
treatment_level <- fread(treatment_filename, stringsAsFactors = F, sep = "\t")
colnames(treatment_level) <- c("chrom", "sPos_minus_1", "ePos", "level")
treatment_level <- treatment_level[treatment_level$ePos > treatment_level$sPos_minus_1, ]
control_level <- fread(control_filename, stringsAsFactors = F, sep = "\t")
colnames(control_level) <- c("chrom", "sPos_minus_1", "ePos", "level")

treat_cont_level_for_each_site <- get_treat_control_level_for_each_site(
    treatment_level = treatment_level,
    control_level = control_level
)
treat_level_each_site <- treat_cont_level_for_each_site$treat$level
treat_pos_each_site <- treat_cont_level_for_each_site$treat$pos
control_level_each_site <- treat_cont_level_for_each_site$control$level
control_pos_each_site <- treat_cont_level_for_each_site$control$pos

peak_submit_flank_pos <- apply(peak_pos, 1, function(x) {
    x <- unlist(x)
    peak_start <- as.numeric(x["sPos_minus_1"]) + 1
    peak_end <- as.numeric(x["ePos"])
    submit_flank_region <- select_submit_flank_position(
        peak_start = peak_start, peak_end = peak_end,
        treatment_level_each_site = treat_level_each_site,
        control_level_each_site = control_level_each_site,
        treatment_pos_each_site = treat_pos_each_site,
        control_pos_each_site = control_pos_each_site,
        flank_length_type = flank_length_type, flank_length = flank_length
    )
})
peak_submit_flank_pos <- do.call(rbind, peak_submit_flank_pos)
collapsed_bed <- get_collapsed_region(
    start = as.numeric(peak_submit_flank_pos[, 2]),
    end = as.numeric(peak_submit_flank_pos[, 3])
)
result <- data.frame(chrom, collapsed_bed)

write.table(result, output, quote = F, col.names = F, row.names = F, sep = "\t")