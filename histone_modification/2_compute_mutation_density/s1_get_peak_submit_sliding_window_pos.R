#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(evobiR)
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

select_submit_sliding_window_position <- function(peak_start = NULL, peak_end = NULL,
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
        sliding_window_genome_index <- lapply(submit, function(x) {
            get_sliding_window(
                index = as.numeric(x), region_start = peak_start,
                region_end = peak_end, flank_specific_length = flank_specific_length
            )
        })
        sliding_window_genome_index <- do.call(rbind, sliding_window_genome_index)
        return(sliding_window_genome_index)
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

get_sliding_window <- function(index = NULL, region_start = NULL,
                               region_end = NULL, flank_specific_length = NULL) {
    left_seq_region <- rev(region_start:(index - 1))
    # right_seq_region <- (index + 1):region_end
    right_seq_region <- index:region_end # 把submit划进右侧第一个窗口

    left_sliding_window <- data.frame(
        left_index = mySlidingWindow_min(data = left_seq_region, window = flank_specific_length, step = flank_specific_length),
        right_index = mySlidingWindow_max(data = left_seq_region, window = flank_specific_length, step = flank_specific_length)
    )
    if (nrow(left_sliding_window) > 0) {
        left_sliding_window$window <- paste0("w", 1:NROW(left_sliding_window))
    }

    right_sliding_window <- data.frame(
        left_index = mySlidingWindow_min(data = right_seq_region, window = flank_specific_length, step = flank_specific_length),
        right_index = mySlidingWindow_max(data = right_seq_region, window = flank_specific_length, step = flank_specific_length)
    )
    if (nrow(right_sliding_window) > 0) {
        right_sliding_window$window <- paste0("w", 1:NROW(right_sliding_window))
    }

    result <- data.frame(
        left_index = c(left_sliding_window$left_index, right_sliding_window$left_index),
        right_index = c(left_sliding_window$right_index, right_sliding_window$right_index),
        window = c(left_sliding_window$window, right_sliding_window$window)
    )
    result <- result[order(result$window), ]
    return(result)
}

mySlidingWindow_min <- function(data = NULL, window = NULL, step = NULL) {
    if (length(data) <= window) {
        return(NULL)
    } else {
        return(SlidingWindow("min", data, window, step))
    }
}

mySlidingWindow_max <- function(data = NULL, window = NULL, step = NULL) {
    if (length(data) <= window) {
        return(NULL)
    } else {
        return(SlidingWindow("max", data, window, step))
    }
}

get_collapsed_region_and_write <- function(peak_submit_sliding_window_pos = NULL, chrom = NULL, output_dir = NULL) {
    data_by_window <- split(peak_submit_sliding_window_pos, f = as.factor(peak_submit_sliding_window_pos$window), drop = T)
    a <- sapply(data_by_window, function(x) {
        start <- x[, 1]
        end <- x[, 2]
        window <- x[1, 3]
        a <- sort(unique(unlist(Map(`:`, start, end))))
        rundiff <- c(1, diff(a))
        difflist <- split(a, cumsum(rundiff != 1))
        collapsed_region <- t(sapply(difflist, range))
        collapsed_region[, 1] <- collapsed_region[, 1] - 1
        result <- data.frame(chrom = chrom, start = collapsed_region[, 1], end = collapsed_region[, 2])
        window <- x[1, 3]
        write.table(result, paste0(output_dir, "/", window, ".bed"), quote = F, col.names = F, row.names = F, sep = "\t")
    })
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
working_dir <- Args[6]
chrom <- Args[7]
peak_pos_filename <- Args[8]
flank_length_type <- Args[9]
flank_length <- Args[10]
treatment_filename <- Args[11]
control_filename <- Args[12]
output_dir <- Args[13]
# 4. variable setting of test module--------------------------------------- TODO:
# working_dir<-"/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/histone_methy/WT-42-46h-H3K4me3/BAM_broad_nomodel_dup1"
# chrom <- "Pf3D7_01_v3"
# peak_pos_filename <- "chrom/Pf3D7_01_v3.peak.bed"
# flank_length_type <- "relative"
# flank_length <- "0.05"
# treatment_filename <- "modifi_level/Pf3D7_01_v3.treat.bed"
# control_filename <- "modifi_level/Pf3D7_01_v3.cont.bed"
# output_dir<-"chrom/sliding_window/relative.0.05/Pf3D7_01_v3"
# 5. process -------------------------------------------------------------- TODO:
setwd(working_dir)
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

peak_submit_sliding_window_pos <- apply(peak_pos, 1, function(x) {
    x <- unlist(x)
    peak_start <- as.numeric(x["sPos_minus_1"]) + 1
    peak_end <- as.numeric(x["ePos"])
    submit_sliding_window_region <- select_submit_sliding_window_position(
        peak_start = peak_start, peak_end = peak_end,
        treatment_level_each_site = treat_level_each_site,
        control_level_each_site = control_level_each_site,
        treatment_pos_each_site = treat_pos_each_site,
        control_pos_each_site = control_pos_each_site,
        flank_length_type = flank_length_type, flank_length = flank_length
    )
})

# for (i in 1:nrow(peak_pos)) {
#     x <- peak_pos[i, ]
#     x <- unlist(x)
#     peak_start <- as.numeric(x["sPos_minus_1"]) + 1
#     peak_end <- as.numeric(x["ePos"])
#     submit_sliding_window_region <- select_submit_sliding_window_position(
#         peak_start = peak_start, peak_end = peak_end,
#         treatment_level_each_site = treat_level_each_site,
#         control_level_each_site = control_level_each_site,
#         treatment_pos_each_site = treat_pos_each_site,
#         control_pos_each_site = control_pos_each_site,
#         flank_length_type = flank_length_type, flank_length = flank_length
#     )
# }
peak_submit_sliding_window_pos <- do.call(rbind, peak_submit_sliding_window_pos)
get_collapsed_region_and_write(peak_submit_sliding_window_pos = peak_submit_sliding_window_pos, chrom = chrom, output_dir = output_dir)
