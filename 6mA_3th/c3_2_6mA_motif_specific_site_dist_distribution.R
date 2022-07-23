#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(ggplot2)
library(patchwork)

# 2. functions ------------------------------------------------------------ TODO:
plot_motif_specific_site_6mA_dist_distribution <- function(site = NULL, direction = NULL, lower_limit = NULL, upper_limit = NULL, xintercept = NULL) {
    cmd <- ""
    if (direction) {
        cmd <- paste0("cat ", motif_site_dir, "/site", site, "/*.peak.bed", "|", bedtools, " intersect -a stdin -b ", m6mA_bed, "|awk '{print $1\"\t\"$2\"\t\"$3\"\t.\t.\t\"$4}'|", bedtools, " closest -D a -io -a stdin -b ", m6mA_bed, " | awk '{print $NF}'")
    } else {
        cmd <- paste0("cat ", motif_site_dir, "/site", site, "/*.peak.bed", "|", bedtools, " intersect -a stdin -b ", m6mA_bed, "|awk '{print $1\"\t\"$2\"\t\"$3\"\t.\t.\t\"$4}'|", bedtools, " closest -d -io -a stdin -b ", m6mA_bed, " | awk '{print $NF}'")
    }

    dist <- unlist(fread(cmd = cmd, stringsAsFactors = FALSE))
    freq <- as.data.frame(table(dist))
    global <- ggplot(data = freq, aes(x = dist, y = Freq)) +
        geom_line(aes(group = 1)) +
        geom_vline(xintercept = as.character(xintercept), color = "red") +
        geom_vline(xintercept = as.character(0 - xintercept), color = "red") +
        theme(axis.text.x = element_text(size = 5, color = "black", angle = 45, hjust = 0.5, vjust = 0.5))
    subset_p <- ggplot(data = freq, aes(x = dist, y = Freq)) +
        geom_point() +
        geom_line(aes(group = 1)) +
        geom_vline(xintercept = as.character(xintercept), color = "red") +
        geom_vline(xintercept = as.character(0 - xintercept), color = "red") +
        scale_x_discrete(limits = factor(seq(lower_limit, upper_limit)))
    p <- (global | subset_p) + plot_annotation(title = paste0("site", site)) &
        theme(plot.title = element_text(size = 20))
    return(p)
}

plot_overlap_motif_specific_site_6mA_dist_distribution <- function(site = NULL, direction = NULL, lower_limit = NULL, upper_limit = NULL, xintercept = NULL) {
    cmd <- ""
    if (direction) {
        cmd <- paste0("cat ", motif_site_dir, "/overlap_site", site, "/*.peak.bed", "|", bedtools, " intersect -a stdin -b ", m6mA_bed, "|awk '{print $1\"\t\"$2\"\t\"$3\"\t.\t.\t\"$4}'|", bedtools, " closest -D a -io -a stdin -b ", m6mA_bed, " | awk '{print $NF}'")
    } else {
        cmd <- paste0("cat ", motif_site_dir, "/overlap_site", site, "/*.peak.bed", "|", bedtools, " intersect -a stdin -b ", m6mA_bed, "|awk '{print $1\"\t\"$2\"\t\"$3\"\t.\t.\t\"$4}'|", bedtools, " closest -d -io -a stdin -b ", m6mA_bed, " | awk '{print $NF}'")
    }

    dist <- unlist(fread(cmd = cmd, stringsAsFactors = FALSE))
    freq <- as.data.frame(table(dist))
    global <- ggplot(data = freq, aes(x = dist, y = Freq)) +
        geom_line(aes(group = 1)) +
        geom_vline(xintercept = as.character(xintercept), color = "red") +
        geom_vline(xintercept = as.character(0 - xintercept), color = "red") +
        theme(axis.text.x = element_text(size = 5, color = "black", angle = 45, hjust = 0.5, vjust = 0.5))
    subset_p <- ggplot(data = freq, aes(x = dist, y = Freq)) +
        geom_point() +
        geom_line(aes(group = 1)) +
        geom_vline(xintercept = as.character(xintercept), color = "red") +
        geom_vline(xintercept = as.character(0 - xintercept), color = "red") +
        scale_x_discrete(limits = factor(seq(lower_limit, upper_limit)))
    p <- (global | subset_p) + plot_annotation(title = paste0("site", site)) &
        theme(plot.title = element_text(size = 20))
    return(p)
}
# 3. input ---------------------------------------------------------------- TODO:
bedtools <- "/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
m6mA_bed <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/3D7-6mA-9Cell_FRACgt20_COVgt25.bed"
motif_site_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW"

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/m6mA_dist_distri")

pdf("m6mA_site_nonoverlap_2_3_5_6_dist_distri.pdf", width = 20, height = 7)
for (site in c(2, 5, 3, 6)) {
    p <- plot_motif_specific_site_6mA_dist_distribution(
        site = site, direction = TRUE, xintercept = 3,
        lower_limit = -15, upper_limit = 15
    )
    print(p)
}
dev.off()

pdf("m6mA_site_2_3_5_6_dist_distri_no_direc.pdf", width = 20, height = 7)
for (site in c(2, 5, 3, 6)) {
    p <- plot_motif_specific_site_6mA_dist_distribution(
        site = site, direction = FALSE, xintercept = 3,
        lower_limit = -15, upper_limit = 15
    )
    print(p)
}
dev.off()

pdf("m6mA_site_overlap_motif_2_3_5_6_dist_distri.pdf", width = 20, height = 7)
for (site in c(2, 5, 3, 6)) {
    p <- plot_overlap_motif_specific_site_6mA_dist_distribution(
        site = site, direction = TRUE, xintercept = 3,
        lower_limit = -15, upper_limit = 15
    )
    print(p)
}
dev.off()
