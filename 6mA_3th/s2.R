#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(pheatmap)
# 2. functions ------------------------------------------------------------ TODO:
read_m6mA_status_peak_data <- function() {
    result <- lapply(autosomes, function(autosome) {
        data <- cbind(autosome, fread(paste0(m6mA_status_dir, "/", autosome, ".peak.txt"), header = F, stringsAsFactors = F))
        colnames(data) <- c("chrom", paste("left_flank_", 10:1, sep = ""), paste("motif_site_", 1:6, sep = ""), paste("right_flank_", 1:10, sep = ""))
        return(data)
    })
    result <- do.call(rbind, result)
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
m6mA_status_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_each_site/m6mA_status_motif_and_flank_site"
motif_bed_dir <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/motif_bed"
autosomes <- c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
m6mA_status_peak_data <- read_m6mA_status_peak_data()
m6mA_status_peak_data[m6mA_status_peak_data == -1] <- NA
a <- m6mA_status_peak_data[, -1]
m6mA_exist_data <- m6mA_status_peak_data[apply(a, 1, function(x) {
    x <- unlist(x)
    return(any(x[!is.na(x)] == 1))
}), ]


data <- m6mA_exist_data[, -1]

annotation_col = data.frame(SiteClass = factor(rep(c("Left_flank", "Motif", "Right_flank"), c(10, 6, 10))))
rownames(annotation_col) <- colnames(data)

pheatmap(data,
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = F, color = c("white", "red"),
    annotation_col = annotation_col, legend = FALSE,
    annotation_colors = list(SiteClass = c(Left_flank = "#31b679", Motif = "#ffe725", Right_flank = "#31b679")),
    filename = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_each_site/m6mA_status_motif_and_flank_site.pdf"
)



data_f3 <- m6mA_status_peak_data[m6mA_status_peak_data$left_flank_3 == 1 | m6mA_status_peak_data$right_flank_3 == 1 | left_flank_3 == 1 | m6mA_status_peak_data$right_flank_3 == 1, -1]
pheatmap(data_f3,
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = F, color = c("white", "red"),
    annotation_col = annotation_col, legend = FALSE,
    annotation_colors = list(SiteClass = c(Left_flank = "#31b679", Motif = "#ffe725", Right_flank = "#31b679")),
    filename = "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public/motif_each_site/m6mA_status_f3.pdf"
)

data_f3_noNA <- m6mA_exist_data[!is.na(m6mA_exist_data$left_flank_3) | !is.na(m6mA_exist_data$right_flank_3), -1]
