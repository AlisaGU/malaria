#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(plyr)

# 2. functions ------------------------------------------------------------ TODO:
get_closest_distance <- function(data = NULL) {
    min_index <- which(data$motif_flank_distance == min(data$motif_flank_distance))[1]
    result <- data[min_index, ]
    return(result)
}

write_flank <- function(data = NULL) {
    flank_site <- as.character(data$motif_flank_distance[1])
    direction <- gsub("flank_", "", basename(output_global_dir))
    output_dir <- paste0(output_global_dir, "/flank_", flank_site)
    system(command = paste0("mkdir -p ", output_dir))
    write.table(data[, 7:9], file = paste0(output_dir, "/", peak_type, ".allChroms.bed"), sep = "\t", row.names = F, col.names = F, append = T, quote = F)
}
# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
cmd <- Args[6]
output_global_dir <- Args[7]
peak_type <- Args[8]
# 4. variable setting of test module--------------------------------------- TODO:
# cmd <- "sort -k1,1 -k2,2n /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/GAMSAA/flank_upstream/flank_1/peak.allChroms.upstream.flank1.bed.intermediate | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a stdin -b /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/core_consistent_peak.bed | /picb/evolgen/users/gushanshan/software/bedtools/bedtools subtract -a stdin -b /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/all_motifs_position.peak.bed |/picb/evolgen/users/gushanshan/software/bedtools/bedtools closest -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/all_motifs_position.peak.bed -b stdin -d -D \"a\" -id|awk '$7~/^P/'"
# output_global_dir<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/GAMSAA/flank_upstream"
# peak_type<-"peak"
# 5. process -------------------------------------------------------------- TODO:
cat(cmd)
data <- fread(cmd = cmd, stringsAsFactors = T)
colnames(data) <- c(
    "chrom_flank_site", "flank_index_minus_1", "flank_index", "chrom_motif",
    "motif_start_minus_1", "motif_end", "sth", "sth1", "motif_strand", "motif_flank_distance"
)
data$motif_flank_distance <- abs(data$motif_flank_distance)

result <- ddply(data, .(chrom_flank_site, flank_index_minus_1, flank_index), get_closest_distance)
ddply(result, .(motif_flank_distance), write_flank)
