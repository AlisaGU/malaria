#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
core_bed <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7.core.bed", as.is = T)

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
chr_core_len <- sapply(unique(core_bed[, 1]), function(x) {
    core.x <- core_bed[core_bed$V1 == x, ]
    core.x.len <- sum(apply(core.x, 1, function(y) {
        as.numeric(y[3]) - as.numeric(y[2])
    }))
})
chr_mutation_count <- sapply(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"), function(x) {
    bed <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_PASS_", x, "_polysite.bed")
    cmd <- paste0("awk '{print $4}' ", bed, "|awk '{s+=$1} END {print s}'")
    result <- as.numeric(system(cmd, intern = T))
})
chr_site_count <- sapply(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"), function(x) {
    bed <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_PASS_", x, "_polysite.bed")
    cmd <- paste0("wc -l ", bed, "|awk '{print $1}'")
    result <- as.numeric(system(cmd, intern = T))
})

chr_pass_site_count <- sapply(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"), function(x) {
    filename <- paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_PASS_", x, ".vcf.gz")
    cmd <- paste0("zcat ", filename, "|grep -v \"^#\"|wc -l")
    result <- as.numeric(system(cmd, intern = T))
})

chr_mutation_density <- chr_mutation_count / chr_core_len
chr_site_density <- chr_site_count / chr_core_len
chr_pass_site_density <- chr_pass_site_count / chr_core_len

data_for_plot <- data.frame(
    group = c(rep("mutation", 14), rep("polysite", 14), rep("pass", 14)),
    chr = c(1:14, 1:14, 1:14),
    value = c(chr_mutation_density, chr_site_density, chr_pass_site_density)
)
ggplot(data_for_plot) +
    geom_line(aes(x = chr, y = value, color = group)) +
    scale_x_continuous(breaks = 1:14, labels = 1:14)



before <- as.numeric(unlist(fread(cmd = paste0("zcat /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_snp_01.vcf.gz|grep -v \"^#\"|awk '{print $2}'"), stringsAsFactors = T)))

after <- as.numeric(unlist(fread("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6/WSEA_PASS_01_snp_polysite.bed", select = 3, stringsAsFactors = F)))

after.unique <- setdiff(after, before)