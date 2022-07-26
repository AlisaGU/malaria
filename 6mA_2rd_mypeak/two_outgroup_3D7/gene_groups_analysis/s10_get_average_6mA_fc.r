#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table, quietly = T)
library(dplyr, warn.conflicts = FALSE)
# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
window_treat_intensity_cmd <- Args[6]
window_cont_intensity_cmd <- Args[7]

# 4. variable setting of test module--------------------------------------- TODO:
# window_treat_intensity_cmd <- "grep \"\\_2$\" /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/include_pseudo/DR/DR_genes.window.bed | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_treat_pileup.bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$8}'"
# window_cont_intensity_cmd<-"grep \"\\_2$\" /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/include_pseudo/DR/DR_genes.window.bed | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_control_lambda.bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$8}'"

# window_treat_intensity_cmd <- "grep -e \"\\s1$\" /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/include_pseudo/VAR/VAR_genes.peak.window.50bp.bed | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_treat_pileup.bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$8}'"
# window_cont_intensity_cmd<-"grep -e \"\\s1$\" /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/include_pseudo/VAR/VAR_genes.peak.window.50bp.bed | /picb/evolgen/users/gushanshan/software/bedtools/bedtools intersect -a /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7_2rd_control_lambda.bdg -b stdin -wb | sort -k1,1 -k2,2n | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$8}'"
# 5. process -------------------------------------------------------------- TODO:
# print(window_treat_intensity_cmd)
# print(window_cont_intensity_cmd)
window_treat_intensity <- fread(cmd = window_treat_intensity_cmd, stringsAsFactors = F, header = F)
window_cont_intensity <- fread(cmd = window_cont_intensity_cmd, stringsAsFactors = F, header = F)
colnames(window_treat_intensity) <- colnames(window_cont_intensity) <- c("chrom", "sPos_minus_1", "ePos", "intensity", "gene")

treat.height.len <- window_treat_intensity$ePos - window_treat_intensity$sPos_minus_1
treat.height.value <- window_treat_intensity$intensity
treat.height <- rep(treat.height.value, treat.height.len)
treat.height_genome_position <- unlist(Map(`:`, window_treat_intensity$sPos_minus_1 + 1, window_treat_intensity$ePos))


cont.height.len <- window_cont_intensity$ePos - window_cont_intensity$sPos_minus_1
cont.height.value <- window_cont_intensity$intensity
cont.height <- rep(cont.height.value, cont.height.len)
cont.height_genome_position <- unlist(Map(`:`, window_cont_intensity$sPos_minus_1 + 1, window_cont_intensity$ePos))

treat_level_each_site <- treat.height
treat_pos_each_site <- treat.height_genome_position
control_level_each_site <- cont.height
control_pos_each_site <- cont.height_genome_position

if (all(control_pos_each_site == treat_pos_each_site)) {
    fc_genes <- sapply(unique(window_treat_intensity$gene), function(x) {
        window_gene_data <- window_treat_intensity %>% filter(gene == x)
        window_gene_min <- unlist(window_gene_data[1, 2] + 1)
        window_gene_max <- unlist(window_gene_data[nrow(window_gene_data), 3])
        window_gene_len <- sum(window_gene_data[, ePos] - window_gene_data[, sPos_minus_1])
        window_gene_index <- which(treat_pos_each_site >= window_gene_min & treat_pos_each_site <= window_gene_max)
        window_gene_fc_average <- mean(treat_level_each_site[window_gene_index] / control_level_each_site[window_gene_index])
    })
    # fc_window <- mean(fc_genes)
    # cat(fc_window)

    cat(fc_genes)
} else {
    cat("error")
}
