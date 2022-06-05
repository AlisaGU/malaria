#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------


# 2. functions ------------------------------------------------------------
plot_specific_variants <- function(variants = NULL, region = NULL) {
    data_for_plot <- lapply(variants, function(variant) {
        r0.1 <- c()
        r0.2 <- c()
        if (region == "whole") {
            r0.1 <- read.table(paste0("central_motif_relative.0.1_mutation_density/", variant, "_mean_variant_count"), header = F, as.is = T)
            r0.2 <- read.table(paste0("central_motif_relative.0.2_mutation_density/", variant, "_mean_variant_count"), header = F, as.is = T)
        } else if (region == "base") {
            r0.1 <- read.table(paste0("central_motif_relative.0.1_mutation_density_base_count/", variant, "_mean_variant_count"), header = F, as.is = T)
            r0.2 <- read.table(paste0("central_motif_relative.0.2_mutation_density_base_count/", variant, "_mean_variant_count"), header = F, as.is = T)
        }

        a <- data.frame(variant = rep(variant, 28), x = c(1:14, 1:14), value = c(r0.1[, 2], r0.2[, 2]), source = c(rep("10%", 14), rep("20%", 14)))
    })
    data_for_plot <- do.call(rbind, data_for_plot)


    p <- ggplot(data_for_plot, aes(x = variant, y = value, fill = factor(source))) +
        geom_boxplot()
    return(p)
}

# 3. input ----------------------------------------------------------------
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/OCE/motif/")

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------

plot_specific_variants(variants = c(
    "A2G", "G2A", "T2C", "C2T",
    "A2C", "A2T", "G2C", "G2T",
    "C2A", "C2G", "T2A", "T2G"
), region = "whole")
plot_specific_variants(variants = c("snps"), region = "whole")
plot_specific_variants(variants = c(
    "A2G", "G2A", "T2C", "C2T",
    "A2C", "A2T", "G2C", "G2T",
    "C2A", "C2G", "T2A", "T2G"
), region = "base")
