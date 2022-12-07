#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/motif/mutation_dentisy_two_outgroup_consistent/each_site_all_motifs_asWhole_LI_formula")
motifs <- c("GAAGAA", "GAAGAT", "GATGAA", "GATGAT")
sites <- c("site1", "site2", "site3", "site4", "site5", "site6")
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- sapply(motifs, function(motif) {
    site_mutation_average_count <- sapply(sites, function(site) {
        mutation_count <- read.table(paste0(motif, "/", site, ".snps_mean_variant_count"))
        mutation_count <- mutation_count[, c(1:8)]
        colnames(mutation_count) <- c(
            "A_mutation_count", "A_base_total_count",
            "T_mutation_count", "T_base_total_count",
            "C_mutation_count", "C_base_total_count",
            "G_mutation_count", "G_base_total_count"
        )

        just_mutation_count <- apply(mutation_count[, c(1, 3, 5, 7)], 1, sum)
        base_count <- apply(mutation_count[, c(2, 4, 6, 8)], 1, sum)
        return(sum(just_mutation_count))
    })
})
apply(data, 2, mean)
