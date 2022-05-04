#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(data.table)


# 2. functions ------------------------------------------------------------
get_inconsistent_region <- function(inputFile = NULL) {
    GT_distant <- fread(inputFile, sep = "\t", stringsAsFactors = F, header = F, select = 2048)
    # posi <- fread(cmd = paste0(bcftools, " view -H ", inputFile, ".vcf.gz|awk '{print $1,$2}'"), stringsAsFactors = F)
    posi <- fread(cmd = paste0(bcftools, " view -H ", inputFile, ".vcf.gz"), select = c(1, 2), stringsAsFactors = F)

    inconsistent_region <- posi[which(GT_distant != 0), ]
    return(inconsistent_region)
}

write_pos <- function(inconsistent_region = NULL, inconsistent_region_file = NULL) {
    if (file.exists(inconsistent_region_file)) {
        file.remove(inconsistent_region_file)
    }

    for (i in 1:length(inconsistent_region)) {
        a <- inconsistent_region[[i]]
        b <- data.frame(a[, 1], a[, 2] - 1, a[, 2])
        write.table(format(b, scientific = FALSE), inconsistent_region_file,
            quote = F, sep = "\t", row.names = F, col.names = F, append = T
        )
    }
}
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
inconsistent_region_file <- Args[6]

# setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/OCE/africa_biallelic_variant")
bcftools <- "/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
# 4. variable setting of test module---------------------------------------
# inconsistent_region_file<-"/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core/inconsistent_region_in_core.bed"

# 5. process --------------------------------------------------------------
chroms <- c(
    "01", "02", "03", "04", "05", "06", "07",
    "08", "09", "10", "11", "12", "13", "14"
)
vcf_files <- paste("WAF_CAF_EAF_biallelic_", chroms, sep = "")
inconsistent_region <- list()
for (i in 1:14) {
    print(i)
    inputFile <- vcf_files[[i]]
    inconsistent_region[[i]] <- get_inconsistent_region(inputFile = inputFile)
}
names(inconsistent_region) <- vcf_files
