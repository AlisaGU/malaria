#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ----------------------------------------
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(grid)
library(patchwork)
library(ggrepel)

# 2. functions ------------------------------------------------------------
read_mutation_density <- function() {
    mutations <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
    data <- lapply(mutations, function(mutation) {
        a <- read.table(paste0(mutation, "_mean_variant_count"), header = F, as.is = T)
        rownames(a) <- c(
            "Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3",
            "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3", "Pf3D7_09_v3",
            "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3"
        )
        colnames(a) <- c(
            "nopeak", "r001", "r005",
            "r010", "r020", "r050", "a010",
            "a030", "a050", "a100", "a200"
        )
        return(a)
    })
    names(data) <- mutations
    return(data)
}

get_mutation_density_for_regions <- function(regions = NULL, mutation_load = NULL) {
    result <- lapply(regions, function(region) {
        a <- lapply(mutation_load, function(x) {
            return(x[, region])
        })
        data <- lapply(1:length(a), function(i) {
            cbind(a[[i]], names(a)[i])
        })
        data <- do.call(rbind, data)
        data <- data.frame(value = as.numeric(data[, 1]), mutation = data[, 2])
    })
    names(result) <- regions
    return(result)
}

my_plot <- function(data = NULL, title = NULL) {
    p <- ggboxplot(data = data, x = "mutation", y = "value", color = "mutation", add = "point", add.params = list(size = 3)) +
        labs(title = title) +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            # axis.text.x = element_text(size = 24, color = "black", angle = 30, hjust = 0.7, vjust = 0.7),
            axis.text.y = element_text(size = 24, color = "black"),
            plot.title = element_text(
                colour = "black",
                size = 28, vjust = 0.5, hjust = 0.5
            ),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 24, color = "black"),
            axis.title.y = element_text(size = 28, color = "black"),
            legend.position = "None",
            plot.margin = margin(0.2, 0.2, 0.3, 0.2, "cm")
        )
    # compare_means(value ~ mutation, data = data)
    return(p)
}
# 3. input ----------------------------------------------------------------
Args <- commandArgs()
popu_symbol <- Args[6]
mutations <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
regions <- c("nopeak", "r001", "r005", "r010", "r020", "r050", "a010", "a030", "a050", "a100", "a200")
setwd(paste0("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic/", popu_symbol, "/mutation_density"))

# 4. variable setting of test module---------------------------------------


# 5. process --------------------------------------------------------------
mutation_load <- read_mutation_density()
mutation_load_for_regions <- get_mutation_density_for_regions(regions = regions, mutation_load = mutation_load)
my_plot(data = mutation_load_for_regions$nopeak, title = "nopeak")

ggsave("/picb/evolgen2/users/gushanshan/projects/malaria/code/nopeak.pdf",width=10,height=10)