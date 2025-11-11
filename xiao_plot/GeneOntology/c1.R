#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:


# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/GeneOntology/")

# 5. process -------------------------------------------------------------- TODO:
filelist <- dir()[grep("tsv", dir())]
for (filename in filelist) {
    data <- read.table(filename, header = T, stringsAsFactors = F, sep = "\t")
    data <- data[, c("Name", "P.value", "Result.count", "Pct.of.bgd")]
    data$sig <- -log10(data$P.value)

    data <- data[order(data$sig, decreasing = F), ]
    data$Name <- factor(data$Name, levels = data$Name)

    ggplot(data, aes(x = Name, y = sig)) +
        geom_point(aes(size = Result.count, color = Pct.of.bgd)) +
        scale_color_gradient(low = "#fc000a", high = "#2900fa") +
        labs(y = expression(-log[10](italic(P)))) +
        coord_flip() +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            plot.title = element_text(
                colour = "black", face = "bold",
                size = 14, vjust = 1, hjust = 0.5
            ),
            axis.title.x = element_text(size = 20, color = "black"),
            axis.title.y = element_blank()
        )

    ggsave(gsub("tsv", "pdf", filename), width = 8, height = 5)
}
