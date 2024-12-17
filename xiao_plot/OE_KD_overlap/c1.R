#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)
library(dplyr)
# 2. functions ------------------------------------------------------------ TODO:


# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/OE_KD_overlap")

# 5. process -------------------------------------------------------------- TODO:
data <- read.table("KD_OE_overlap_GO.txt", header = T, stringsAsFactors = F, as.is = T, sep = "\t")
data$Percentage <- data$Result.count / data$Bgd.count * 100
data <- data[order(data$X.log10.p, decreasing = F), ]
data$Name <- factor(data$Name, levels = data$Name)
colnames(data)[colnames(data) == "Result.count"] <- "Count"

ggplot(data, aes(x = Name, y = X.log10.p)) +
    geom_point(aes(size = Count, color = Percentage)) +
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

ggsave("OE KD overlap GO.pdf", width = 8, height = 5)
