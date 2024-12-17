#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
suppressWarnings(suppressMessages(library(preprocessCore)))
library(ggplot2)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
detect_operating_system <- function() {
    switch(Sys.info()[["sysname"]],
        Windows = {
            return("Windows")
        },
        Linux = {
            return("Linux")
        }
    )
}

path_convert <- function(path = NULL) {
    result <- paste0("I:", gsub("\\picb\\evolgen\\users\\gushanshan", "", gsub("/", "\\", path, fixed = T), fixed = T))
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
server_path <- "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_SachelMok"
GSEid <- "GSE59099"
data_repository <- paste0(GSEid, ".exp_meta_probe.Rdata")
probes <- c("aMAL13P1.255_0", "aMAL13P1.255_1") # 在GPL18893平台上，PF3D7_1350700基因对应两个探针：aMAL13P1.255_0和aMAL13P1.255_1
argine_probes <- c("aPFI0680c_1", "aPFI0680c_0") # argine的基因id为PF3D7_0913900
# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
setwd(ifelse(detect_operating_system() == "Windows", path_convert(server_path), server_path))
load(data_repository)
metainfo <- GSEdata[[1]]$metaInfo
sample_clearanceTime <- metainfo[, c("geo_accession", "characteristics_ch1.6")]
colnames(sample_clearanceTime) <- c("sample", "clearanceTime")

## 排除体外样本
sample_clearanceTime <- sample_clearanceTime[grep("parasite clearance halflife upon artemisinin treatment (h): ",
    sample_clearanceTime$clearanceTime,
    fixed = T
), ]

a <- gsub("parasite clearance halflife upon artemisinin treatment (h): ", "", sample_clearanceTime$clearanceTime, fixed = T)
sample_clearanceTime$clearanceTime <- gsub("h", "", a, fixed = T)
na_index <- which(sample_clearanceTime$clearanceTime == "NA") # 18个样本的清除时间为NA，排除这18个样本
sample_clearanceTime <- sample_clearanceTime[-na_index, ]
sample_clearanceTime$clearanceTime <- as.numeric(sample_clearanceTime$clearanceTime)
sample_clearanceTime$reponse <- sapply(sample_clearanceTime$clearanceTime, function(x) {
    ifelse(x <= 5, "sensitive", "resistant")
})

# 1. 恢复为原始数据
exp <- GSEdata[[1]]$exp
newExp <- 2^(exp)

normExp <- newExp
# 2. 按照XB的操作重新做一遍


namt_normExp <- normExp[match(probes, rownames(normExp)), match(sample_clearanceTime$sample, colnames(normExp))]
namt_normExp_means <- colMeans(namt_normExp, na.rm = TRUE)

argine_normExp <- normExp[match(argine_probes, rownames(normExp)), match(sample_clearanceTime$sample, colnames(normExp))]
argine_normExp_means <- colMeans(argine_normExp, na.rm = TRUE)


data_for_plot <- data.frame(t(rbind(namt_normExp, namt_normExp_means, argine_normExp, argine_normExp_means)))
data_for_plot$response <- sample_clearanceTime$reponse[match(rownames(data_for_plot), sample_clearanceTime$sample)]
data_for_plot$response <- factor(data_for_plot$response, levels = c("resistant", "sensitive"))

data_for_plot$ratio <- data_for_plot$namt_normExp_means / data_for_plot$argine_normExp_means

namt <- ggplot(data_for_plot, aes(x = response, y = namt_normExp_means)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    labs(y = "namt (Cy5/Cy3)") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
namt_p_t <- t.test(data_for_plot$namt_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$namt_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value
namt_p_wilcox <- wilcox.test(data_for_plot$namt_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$namt_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value


argine <- ggplot(data_for_plot, aes(x = response, y = argine_normExp_means)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    labs(y = "Argine (Cy5/Cy3)") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
argine_p_t <- t.test(data_for_plot$argine_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$argine_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value
argine_p_wilcox <- wilcox.test(data_for_plot$argine_normExp_means[data_for_plot$response == "resistant"],
    data_for_plot$argine_normExp_means[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value



ratio <- ggplot(data_for_plot, aes(x = response, y = ratio)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun.y = median, geom = "point", size = 2, color = "red") +
    labs(y = "Fold Change\n(namt/argine)") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
ratio_p_t <- t.test(data_for_plot$ratio[data_for_plot$response == "resistant"],
    data_for_plot$ratio[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value
ratio_p_wilcox <- wilcox.test(data_for_plot$ratio[data_for_plot$response == "resistant"],
    data_for_plot$ratio[data_for_plot$response == "sensitive"],
    alternative = "greater"
)$p.value

namt | argine | ratio

p_value <- matrix(c(namt_p_t, argine_p_t, ratio_p_t, namt_p_wilcox, argine_p_wilcox, ratio_p_wilcox), byrow = T, nrow = 2)
rownames(p_value) <- c("t.test", "wilcox.test")
colnames(p_value) <- c("namt", "argine", "namt/argine")
write.table(p_value, "p.value.txt", sep = "\t", quote = F)


## ratio plot for publication
ggplot(data_for_plot, aes(x = response, y = ratio, fill = response)) +
    geom_violin(color = "black", size = 1.5, width = 0.4) +
    geom_boxplot(width = 0.2, fill = "white", size = 2, color = "black") +
    scale_fill_manual(values = c("resistant" = "#99c192", "sensitive" = "#ca7b42")) +
    labs(y = "Fold Change\n(namt/argine)") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 2.2)) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36, color = "black"),
        legend.position = "none"
    )
ggsave("/picb/evolgen2/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_SachelMok/2015_namt_argine.pdf", width = 6, height = 7)


# 数字统计
ratio_sensitive <- data_for_plot$ratio[data_for_plot$response == "sensitive"]
median(ratio_sensitive[!is.na(ratio_sensitive)])
length(ratio_sensitive[!is.na(ratio_sensitive)])


ratio_resistant <- data_for_plot$ratio[data_for_plot$response == "resistant"]
median(ratio_resistant[!is.na(ratio_resistant)])
length(ratio_resistant[!is.na(ratio_resistant)])
