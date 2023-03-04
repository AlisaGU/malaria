#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(ggplot2)

# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/heatmap")

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
data <- read.table("Rhoptry_list_RNA-seq_V2.txt", header = T, sep = "\t")
data_R <- data[, 4:10]
data_T <- data[, 11:17]
data_S <- data[, 18:24]

data_R$Rstage_WT_TPM_AV <- apply(data_R[, 3:4], 1, mean)
data_R$Rstage_KD_TPM_AV <- apply(data_R[, 5:6], 1, mean)
data_R$Stage <- "R"

data_T$Tstage_WT_TPM_AV <- apply(data_T[, 3:4], 1, mean)
data_T$Tstage_KD_TPM_AV <- apply(data_T[, 5:6], 1, mean)
data_T$Stage <- "T"

data_S$Sstage_WT_TPM_AV <- apply(data_S[, 3:4], 1, mean)
data_S$Sstage_KD_TPM_AV <- apply(data_S[, 5:6], 1, mean)
data_S$Stage <- "S"

data_R <- data_R[, c(1, 2, 8, 9, 10)]
data_S <- data_S[, c(1, 2, 8, 9, 10)]
data_T <- data_T[, c(1, 2, 8, 9, 10)]
colnames(data_R) <- colnames(data_S) <- colnames(data_T) <- c("logfc", "pvalue", "wt_tpm", "kd_tpm", "stage")

data_for_plot <- data.frame(
    gene_symbol = rep(data$Symbol, times = 3),
    rbind(rbind(data_R, data_S), data_T)
)
data_for_plot$stage <- factor(data_for_plot$stage, levels = c("R", "T", "S"))
data_for_plot$pvalue <- -log(data_for_plot$pvalue, base = 10)
data_for_plot$direc <- ifelse(data_for_plot$logfc > 0, "Up", "Down")
data_for_plot$direc <- factor(data_for_plot$direc, levels = c("Down", "Up"))
data_for_plot$gene_symbol <- factor(data_for_plot$gene_symbol, levels = rev(data$Symbol))

data_for_plot$logfc[data_for_plot$logfc <= -3] <- -3
## FC+p value
ggplot(data_for_plot, aes(x = stage, y = gene_symbol, color = logfc, size = pvalue)) +
    geom_point() +
    scale_size(range = c(1, 5)) +
    # scale_shape_manual(values = c(25,24))+
    # scale_color_gradient2(low = "#8c510a",mid="white",midpoint=0, high = "#01665e",limits=c(-3,3)) +
    scale_color_gradientn(colors = rev(c(colorRampPalette(c("#01665e", "white"))(300), "white", colorRampPalette(c("white", "#eeddc9"))(100), colorRampPalette(c("#eeddc9", "#8c510a"))(200))), limits = c(-3, 3)) +
    labs(size = expression(paste("-", log[10], " p value"), parsed = T), color = expression(paste(log[2], " fold change"), parsed = T)) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12, color = "black")
    )
ggsave("logfc_pvalue.pdf", height = 10, width = 5)


data_for_plot$gene_symbol <- factor(data_for_plot$gene_symbol,
    levels = rev(levels(data_for_plot$gene_symbol))
)
ggplot(data_for_plot, aes(x = gene_symbol, y = stage, color = logfc, size = pvalue)) +
    geom_point() +
    scale_size(range = c(1, 5)) +
    # scale_shape_manual(values = c(25,24))+
    # scale_color_gradient2(low = "#8c510a",mid="white",midpoint=0, high = "#01665e",limits=c(-3,3)) +
    scale_color_gradientn(colors = rev(c(colorRampPalette(c("#01665e", "white"))(300), "white", colorRampPalette(c("white", "#eeddc9"))(100), colorRampPalette(c("#eeddc9", "#8c510a"))(200))), limits = c(-3, 3)) +
    labs(size = expression(paste("-", Log[10], "(", italic(P), " value", ")"), parsed = T), color = expression(paste(Log[2], "(fold change)"), parsed = T)) +
    theme_bw() +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15, color = "black"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12, color = "black")
    )
ggsave("logfc_pvalue_horizontal.pdf", height = 3, width = 12)
ggsave("logfc_pvalue_horizontal_ai.pdf", height = 3, width = 12)
