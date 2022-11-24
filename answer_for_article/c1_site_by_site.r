ggplot(data_for_plot, aes(x = site_f, y = value)) +
    geom_boxplot(position = position_dodge(), outlier.shape = NA, width = 0.8, size = 1) +
    geom_point(aes(color = chrom), size = 4) +
    scale_color_manual(values = c(
        "chrom1" = rgb(5, 48, 97, maxColorValue = 255),
        "chrom2" = rgb(66, 147, 195, maxColorValue = 255),
        "chrom3" = rgb(116, 42, 131, maxColorValue = 255),
        "chrom4" = rgb(251, 77, 41, maxColorValue = 255),
        "chrom5" = rgb(0, 69, 27, maxColorValue = 255),
        "chrom6" = rgb(180, 24, 45, maxColorValue = 255),
        "chrom10" = rgb(30, 33, 223, maxColorValue = 255),
        "chrom8" = rgb(254, 204, 59, maxColorValue = 255),
        "chrom9" = rgb(166, 112, 0, maxColorValue = 255),
        "chrom7" = rgb(84, 48, 4, maxColorValue = 255),
        "chrom11" = rgb(245, 150, 94, maxColorValue = 255),
        "chrom12" = rgb(45, 125, 180, maxColorValue = 255),
        "chrom13" = rgb(227, 66, 47, maxColorValue = 255),
        "chrom14" = rgb(168, 132, 189, maxColorValue = 255)
    )) +
    labs(y = "Relative mutation density") +
    scale_x_discrete(labels = c(paste("S", 1:6, sep = ""), paste("F", 1:10, sep = ""))) +
    theme_bw() +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 44, color = "black", face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(
            size = 36, color = "black",
        ),
        axis.text.y = element_text(size = 36, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, color = "black"),
        legend.position = "bottom",
        plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"),
        panel.spacing = unit(3, "lines")
    )



f1 <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/flank_nomatter_strand.remove_close_overlap/flank_1/flank1.peak.base_count", header = F, as.is = T)

f3 <- read.table("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/flank_nomatter_strand.remove_close_overlap/flank_3/flank3.peak.base_count", header = F, as.is = T)

colnames(f1) <- colnames(f3) <- c("A", "T", "C", "G", "all")
rownames(f1) <- rownames(f3) <- paste("chrom", 1:14, sep = "")

## sort by A
f1 <- as.data.frame(t(apply(f1[, 1:4], 1, function(x) {
    x / sum(x)
})))

f1_1 <- data.frame(value = unlist(f1), base = rep(colnames(f1), each = 14), chrom = rep(rownames(f1), times = ncol(f1)))
f1_1$base <- factor(f1_1$base, levels = c("A", "T", "C", "G"))
f1_1 <- f1_1 %>% arrange(base, value)
f1_1$chrom <- factor(f1_1$chrom, levels = f1_1$chrom[1:14])

ggplot(f1_1, aes(x = chrom, y = value, fill = base)) +
    geom_col(position = "stack", width = 0.6)



f3 <- as.data.frame(t(apply(f3[, 1:4], 1, function(x) {
    x / sum(x)
})))

f3_1 <- data.frame(value = unlist(f3), base = rep(colnames(f3), each = 14), chrom = rep(rownames(f3), times = ncol(f3)))
f3_1$base <- factor(f3_1$base, levels = c("A", "T", "C", "G"))
f3_1 <- f3_1 %>% arrange(base, value)
f3_1$chrom <- factor(f3_1$chrom, levels = f3_1$chrom[1:14])

ggplot(f3_1, aes(x = chrom, y = value, fill = base)) +
    geom_col(position = "stack", width = 0.6)

## sort by A+T
f1 <- as.data.frame(t(apply(f1[, 1:4], 1, function(x) {
    sum(x[c(1, 2)]) / sum(x)
})))

data <- data.frame(value = unlist(f1), chrom = names(f1))
data <- data[order(data$value), ]
data$chrom <- factor(data$chrom, levels = data$chrom)
ggplot(data, aes(x = chrom, y = value)) +
    geom_bar(stat = "identity", fill = "blue")


f3 <- as.data.frame(t(apply(f3[, 1:4], 1, function(x) {
    sum(x[c(1, 2)]) / sum(x)
})))

data <- data.frame(value = unlist(f3), chrom = names(f3))
data <- data[order(data$value), ]
data$chrom <- factor(data$chrom, levels = data$chrom)
ggplot(data, aes(x = chrom, y = value)) +
    geom_bar(stat = "identity", fill = "blue")
