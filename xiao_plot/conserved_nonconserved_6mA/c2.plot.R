#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(dplyr)
library(ggplot2)
# 2. functions ------------------------------------------------------------ TODO:
readData <- function() {
    filenames <- dir()[grep("chip|input", dir(), perl = T)]
    result <- lapply(filenames, function(x) {
        genomicRegion <- unlist(strsplit(x, ".", fixed = T))[1]
        group <- unlist(strsplit(x, ".", fixed = T))[3]
        data <- read.table(x, header = F, as.is = T)
        colnames(data) <- c("chrom", "start", "end", "readCount")
        data$len <- data$end - data$start
        result <- data %>%
            group_by(chrom) %>%
            mutate(totalReadCount = sum(readCount)) %>%
            select(chrom, totalReadCount) %>%
            distinct()
        result$genomicRegion <- genomicRegion
        result$group <- group
        return(result)
    })
    result <- do.call(rbind, result)

    finalResult <- lapply(c("conservedGenomicRegion", "variableGenomicRegion"), function(x) {
        chipData <- result %>% filter(genomicRegion == x & group == "chip")
        inputData <- result %>% filter(genomicRegion == x & group == "input")
        integratedData <- dplyr::left_join(chipData, inputData, by = c("chrom" = "chrom"), suffix = c(".chip", ".input")) %>% select(chrom, totalReadCount.chip, totalReadCount.input)
        integratedData$normalizedReadCount <- (integratedData$totalReadCount.chip / chipLibrarySize) / (integratedData$totalReadCount.input / inputLibrarySize)
        integratedData$genomicRegion <- x
        return(integratedData)
    })
    finalResult <- do.call(rbind, finalResult)
    return(finalResult)
}

# 3. variable setting of test module--------------------------------------- TODO:


# 4. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/conserved_nonconserved_6mA")
chipLibrarySize <- 36971023
inputLibrarySize <- 21500156
# 5. process -------------------------------------------------------------- TODO:
data_for_plot <- readData()
data_for_plot <- data_for_plot %>% filter(!chrom %in% c("Pf3D7_API_v3", "Pf_M76611"))

t.test(data_for_plot$normalizedReadCount[data_for_plot$genomicRegion == "conservedGenomicRegion"], data_for_plot$normalizedReadCount[data_for_plot$genomicRegion == "variableGenomicRegion"], paired = T)

ggplot(data_for_plot, aes(x = genomicRegion, y = normalizedReadCount)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.1)) +
    scale_x_discrete(labels = c("Conserved\ngenome", "Variable\ngenome")) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "grey"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(colour = "grey"),
        plot.title = element_text(
            colour = "black", face = "bold",
            size = 14, vjust = 1, hjust = 0.5
        ),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
ggsave("conserved_nonconserved_6mA.pdf", height = 5, width = 5)
