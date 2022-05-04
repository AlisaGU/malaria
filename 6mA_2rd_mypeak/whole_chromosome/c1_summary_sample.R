#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:
library(data.table)
library(rworldmap)
library(maps)
library(ggplot2)
library(patchwork)
# 2. functions ------------------------------------------------------------ TODO:
get_unique <- function(data = NULL) {
    countrys <- unique(data$Country)
    result <- c()
    for (i in 1:length(countrys)) {
        country.i <- countrys[i]
        corordinates <- unlist(data[match(country.i, data$Country), 2:4])
        result <- rbind(result, c(length(which(data$Country == country.i)), corordinates))
    }
    colnames(result) <- c("Sample_count", "Long", "Lat", "Population")
    result <- as.data.frame(result)
    result$Sample_count <- as.numeric(result$Sample_count)
    result$Long <- as.numeric(result$Long)
    result$Lat <- as.numeric(result$Lat)
    result$country <- countrys
    return(result)
}

# 3. input ---------------------------------------------------------------- TODO:
setwd("/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd")
sample_metadata <- fread("sample_provenance_and_sequencing_metadata.txt", header = T, stringsAsFactors = F)


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
QC_pass_sample_data <- sample_metadata[sample_metadata$QCpass, ]
data_for_map <- get_unique(data = QC_pass_sample_data[, c("Country", "Long", "Lat", "Population")])
country_distri <- ggplot(data_for_map, aes(x = Long, y = Lat)) +
    borders("world", colour = NA, fill = "wheat1") +
    geom_point(aes(size = Sample_count, color = Population), alpha = 0.5) +
    geom_text(aes(label = country), color = "gray20", fontface = "italic", check_overlap = F, size = 1.5) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
    )
sample_count_distri <- ggplot(data_for_map, aes(x = Long, y = Lat)) +
    borders("world", colour = NA, fill = "wheat1") +
    geom_point(aes(size = Sample_count, color = Population), alpha = 0.5) +
    geom_text(aes(label = Sample_count), color = "gray20", fontface = "italic", check_overlap = F, size = 1.5) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
    )
# (country_distri | sample_count_distri) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# ggsave(filename = "sample_distribution.pdf", width = 30, height = 8)
pdf("sample_distribution.pdf", width = 15, height = 8)
country_distri
sample_count_distri
dev.off()