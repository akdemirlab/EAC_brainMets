
library(ggplot2)
library(dplyr)
library(tidyr)  
library(readxl)


Brainmet_ascat_CN$New_LogCN <- log(Brainmet_ascat_CN$Copy_Number, base = 2)
TCGA_ascat_CN$New_LogCN  <- log(TCGA_ascat_CN$Copy_Number, base = 2)



# Filter for specific genes
library(dplyr)

Brainmet_ascat_CN_specificgenes <- Brainmet_ascat_CN %>%
  filter(Gene %in% c("ERBB2", "KRAS", "MYC", "EGFR")) %>%
  dplyr::select(-LogCN)


# Filter for specific genes
TCGA_ascat_CN_specificgenes <- TCGA_ascat_CN %>%
  filter(Gene %in% c("ERBB2", "KRAS", "MYC", "EGFR")) %>%
  dplyr::select(-LogCN)


# Define the desired order of genes
desired_gene_order <- c("ERBB2", "MYC", "EGFR", "KRAS")  

# Set the Gene variable to be a factor with the desired order
Brainmet_ascat_CN_specificgenes$Gene <- factor(Brainmet_ascat_CN_specificgenes$Gene, levels = desired_gene_order)
TCGA_ascat_CN_specificgenes$Gene <- factor(TCGA_ascat_CN_specificgenes$Gene, levels = desired_gene_order)

Brainmet_ascat_CN_specificgenes$Sample <- "EAC Brain Mets"
TCGA_ascat_CN_specificgenes$Sample <- "TCGA Primary EAC"

library(dplyr)

# Specify the column name you want to remove
column_to_remove <- "Sample ID"

# Remove the specified column
Brainmet_ascat_CN_specificgenes <- Brainmet_ascat_CN_specificgenes %>%
  select(-{{column_to_remove}})

TCGA_ascat_CN_specificgenes <- TCGA_ascat_CN_specificgenes %>%
  select(-{{column_to_remove}})




combined_data_CN_TCGA_brainmet <- bind_rows(Brainmet_ascat_CN_specificgenes, TCGA_ascat_CN_specificgenes)


combined_data_CN_TCGA_brainmet$Sample <- factor(combined_data_CN_TCGA_brainmet$Sample, levels = c("TCGA Primary EAC", "EAC Brain Mets"))
# Calculate sample counts for each gene
sample_counts6 <- combined_data_CN_TCGA_brainmet %>%
  group_by(Gene, Sample) %>%
  summarize(SampleCount = n())

ggplot_CN_TCGA_brainmet <- ggplot(combined_data_CN_TCGA_brainmet, aes(x = Gene, y = New_LogCN, fill = Sample)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), aes(color = Sample), alpha = 0.5, size = 1) +
  labs(title = "Log2 Copy Number of Oncogenes",
       y = "Log2 Copy Number",
       x = "Genes"
  ) +
  scale_fill_manual(values = c("TCGA Primary EAC" = "#aec7e8", "EAC Brain Mets" = "#98df8a")) +  
  scale_color_manual(values = c("TCGA Primary EAC" = "#aec7e8", "EAC Brain Mets" = "#98df8a")) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "right"  # Adjust legend position
  ) +
  facet_grid(. ~ Gene, scales = "free", space = "free") 

ggplot_CN_TCGA_brainmet <- ggplot_CN_TCGA_brainmet + geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
  geom_text(data = sample_counts6, aes(x = Gene, y = 8, label = paste("n=", SampleCount)), position = position_dodge(width = 0.75), size = 3)

print(ggplot_CN_TCGA_brainmet)

library(ggpubr)

# Add significance asterisks using ggpubr
ggplot_CN_TCGA_brainmet_sig <- ggplot_CN_TCGA_brainmet + stat_compare_means(
  method = "wilcox.test",
  method.args = list(exact = FALSE),
  label = "p.signif",
  tip.length = 0.5,
  label.y = c(8.3)
)


# Print the plot
print(ggplot_CN_TCGA_brainmet_sig)
