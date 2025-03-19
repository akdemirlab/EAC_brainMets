
library(readxl)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)  
library(paletteer)
library(ggpubr)
library(ggsignif)

Pt.primary <- "/Primary_cell_types.tsv"
Pt.brainmet <- "/brainmet_cell_types.tsv"


Pt.primary.data <- read.table(Pt.primary, header = TRUE, sep = "\t")
Pt.data.brainmet <- read.table(Pt.brainmet, header = TRUE, sep = "\t")



Pt.primary.data$Sample <-"Pt.primary"
Pt.data.brainmet$Sample <-"Pt.brainmet"


Pt.primary_filter <- Pt.primary.data %>%
  select(Sample, Group) %>%
  group_by(Sample, Group) %>%
  summarize(Count = n())

Pt.data.brainmet_filter <- Pt.data.brainmet  %>%
  select(Sample, Group) %>%
  group_by(Sample, Group) %>%
  summarize(Count = n())


# Combining datasets and converting counts for metastasis to negative
combined_data <- rbind(
  cbind(Pt.primary_filter, Type = 'Primary'),
  cbind(Pt.data.brainmet_filter, Type = 'Brain Met')
)

combined_data$Count[combined_data$Type == 'Brain Met'] <- -combined_data$Count[combined_data$Type == 'Brain Met']

# Calculate total counts
total_primary <- sum(combined_data$Count[combined_data$Type == 'Primary'])
total_metastasis <- sum(abs(combined_data$Count[combined_data$Type == 'Brain Met']))

# Calculate percentages
combined_data$Percentage <- ifelse(combined_data$Type == 'Primary', 
                                   (combined_data$Count / total_primary) * 100, 
                                   ((combined_data$Count - 0) / total_metastasis) * 100)  # Subtract metastasis counts from 0

celltype_order <- c("Fibroblast", "T-cell", "B-cell", "Plasma cell", "Macrophage", "Potential Neuronal Cells", "Malignant")
combined_data <- mutate(combined_data, Group = factor(Group, levels = celltype_order))

combined_data <- na.omit(combined_data)


# Plotting
ggplot(combined_data, aes(x = Group, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Types", y = "Percentage of Cells", title = "Percentage of Cells per Cell Type Pt (Primary vs Brain Met)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Fibroblast" = "#21EFA5", "T-cell" = "#9F70AA", "B-cell" = "#436ae7", "Dendritic cell" = "#87EEF0", "Plasma cell" = "#D9D9D9", "Mast cell" = "#969696", "Macrophage" = "#A6912E",  "Malignant" = "#FF6464")) +
  theme_minimal() +
  scale_y_continuous(limits = c(-100, 100))  # Set limits of y-axis from -10



