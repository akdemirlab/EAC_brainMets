# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(hrbrthemes)
library(viridis)
library(paletteer)
library(scales)  

# Define the input file names
input_file1 <- "ecDNA_counts_BarrettsEso_Luebeck2023.xlsx"
input_file2 <- "ecDNA_counts_mydata.xlsx"
input_file3 <- "histology_ecDNA_counts_Luebeck.xlsx"

# Read the input Excel files
data1_ecdna <- read_excel(input_file1)
data2_ecdna <- read_excel(input_file2)
histology_data <- read_excel(input_file3)


# Set Sample names
data2_ecdna$Sample <- "EAC Brain Mets"


# Summarize ecDNA counts by 'sample_id'
data1_ecdna <- data1_ecdna %>%
  group_by(sample_id) %>%
  summarize(ecDNA = sum(ecDNA, na.rm = TRUE))

# Merge the data frames based on the sample_id column
merged_data <- merge(data1_ecdna, histology_data, by = "sample_id", all.x = TRUE)

# The `all.x = TRUE` argument keeps all rows from ecDNA_data and adds histology from histology_data

names(merged_data)[names(merged_data) == "Histology"] <- "Sample"

# Convert the sample_id column to character in data2_ecdna
data2_ecdna$sample_id <- as.character(data2_ecdna$sample_id)


# Combine the dataframes
combined_data_ecDNA <- bind_rows(merged_data, data2_ecdna)

combined_data_ecDNA$Sample <- factor(combined_data_ecDNA$Sample, levels = c("BE", "EAC", "EAC Brain Mets"))
library(ggplot2)
library(dplyr)

filtered_combined_data_ecDNA <- combined_data_ecDNA[!is.na(combined_data_ecDNA$Sample), ]


# Calculate the percentage of samples with at least one ecDNA for each group
percentage_data <- filtered_combined_data_ecDNA %>%
  group_by(Sample) %>%
  summarize(True = sum(ecDNA > 0) / n() * 100) %>%
  mutate(False = 100 - True)    %>%
  pivot_longer(cols = c(True, False), names_to = "PercentageType", values_to = "Value")


# Create a stacked bar plot
ggplot(percentage_data, aes(x = Sample, y = Value, fill = PercentageType)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of Samples with At Least 1 ecDNA", x = "Sample Group", y = "Percentage") +
  scale_fill_manual(values = c("True" = "#66c2a5", "False" = "#fc8d62")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100)
  )
