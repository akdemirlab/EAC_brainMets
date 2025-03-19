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
input_file1 <- "Hartwig_EACmets_specific_drivers.xlsx"
input_file2 <- "mydata_specific_drivers.xlsx"
input_file3 <- "PCAWG_primaryEAC_specific drivers.xlsx"

# Read the input Excel files
data1_drgene <- read_excel(input_file1)
data2_drgene <- read_excel(input_file2)
data3_drgene <- read_excel(input_file3)

# Set Sample names
data1_drgene$Sample <- "Hartwig EAC Mets"
data2_drgene$Sample <- "EAC Brain Mets"
data3_drgene$Sample <- "PCAWG Primary EAC"

data1_drgene_filter <- data1_drgene %>%
  select(Sample, gene, mut, total)

data2_drgene_filter <- data2_drgene %>%
  select(Sample, gene, mut, total)

data3_drgene_filter <- data3_drgene %>%
  select(Sample, gene, mut, total)

# Combine the dataframes
combined_data_drgene <- bind_rows(data1_drgene, data2_drgene, data3_drgene)

combined_data_drgene$Sample <- factor(combined_data_drgene$Sample, levels = c("PCAWG Primary EAC", "Hartwig EAC Mets", "EAC Brain Mets"))

# Select only the relevant columns
genes_of_interest <- c("TP53", "CDKN2A", "KRAS", "MYC", "ERBB2", "CCND1", "FHIT")

combined_data_drgene_filter <- combined_data_drgene %>%
  filter(gene %in% genes_of_interest) %>%
  select(Sample, driver, gene, mut, total)

# Calculate the percentage of samples with driver genes for each gene and sample group
driver_percentage_data <- combined_data_drgene_filter %>%
  group_by(Sample, gene, driver) %>%
  summarize(
    mut = sum(mut),
    total = first(total)
  ) %>%
  ungroup() %>%
  group_by(Sample, gene) %>%
  mutate(
    Driver = ifelse(driver %in% c("AMP", "DEL", "MUTATION"), mut / total * 100, 0),
    NoDriver = ifelse(driver == "No driver", (total - mut) / total * 100, 0)
  ) %>%
  pivot_longer(cols = c(Driver, NoDriver), names_to = "PercentageType", values_to = "Value")

# Define the order of gene levels as desired
gene_order <- c("TP53", "ERBB2", "CDKN2A", "KRAS", "MYC", "CCND1", "FHIT")

# Set gene as a factor with the specified order
driver_percentage_data$gene <- factor(driver_percentage_data$gene, levels = gene_order)

ggplot(driver_percentage_data, aes(x = Sample, y = Value, fill = driver)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of Samples with Specific Driver Genes", x = "Sample Group", y = "Percentage") +
  scale_fill_manual(
    values = c("AMP" = "#66c2a5", "DEL" = "#fc8d62", "MUTATION" = "#8da0cb", "No driver" = "pink")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 14),
    text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100)  
  ) +
  guides(fill = guide_legend(title = "Driver Type")) +
  facet_grid(. ~ gene)
