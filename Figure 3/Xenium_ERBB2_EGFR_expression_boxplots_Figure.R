
library(readxl)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)  
library(paletteer)
library(ggpubr)
library(ggsignif)


Pt.primary <- "primary_transcriptcount_bycelltype.tsv"
Pt.brainmet <- "brainmet_transcriptcount_bycelltype.tsv"


# Read the TSV file into a data frame

Pt.primary.data <- read.table(Pt.primary, header = TRUE, sep = "\t")
Pt.data.brainmet <- read.table(Pt.brainmet, header = TRUE, sep = "\t")



Pt.primary.data$Sample <-"Pt.primary"
Pt.data.brainmet$Sample <-"Pt.brainmet"



Pt.primary_filter <- Pt.primary.data %>%
  select(Sample, ERBB2_counts, EGFR_counts, ERBB3_counts, Group) %>%
  mutate(
    ERBB2_log2 = log2(ERBB2_counts + 1),  
    EGFR_log2 = log2(EGFR_counts + 1),
    ERBB3_log2 = log2(ERBB3_counts + 1)
  )


Pt_filter <- Pt.data.brainmet %>%
  select(Sample, ERBB2_counts, EGFR_counts, ERBB3_counts, Group) %>%
  mutate(
    ERBB2_log2 = log2(ERBB2_counts + 1), 
    EGFR_log2 = log2(EGFR_counts + 1),
    ERBB3_log2 = log2(ERBB3_counts + 1)
  )



# Add an identifier column to each dataset
filtered_datap <- Pt.primary.data %>%
  filter(Group == "Malignant") %>%
  mutate(
    ERBB2_log2 = log2(ERBB2_counts + 1),  # Adding 1 to avoid log2(0)
    EGFR_log2 = log2(EGFR_counts + 1),
    ERBB3_log2 = log2(ERBB3_counts + 1)
  )


filtered_data <- Pt.data.brainmet %>%
  filter(Group == "Malignant") %>%
  mutate(
    ERBB2_log2 = log2(ERBB2_counts + 1),  # Adding 1 to avoid log2(0)
    EGFR_log2 = log2(EGFR_counts + 1),
    ERBB3_log2 = log2(ERBB3_counts + 1)
  )


# Combine the datasets
combined_data <- bind_rows(filtered_datap, filtered_data)

sample_order <- c("Patient.primary", "Patient.brainmet")
combined_data <- mutate(combined_data, Sample = factor(Sample, levels = sample_order))



PTERBB2_box <- ggplot(combined_data, aes(x = Sample, y = ERBB2_log2, fill = Sample)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), aes(color = Sample), alpha = 0.05, size = .5) +
  labs(x = "Sample", y = "ERBB2 log2+1 transcript counts", title = "ERBB2 expression Patient") +
  scale_fill_manual(values = c("Patient.primary" = "#C86050", "Patient.brainmet" = "#6685c2")) +
  scale_color_manual(values = c("Patient.primary" = "#C86050", "Patient.brainmet" = "#6685c2")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "transparent"),  
    panel.grid.minor = element_blank(),  
    plot.title = element_text(size = 14),  
    axis.text = element_text(size =8), 
    legend.position = "bottom", 
    legend.title = element_blank(),  
    legend.text = element_text(size = 8), 
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5)  
  ) 
print(PTERBB2_box)

mycomparisons <- list(c("Patient.primary", "Patient.brainmet"))


# Create the plot with Wilcoxon rank-sum test and multiple testing correction
PTERBB2_box_sig <- PTERBB2_box +
  stat_compare_means(
    comparisons = mycomparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    p.adjust.method = "BH",  # Benjamini-Hochberg adjustment
    label = "p.signif",      # Annotate plot with significance labels
    #label = "p.value",
    tip.length = 0.01,
    label.y = c(5.5),          #
    vjust = 0.7
  )

# Display the plot
print(PTERBB2_box_sig)

library(dplyr)

# Calculate the median for each group
median_values <- combined_data %>%
  group_by(Sample) %>%
  summarise(median_ERBB2 = median(ERBB2_counts, na.rm = TRUE))

# Print the median values
print(median_values)


egfr_box <- ggplot(combined_data, aes(x = Sample, y = EGFR_log2, fill = Sample)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6,outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), aes(color = Sample), alpha = 0.05, size = .5) +
  labs(x = "Sample", y = "EGFR log2+1 transcript counts", title = "EGFR expression Patient") +
  scale_fill_manual(values = c("Patient.primary" = "#C86050", "Patient.brainmet" = "#6685c2")) +
  scale_color_manual(values = c("Patient.primary" = "#C86050", "Patient.brainmet" = "#6685c2")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "transparent"), 
    panel.grid.minor = element_blank(), 
    plot.title = element_text(size = 14),  
    axis.text = element_text(size =8), 
    legend.position = "bottom",  
    legend.title = element_blank(),  
    legend.text = element_text(size = 8), 
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5)  
  )  +
  scale_y_continuous(limits = c(0, 4)) 
print(egfr_box)

mycomparisons <- list(c("Patient.primary", "Patient.brainmet"))

egfr_box_sig <- egfr_box +   stat_compare_means(
  comparisons = mycomparisons,
  method = "wilcox.test",
  method.args = list(exact = FALSE),
  p.adjust.method = "BH",  # Benjamini-Hochberg adjustment
  label = "p.signif",      # Annotate plot with significance labels
 # label = "p.value",
  tip.length = 0.01,
  label.y = c(3.5),          
  vjust = 0.7

)

print(egfr_box_sig)

library(dplyr)

# Calculate the median for each group
median_values <- combined_data %>%
  group_by(Sample) %>%
  summarise(median_EGFR = median(EGFR_counts, na.rm = TRUE))

# Print the median values
print(median_values)

