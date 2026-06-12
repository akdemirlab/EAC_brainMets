# run the heatmap script first as this script requires ouputs from the scWGS_heatmap script

primary_cells1 <- dat_seg_filtered[Group == "Patient-Primary", ]
metastasis_cells1 <- dat_seg_filtered[Group == "Patient-BrainMet", ]

library(ggplot2)
library(dplyr)


primary_cells1_df <- as.data.frame(primary_cells1)

library(reshape2)
primary_plot_data1 <- melt(
  primary_cells1_df,
  variable.name = "Bin",       
  value.name = "Ratio"        
)

primary_plot_data1$Chromosome <- rep(v_adjusted, each = nrow(primary_cells1_df))

rownames(primary_plot_data1) <- paste0("V", seq_len(nrow(primary_plot_data1)))

metastasis_cells1_df <- as.data.frame(metastasis_cells1)

library(reshape2)
metastasis_plot_data1 <- melt(
  metastasis_cells1_df,
  variable.name = "Bin",        
  value.name = "Ratio"       
)

metastasis_plot_data1$Chromosome <- rep(v_adjusted, each = nrow(metastasis_cells1_df))

rownames(metastasis_plot_data1) <- paste0("V", seq_len(nrow(metastasis_plot_data1)))


df <- primary_plot_data1

filter_average_ratio_per_bin <- function(df, Ratio) {
  if (!is.data.frame(df)) {
    stop("Input df must be a data frame.")
  }
    if (!(Ratio %in% names(df))) {
    stop("Specified Ratio column does not exist in the dataframe.")
  }
  if (!("Bin" %in% names(df))) {
    stop("The column 'Bin' does not exist in the dataframe.")
  }
  
  df <- df %>%
    mutate(High_Ratio = ifelse(!!sym(Ratio) > 3, TRUE, FALSE))  
  
  #Calculate the average Ratio for each bin (all rows included)
  avg_result <- df %>%
    group_by(Bin) %>%
    summarize(
      Average_Ratio = mean(!!sym(Ratio), na.rm = TRUE),  
      .groups = 'drop'
    )
  
  #Merge the averages back into the original dataframe
  df_with_avg <- df %>%
    left_join(avg_result, by = "Bin") %>%
    mutate(
      Final_Ratio = Average_Ratio  # Use the bin average for all rows
    )
  
  return(df_with_avg)
}




df_with_avg <- filter_average_ratio_per_bin(df, "Ratio")

# Plotting code
primary_plot <- ggplot(df_with_avg, aes(x = Bin)) + 
  geom_point(aes(y = Final_Ratio, color = ifelse(Final_Ratio >= 1.25 & Final_Ratio <= 6, "red", "gray")), 
             size = 0.1, alpha = 0.4) +
  scale_color_identity() +  
  geom_point(data = df_with_avg %>% filter(High_Ratio), 
             aes(y = Ratio - 1.75), color = "red", size = 0.2, alpha = 0.8) +
  
  scale_y_continuous(limits = c(0.5, 6)) + 
  labs(
    title = "Coverage Ratio: Primary",
    x = "Genomic Bins",
    y = "Coverage Ratio"
  ) + 
  theme_minimal(base_size = 14) + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black")
  )

ggsave("coverage_Ratio_primary.png", plot = primary_plot, width = 12, height = 3, dpi = 300)



df <- metastasis_plot_data1

filter_average_ratio_per_bin <- function(df, Ratio) {
  if (!is.data.frame(df)) {
    stop("Input df must be a data frame.")
  }
  
  if (!(Ratio %in% names(df))) {
    stop("Specified Ratio column does not exist in the dataframe.")
  }
  if (!("Bin" %in% names(df))) {
    stop("The column 'Bin' does not exist in the dataframe.")
  }
  
  df <- df %>%
    mutate(High_Ratio = ifelse(!!sym(Ratio) > 3, TRUE, FALSE))
  
  # Calculate the average Ratio for each bin (all rows included)
  avg_result <- df %>%
    group_by(Bin) %>%
    summarize(
      Average_Ratio = mean(!!sym(Ratio), na.rm = TRUE), 
      .groups = 'drop'
    )
  
  # Merge the averages back into the original dataframe
  df_with_avg <- df %>%
    left_join(avg_result, by = "Bin") %>%
    mutate(
      Final_Ratio = Average_Ratio  
    )
  
  return(df_with_avg)
}

df_with_avg <- filter_average_ratio_per_bin(df, "Ratio")


metastasis_plot <- ggplot(df_with_avg, aes(x = Bin)) + 
  geom_point(aes(y = Final_Ratio, color = ifelse(Final_Ratio >= 1.25 & Final_Ratio <= 3, "red", "gray")), 
             size = 0.1, alpha = 0.4) +
  scale_color_identity() + 
    geom_point(data = df_with_avg %>% filter(High_Ratio), 
             aes(y = Ratio - 1.75), color = "red", size = 0.2, alpha = 0.8) +
  
  scale_y_continuous(limits = c(0.5, 6)) +  
  labs(
    title = "Coverage Ratio: Metastasis",
    x = "Genomic Bins",
    y = "Coverage Ratio"
  ) + 
  theme_minimal(base_size = 14) + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black")
  )

ggsave("coverage_Ratio_metastasis.png", plot = metastasis_plot, width = 12, height = 3, dpi = 300)
