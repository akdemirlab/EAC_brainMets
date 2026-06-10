
library(ggplot2)
library(swimplot)
library(ggtext)
library(dplyr, warn.conflicts=FALSE)
library(reshape2)
library(grid)
library(plotly)
library(knitr)

# Basic plot
df <- Extended_brainmets_survival_datatypes_swimmerplot
df <- as.data.frame(df)

df.basic <- df %>%
  dplyr::select(`ID-dec`, End) %>%
  distinct(`ID-dec`, .keep_all = TRUE)

df.basic <- df.basic %>%
  mutate(id = as.character(`ID-dec`))
df.basic <- as.data.frame(df.basic)

# -----------------------------
# Clean + prepare data
# -----------------------------
# Standardize ID formatting
df$ID_dec <- gsub("\\s*\\+\\s*", " + ", df$`ID-dec`)
df$ID_dec <- as.character(df$ID_dec)

# Define desired order
id_order <- sprintf("P-%02d", 65:1)

# Apply factor levels
df$ID_dec <- factor(df$ID_dec, levels = id_order)

# -----------------------------
# Build swimmer plot
# -----------------------------
arm_plot <- swimmer_plot(
  df = df,
  id = "ID_dec",
  end = "End",
  id_order = id_order,
  name_fill = "LMD_status",
  name_col = "LMD_status",  # Changed: use name_col parameter
  alpha = 0.75,
  width = 0.8
)

# -----------------------------
# Apply fill AND outline colors
# -----------------------------
arm_plot <- arm_plot +
  scale_fill_manual(values = c("No" = "gray", "Yes" = "darkred")) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "darkred"))

# -----------------------------
# Create colored axis labels
# -----------------------------
id_status <- df %>%
  distinct(ID_dec, LMD_status) %>%
  mutate(
    label = ifelse(
      LMD_status == "Yes",
      paste0("<span style='color:darkred'>", ID_dec, "</span>"),
      paste0("<span style='color:gray'>", ID_dec, "</span>")
    )
  )

label_vec <- setNames(id_status$label, id_status$ID_dec)


AE_plot_dash <- arm_plot +
  geom_hline(yintercept = 25, linetype = "dashed", color = "#A883FD", size = 0.5)
AE_plot_dash

AE_plot_dash1 <- AE_plot_dash +
  geom_hline(yintercept = 6, linetype = "dashed", color = "#DC7454", size = 0.5)
AE_plot_dash1


AE_plot_1 <- AE_plot_dash1 +  # Or geom_segment(), etc. for swimmer plot
  theme_classic() +  # Use a clean theme first
  theme(
    panel.border = element_blank(),  # Removes box around plot
    axis.line = element_blank()      # Removes x and y axis lines
  )

AE_plot_1

AE_plot_1_yaxis <- AE_plot_1 +  scale_y_continuous(name = "Time since BM Dx (months)", breaks = seq(0,210,by=20))

library(ggplot2)

AE_plot_1_yaxis <- AE_plot_1 +
  scale_y_continuous(
    name = "Time since BM Dx (months)",
    breaks = seq(0, 210, by = 20)
  ) +
  theme(
    axis.title.x = element_text(size = 14),    # X-axis title
    axis.title.y = element_text(size = 14),    # Y-axis title
    axis.text.x  = element_text(size = 12),    # X-axis labels
    axis.text.y  = element_text(size = 12),    # Y-axis labels
    legend.title = element_text(size = 12),    # Legend title
    legend.text  = element_text(size = 12),    # Legend labels
    strip.text  = element_text(size = 13)      # Facet labels (if any)
  )

# Display the plot
AE_plot_1_yaxis
ggsave("swimmer_plot-LMD.png", plot = AE_plot_1_yaxis, width = 8, height = 10, dpi = 300)

# -----------------------------
# Apply colored labels to plot
# -----------------------------
final_plot <- AE_plot_1_yaxis +
  scale_y_discrete(labels = label_vec) +
  theme(
    axis.text.y = element_markdown()
  )

# -----------------------------
# Render plot
# -----------------------------
final_plot


