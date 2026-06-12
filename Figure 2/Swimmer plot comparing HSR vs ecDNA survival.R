
library(swimplot)
library(ggplot2)
library(dplyr, warn.conflicts=FALSE)   
library(reshape2)
library(grid)
library(plotly) 
library(knitr)


df <- ERBB2_amptype_survivaldata
df <- as.data.frame(df)

df.basic <- df %>%
  # Get just the subject and response time columns
  dplyr::select(ID, End) %>%
  # Remove duplicate rows based on 'id'
  distinct(ID, .keep_all = TRUE)

df.basic <- df.basic %>%
  mutate(id = as.character(ID))
df.basic <- as.data.frame(df.basic)

df$ERBB2_Amplification <- gsub("\\s*\\+\\s*", " + ", df$`ERBB2 Amplification`)

df$ERBB2_Amplification_Type <- gsub("\\s*\\+\\s*", " + ", df$`ERBB2 Amplification Type`)


Swim <- swimmer_plot(df = df, id = "ID", end = "End", name_fill = "ERBB2_Amplification_Type", width = 0.85, col= "black")

Swim

#change ID order

id_order = c('P-63', 'P-62', 'P-57', 'P-52', 'P-51', 'P-48', 'P-43', 'P-41', 'P-40', 'P-39', 'P-36', 'P-24', 'P-21', 'P-10', 'P-06', 'P-01')

arm_plot <- swimmer_plot(df=df,id='ID',end='End',name_fill='ERBB2_Amplification_Type',
                       id_order= id_order,col="black",width=.8)
arm_plot


AE_plot1 <- arm_plot + swimmer_points(df_points=
   df,id='ID',time='End', name_fill = "Alive", size=4.5,fill='white')

AE_plot1

AE_plot <-  arm_plot +
  scale_fill_manual(name="Status",values=c("Alive" = "lightblue", "Death"="gray"))+
  scale_color_manual(name="Status",values=c("Alive" = "lightblue", "Death"="gray"))
AE_plot

AE_plot <-  arm_plot +
  scale_fill_manual(name="ERBB2_Amplification_Type",values=c("HSR" = "#BEA2FE", "ecDNA"="#b5585b"))+
  scale_color_manual(name="ERBB2_Amplification_Type",values=c("HSR" = "#BEA2FE", "ecDNA"="#b5585b"))
AE_plot

AE_plot_dash <- AE_plot +
  geom_hline(yintercept = 60, linetype = "dashed", color = "#BEA2FE", size = 0.5)
AE_plot_dash

AE_plot_dash1 <- AE_plot_dash +
  geom_hline(yintercept = 8, linetype = "dashed", color = "#b5585b", size = 0.5)
AE_plot_dash1



AE_plot_dash_clean <- AE_plot_dash1 +  t
  theme_classic() +  
  theme(
    panel.border = element_blank(),  
    axis.line = element_blank()      
  )

AE_plot_dash_clean

AE_plot_dash_clean_yaxis <- AE_plot_dash_clean +  scale_y_continuous(name = "Time since BM Dx (months)", breaks = seq(0,210,by=20))
                      
library(ggplot2)

AE_plot_dash_clean_yaxis <- AE_plot_dash_clean +
  scale_y_continuous(
    name = "Time since BM Dx (months)",
    breaks = seq(0, 210, by = 20)
  ) +
  theme(
    axis.title.x = element_text(size = 14),    # X-axis title
    axis.title.y = element_text(size = 14),    # Y-axis title
    axis.text.x  = element_text(size = 14),    # X-axis labels
    axis.text.y  = element_text(size = 16),    # Y-axis labels
    legend.title = element_text(size = 12),    # Legend title
    legend.text  = element_text(size = 12),    # Legend labels
    strip.text  = element_text(size = 13)      # Facet labels (if any)
  )

# Display the plot
AE_plot_dash_clean_yaxis


                                                                                                        

ggsave("ecDNAvHSR_survival_swimmerplot.png", plot = AE_plot_dash_clean_yaxis, width = 8, height = 10, dpi = 300)

