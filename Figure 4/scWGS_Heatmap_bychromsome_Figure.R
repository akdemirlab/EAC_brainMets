

packages <- c("tidyverse", "fs", "shiny",
              "here", "flexdashboard", "cowplot",
              "janitor", "ape", "ggsci", "plotly",
              "amap", "paletteer", "scales", "umap","uwot",
              "flsa", "BiocManager", "DT", "devtools", "Metrics","scquantum")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
    suppressWarnings(library(x, character.only = TRUE))
  }
})

bioc_packages <- c("ggtree", "ComplexHeatmap")

lapply(bioc_packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, update = FALSE)
    suppressWarnings(library(x, character.only = TRUE))
  }
})

git_packages<- c("Rphenograph")
lapply(git_packages, FUN = function(x) {
  if(!require(x, character.only = TRUE)){
    devtools::install_github("JinmiaoChenLab/Rphenograph")
    library(x,character.only = TRUE)
  }
  
})

options(bitmapType='cairo')
theme_set(theme_cowplot())

#core script for generating heatmap
# reading data

dat_seg <- uber_cnv_1_seg
# Setting up dat_seg_cp by removing specific columns
dat_seg_cp <- dat_seg %>% dplyr::select(-chrom, -chrompos, -abspos)

# Setup dat_seg_s and dat_seg_t for UMAP calculation
dat_seg_s <- dat_seg %>% dplyr::select(-chrom, -chrompos, -abspos)
dat_seg_t <- as.data.frame(t(dat_seg_s))

# Number of cells in the original dataset
n_cells <- dat_seg %>% dplyr::select(-abspos, -chrompos, -chrom) %>% ncol()

# UMAP calculation with dynamic n_neighbors
umap_n_neighbors <- if (n_cells < 30) { n_cells - 1 } else { 30 }
set.seed(42)
dat_umap <- uwot::umap(dat_seg_t, metric = "manhattan", n_threads = 20, min_dist = 0, n_neighbors = umap_n_neighbors)

# Prepare UMAP dataframe
umap_df <- as.data.frame(dat_umap) %>% dplyr::rename("umap1" = "V1", "umap2" = "V2")
rownames(umap_df) <- rownames(dat_seg_t)

# Set k for Rphenograph clustering
k_param <- if (n_cells < 30) { n_cells - 1 } else { 30 }
pheno <- Rphenograph::Rphenograph(umap_df, k = k_param)
cl <- tibble(cell = rownames(umap_df), cluster = igraph::membership(pheno[[2]])) %>% arrange(cluster)

# Order cells by clustering result and filter out 500 cells
dat_seg_order <- dat_seg_t[cl$cell,]
dat_seg_filtered <- dat_seg_order[1:(nrow(dat_seg_order) - 500), ]  

# Number of cells after filtering
n_cells_filtered <- ncol(dat_seg_filtered)


# List of chromosomes to plot
chromosomes_to_plot <- c(7, 8, 17, 1, 5, 11, 20)

# Loop through each chromosome and generate heatmaps
for (chr in chromosomes_to_plot) {
  
  # Subset data for the current chromosome
  chr_indices <- which(dat_seg$chrom == chr)  
  
  # Filter the data to include only the selected chromosome
  dat_chr_filtered <- dat_seg_filtered[, chr_indices]  
  
  # Calculate chromosome lengths 
  chr_lengths_chr <- ncol(dat_chr_filtered)  
  
  # Create positions for the chromosome labels 
  chr_l_means <- round(seq(1, chr_lengths_chr, length.out = chr_lengths_chr))
  
  # Create the vector `v` based on chromosome labels
  chrom_names <- as.character(chr)  
  v_chr <- rep("", chr_lengths_chr)
  v_chr[chr_l_means] <- chrom_names  
  
  v_chr[is.na(v_chr)] <- ""
  
  # Ensure the adjusted annotation vector matches the filtered data
  v_adjusted <- v_chr  
  
  # generate the binary chromosome annotation for the filtered subset
  chr_binary_chr <- rep(c(2, 1), length.out = chr_lengths_chr)  
  chr_df_chr <- data.frame(chr = chr_binary_chr) 
  
  # Create the HeatmapAnnotation for the specific chromosome
  chr_bar_chr <- HeatmapAnnotation(
    chr_text = anno_text(v_adjusted, gp = gpar(fontsize = 8)), 
    df = chr_df_chr, 
    show_legend = FALSE, 
    which = "column", 
    col = list(chr = c("1" = "grey88", "2" = "black"))
  )
  
  ####
  # Define groups based on cell name patterns
  Group <- rep("Patient-Primary", nrow(dat_seg_filtered))
  Group[grep("2", rownames(dat_seg_filtered), ignore.case = TRUE)] <- "Patient-BrainMet"
  
  # Define the colors for each group
  group_colors <- c("Patient-Primary" = "purple", "Patient-BrainMet" = "yellow")
  
  # Create the annotation using the colors
  heat_row_col <- rowAnnotation(
    df = data.frame(Group = Group),
    col = list(Group = group_colors)
  )
  ####
  
  # Generate the heatmap for the specific chromosome
  ht_chr <- Heatmap(as.matrix(log2(dat_chr_filtered + 1e-3)),
                    cluster_columns = FALSE,
                    cluster_rows = TRUE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    use_raster = TRUE,
                    left_annotation = heat_row_col,
                    col = circlize::colorRamp2(breaks = c(-2, 0.1, 2), c("dodgerblue3", "white", "firebrick3")),
                    heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                                title_gp = gpar(fontsize = 15, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 12)),
                    column_title = paste("Chromosome", chr)  
  )
  
  pdf(paste0("Heatmap_chr", chr, ".pdf"))
  print(ht_chr)
  dev.off()
}