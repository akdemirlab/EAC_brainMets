

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

# Set up the chromosome annotation
chr_lengths <- dat_seg %>% dplyr::select(abspos, chrom, chrompos) %>% 
  group_by(chrom) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::pull(n)

chr_binary <- rep(c(2, 1), 12)
chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))

# Calculate the cumulative sum of chromosome lengths and rowMeans
chr_rl_c <- c(1, cumsum(chr_lengths))
chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], b = chr_rl_c[2:length(chr_rl_c)])
chr_l_means <- round(rowMeans(chr_df))

# Define the chromosome names
chrom.names <- c(1:22, "X", "Y")

# Create a vector for chromosome labels
v <- vector(mode = "character", length = sum(chr_lengths))
v[chr_l_means] <- chrom.names
v[is.na(v)] <- ""

# Adjust chromosome label positions for filtered data
# Calculate new positions for filtered data based on the cumulative lengths of chromosomes
position_indices <- round(cumsum(chr_lengths) * n_cells_filtered / sum(chr_lengths))

# Create the adjusted vector of chromosome labels for filtered data
v_adjusted <- rep("", n_cells_filtered)
for (i in seq_along(chrom.names)) {
  if (i <= length(position_indices)) {
    v_adjusted[position_indices[i]] <- chrom.names[i]
  }
}

# Ensure there are no NAs or empty values
v_adjusted[is.na(v_adjusted)] <- ""

# Create the HeatmapAnnotation for chromosome labels
chr_bar <- HeatmapAnnotation(chr_text = anno_text(v_adjusted, gp = gpar(fontsize = 8)), 
                             df = chr, 
                             show_legend = FALSE, 
                             which = "column", 
                             col = list(chr = c("1" = "grey88", "2" = "black")))

# Peak and Cells annotations (no changes)
Peak <- rep("Peak_D", nrow(dat_seg_filtered))
Peak[grep(rownames(dat_seg_filtered), pattern = "_A_", ignore.case = TRUE)] <- "Peak_A"
Peak[grep(rownames(dat_seg_filtered), pattern = "Ctrl", ignore.case = TRUE)] <- "Ctrl"

Cells <- rep("preTX", nrow(dat_seg_filtered))
Cells[grep(rownames(dat_seg_filtered), pattern = "preTX", ignore.case = TRUE)] <- "MIDTX"

col <- c("blue", "orange")
names(col) <- c("MIDTX", "preTX")

pcol <- c("red", "purple", "green")
names(pcol) <- c("Ctrl", "Peak_D", "Peak_A")


# Create heatmap annotation for Peak
heat_row_col <- rowAnnotation(df = as.data.frame(Peak), col = list(Peak = pcol))

# Generate the heatmap
ht <- Heatmap(as.matrix(log2(dat_seg_filtered + 1e-3)),
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              use_raster = TRUE,
              left_annotation = heat_row_col,
              top_annotation = chr_bar,
              col = circlize::colorRamp2(breaks = c(-2, 0.1, 2), c("dodgerblue3", "white", "firebrick3")),
              heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 12)))





png("heatmap_scWGS_wholegenome.png", width = 2200, height = 2000, res = 300)

draw(ht,heatmap_legend_side = "right")


# Close the device to save the file
dev.off()