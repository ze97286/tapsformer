install_and_load <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) {
        install.packages(new.pkg, dependencies = TRUE)
    }
    sapply(pkg, require, character.only = TRUE)
}

bioc_install_and_load <- function(pkg) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) {
        BiocManager::install(new.pkg, update = FALSE)
    }
    sapply(pkg, require, character.only = TRUE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))
bioc_packages <- c("ComplexHeatmap","bsseq")
cran_packages <- c("cluster", "Rtsne","ggplot2")
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)
sessionInfo()

base_dir <- "/Users/zoharetzioni/Downloads/methylation_clustering"
tumour_bsseq <- load_and_create_bsseq(base_dir, "tumour")
methylation_levels <- bsseq::getMeth(tumour_bsseq, type = "raw")

print("tumour data loaded")

# Step 2: Perform hierarchical clustering
# Compute distances and clustering
dist_matrix <- dist(t(methylation_levels))  # Distance matrix
hclust_res <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering

# Step 3: Save the heatmap to an SVG file
heatmap_file <- file.path(base_dir,"tumour_samples_heatmap.svg")
svg(heatmap_file, width = 8, height = 10)
print("heatmap saved")


# Optional: Subset data for better visualization if too many CpGs
# methylation_levels_subset <- methylation_levels[sample(1:nrow(methylation_levels), 1000),]

heatmap <- Heatmap(methylation_levels,
                   cluster_rows = hclust_res, 
                   show_row_dend = FALSE, 
                   show_column_names = FALSE,
                   column_title = "Hierarchical Clustering of Tumour Samples")

draw(heatmap)
dev.off()

print("Hierarchical Clustering of Tumour Samples saved")

# Step 4: Optional - Dimensionality Reduction and Visualization
# Save PCA plot to SVG
pca_res <- prcomp(t(methylation_levels), scale. = TRUE)
pca_df <- data.frame(pca_res$x, Sample = colnames(methylation_levels))
pca_file <- file.path(base_dir,"tumour_samples_pca.svg")
svg(pca_file, width = 8, height = 8)
ggplot(pca_df, aes(PC1, PC2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(hjust = 1.5, vjust = 1.5) +
    theme_minimal() +
    labs(title = "PCA of Tumour Samples")
dev.off()

print("pca saved")

# Save t-SNE plot to SVG
tsne_res <- Rtsne(t(methylation_levels), dims = 2, perplexity = 30)
tsne_df <- data.frame(tsne_res$Y, Sample = colnames(methylation_levels))
tsne_file <- file.path(base_dir,"tumour_samples_tsne.svg")
svg(tsne_file, width = 8, height = 8)
ggplot(tsne_df, aes(X1, X2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(hjust = 1.5, vjust = 1.5) +
    theme_minimal() +
    labs(title = "t-SNE of Tumour Samples")
dev.off()

print("t-SNE saved")