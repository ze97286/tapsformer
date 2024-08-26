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
bioc_packages <- c("ComplexHeatmap", "bsseq")
cran_packages <- c("cluster", "Rtsne", "ggplot2")
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)
sessionInfo()

load_and_create_bsseq <- function(base_dir, prefix) {
    sample_files <- list.files(path = base_dir, pattern = paste0("^", prefix, "_.*\\.rds$"), full.names = TRUE)
    if (length(sample_files) == 0) {
        stop("Error: No RDS files found with the given prefix.")
    }
    bsseq_list <- lapply(sample_files, function(file_path) {
        sample_data <- readRDS(file_path)
        sample_name <- gsub("\\.rds$", "", basename(file_path))

        # Removing or collapsing duplicate loci
        unique_loci <- !duplicated(paste(sample_data$chr, sample_data$pos, sep = ":"))
        sample_data <- sample_data[unique_loci, ]

        BSseq(
            chr = sample_data$chr,
            pos = sample_data$pos,
            M = as.matrix(sample_data$X),
            Cov = as.matrix(sample_data$N),
            sampleNames = sample_name
        )
    })

    combined_bsseq <- do.call(combineList, bsseq_list)
    return(combined_bsseq)
}

base_dir <- "/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg/oac"
tumour_bsseq <- load_and_create_bsseq(base_dir, "tumour")
print("tumour data loaded")

# Step 2: Subsample or filter data to reduce size
methylation_levels <- getMeth(tumour_bsseq, type = "raw")

# Remove rows with NA or NaN values
methylation_levels <- na.omit(methylation_levels)

# Optionally, filter out rows (CpG sites) with low variance
var_filter <- apply(methylation_levels, 1, var, na.rm = TRUE)

# Handle any potential NaN values in var_filter
var_filter <- var_filter[!is.na(var_filter)]

high_var_sites <- var_filter > quantile(var_filter, 0.9, na.rm = TRUE)  # Keep top 10% of most variable sites
methylation_levels_filtered <- methylation_levels[high_var_sites, ]

# Step 3: Proceed with clustering and visualization
dist_matrix <- dist(t(methylation_levels_filtered))  # Distance matrix for samples
hclust_res <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering

# Step 5: Save the heatmap to an SVG file
heatmap_file <- file.path(base_dir, "tumour_samples_heatmap.svg")
svg(heatmap_file, width = 8, height = 10)

heatmap <- Heatmap(methylation_levels_filtered,
                   cluster_columns = hclust_res, # Clustering applied to columns (samples)
                   show_row_dend = FALSE, 
                   show_column_names = TRUE,
                   column_title = "Hierarchical Clustering of Tumour Samples")

draw(heatmap)
dev.off()

print("heatmap saved")


# # Optional: Subset data for better visualization if too many CpGs
# # methylation_levels_subset <- methylation_levels[sample(1:nrow(methylation_levels), 1000),]

# heatmap <- Heatmap(methylation_levels,
#                    cluster_rows = hclust_res,
#                    show_row_dend = FALSE,
#                    show_column_names = FALSE,
#                    column_title = "Hierarchical Clustering of Tumour Samples")

# draw(heatmap)
# dev.off()

# print("Hierarchical Clustering of Tumour Samples saved")

# # Step 4: Optional - Dimensionality Reduction and Visualization
# # Save PCA plot to SVG
# pca_res <- prcomp(t(methylation_levels), scale. = TRUE)
# pca_df <- data.frame(pca_res$x, Sample = colnames(methylation_levels))
# pca_file <- file.path(base_dir,"tumour_samples_pca.svg")
# svg(pca_file, width = 8, height = 8)
# ggplot(pca_df, aes(PC1, PC2, label = Sample)) +
#     geom_point(size = 3) +
#     geom_text(hjust = 1.5, vjust = 1.5) +
#     theme_minimal() +
#     labs(title = "PCA of Tumour Samples")
# dev.off()

# print("pca saved")

# # Save t-SNE plot to SVG
# tsne_res <- Rtsne(t(methylation_levels), dims = 2, perplexity = 30)
# tsne_df <- data.frame(tsne_res$Y, Sample = colnames(methylation_levels))
# tsne_file <- file.path(base_dir,"tumour_samples_tsne.svg")
# svg(tsne_file, width = 8, height = 8)
# ggplot(tsne_df, aes(X1, X2, label = Sample)) +
#     geom_point(size = 3) +
#     geom_text(hjust = 1.5, vjust = 1.5) +
#     theme_minimal() +
#     labs(title = "t-SNE of Tumour Samples")
# dev.off()

# print("t-SNE saved")
