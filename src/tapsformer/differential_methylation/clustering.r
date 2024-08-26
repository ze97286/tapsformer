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
cran_packages <- c("cluster", "Rtsne", "ggplot2", "ggrepel")
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)
sessionInfo()

args <- commandArgs(trailingOnly = TRUE)
suffix <- args[1]
type_prefix <- args[2]

base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg",suffix)
subsampled_file <- file.path(base_dir,paste(type_prefix,"_subsampled_methylation_levels.rds",sep=""))
print(subsampled_file)

if (file.exists(subsampled_file)) {
    # Load the subsampled data
    print("loading from subsampled file")
    methylation_levels_subset <- readRDS(subsampled_file)
    print("Loaded subsampled methylation data from file.")
} else {
    # Step 1: Modify the function to detect and handle duplicate loci before creating BSseq object
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

    print("creating data set")
    type_bsseq <- load_and_create_bsseq(base_dir, type_prefix)
    print(type_prefix,"data loaded")

    # Step 2: Subsample CpG sites to reduce the size of the dataset
    set.seed(123)  # For reproducibility
    num_sites_to_keep <- 100000  

    # Randomly select a subset of CpG sites
    total_sites <- nrow(type_bsseq)
    subset_indices <- sample(1:total_sites, min(num_sites_to_keep, total_sites))

    methylation_levels <- getMeth(type_bsseq, type = "raw")
    methylation_levels_subset <- methylation_levels[subset_indices, ]

    # Step 3: Filter out rows with NA or NaN values
    methylation_levels_subset <- na.omit(methylation_levels_subset)

    # Save the subsampled data to disk
    saveRDS(methylation_levels_subset, subsampled_file)
    print("Subsampled methylation data saved to file.")
}

# Step 4: Proceed with clustering and visualization on the subset
dist_matrix <- dist(t(methylation_levels_subset))  # Distance matrix for samples
hclust_res <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering

# Step 5: Save the heatmap to an SVG file
heatmap_file <- file.path(base_dir,paste(type_prefix,"_samples_heatmap.svg",sep=""))
svg(heatmap_file, width = 8, height = 10)

heatmap <- Heatmap(methylation_levels_subset,
                   cluster_columns = hclust_res, # Clustering applied to columns (samples)
                   show_row_dend = FALSE, 
                   show_column_names = TRUE,
                   column_title = sprintf("Hierarchical Clustering of %s Samples",type_prefix))

draw(heatmap)
dev.off()
print("heatmap saved")

# Step 6: PCA plot with ggrepel to avoid label overlap
pca_res <- prcomp(t(methylation_levels_subset), scale. = FALSE)
pca_df <- data.frame(pca_res$x, Sample = colnames(methylation_levels_subset))
pca_file <- file.path(base_dir,paste(type_prefix,"_samples_pca.svg",sep=""))
svg(pca_file, width = 8, height = 8)
ggplot(pca_df, aes(PC1, PC2, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel() +  # Use geom_text_repel to avoid overlap
    theme_minimal() +
    labs(title = sprintf("PCA of %s Samples",type_prefix))
dev.off()
print("pca saved")

# Step 7: t-SNE with adjusted perplexity and ggrepel to avoid label overlap
# Get the number of samples
num_samples <- ncol(methylation_levels_subset)

# Set a very conservative perplexity value based on the number of samples
perplexity_value <- min(10, floor((num_samples - 1) / 3))  # A safer approach, ensuring perplexity is much smaller than the number of samples

# Run t-SNE with the adjusted perplexity
tsne_res <- Rtsne(t(methylation_levels_subset), dims = 2, perplexity = perplexity_value)
tsne_df <- data.frame(tsne_res$Y, Sample = colnames(methylation_levels_subset))
tsne_file <- file.path(base_dir,paste(type_prefix,"_samples_tsne.svg", sep=""))
svg(tsne_file, width = 8, height = 8)
ggplot(tsne_df, aes(X1, X2, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel() + 
    theme_minimal() +
    labs(title = sprintf("t-SNE of %s Samples",type_prefix))
dev.off()

print("t-SNE saved")