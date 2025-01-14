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
clusters <- as.numeric(args[3])

base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg", suffix)
subsampled_file <- file.path(base_dir, paste("clustering_",type_prefix, "_subsampled_methylation_levels.rds", sep = ""))
print(subsampled_file)

load_and_create_bsseq <- function(base_dir, prefix) {
    sample_files <- list.files(path = base_dir, pattern = paste0("^", prefix, "_.*\\.rds$"), full.names = TRUE)
    if (length(sample_files) == 0) {
        stop("Error: No RDS files found with the given prefix.")
    }

    bsseq_list <- lapply(sample_files, function(file_path) {
        sample_data <- readRDS(file_path)
        if (is.list(sample_data)) {
            if (all(c("chr", "pos", "X", "N") %in% names(sample_data))) {
                sample_name <- gsub("\\.rds$", "", basename(file_path))
                print(sample_name)
                unique_loci <- !duplicated(paste(sample_data$chr, sample_data$pos, sep = ":"))
                sample_data <- sample_data[unique_loci, ]

                return(BSseq(
                    chr = sample_data$chr,
                    pos = sample_data$pos,
                    M = as.matrix(sample_data$X),
                    Cov = as.matrix(sample_data$N),
                    sampleNames = sample_name
                ))
            } else {
                stop("Error: The loaded sample_data does not contain the required fields ('chr', 'pos', 'X', 'N').")
            }
        } else {
            stop("Error: sample_data is not a list. Check the contents of the RDS files.")
        }
    })

    combined_bsseq <- do.call(combineList, bsseq_list)
    return(combined_bsseq)
}

load_consistent_dmls <- function(dml_positions_file) {
    consistent_dmls <- readRDS(dml_positions_file)
    parsed_dmls <- data.frame(
        chr = sapply(strsplit(consistent_dmls, " "), function(x) x[1]),
        pos = as.numeric(sapply(strsplit(consistent_dmls, " "), function(x) x[2])),
        stringsAsFactors = FALSE
    )

    return(parsed_dmls)
}

if (file.exists(subsampled_file)) {
    print("loading from subsampled file")
    methylation_levels_subset <- readRDS(subsampled_file)
    print("Loaded subsampled methylation data from file.")
} else {
    dml_positions_file <- file.path(base_dir, "consistent_dmls.rds")
    parsed_dmls <- load_consistent_dmls(dml_positions_file)
    print("Loaded consistent DML positions")

    print("creating data set")
    type_bsseq <- load_and_create_bsseq(base_dir, type_prefix)
    print("data loaded")

    bsseq_chr <- as.character(seqnames(type_bsseq))
    bsseq_pos <- start(type_bsseq)
    matching_indices <- which(bsseq_chr %in% parsed_dmls$chr & bsseq_pos %in% parsed_dmls$pos)
    methylation_levels_subset <- getMeth(type_bsseq[matching_indices, ], type = "raw")
    methylation_levels_subset <- na.omit(methylation_levels_subset)
    saveRDS(methylation_levels_subset, subsampled_file)
    print("Subsampled methylation data saved to file.")
}

max_sites_to_plot <- 10000

if (nrow(methylation_levels_subset) > max_sites_to_plot) {
    subset_indices <- sample(1:nrow(methylation_levels_subset), max_sites_to_plot)
    methylation_levels_for_heatmap <- methylation_levels_subset[subset_indices, ]
} else {
    methylation_levels_for_heatmap <- methylation_levels_subset
}

dist_matrix <- dist(t(methylation_levels_for_heatmap))
hclust_res <- hclust(dist_matrix, method = "ward.D2")


heatmap_file <- file.path(base_dir, paste(type_prefix, "_samples_heatmap.svg", sep = ""))
svg(heatmap_file, width = 8, height = 10)
heatmap <- Heatmap(methylation_levels_for_heatmap,
    cluster_rows = FALSE,
    cluster_columns = hclust_res,
    show_row_dend = FALSE,
    show_column_names = TRUE,
    column_title = sprintf("Heatmap of %s Samples", type_prefix)
)

draw(heatmap)
dev.off()
print("Heatmap saved")

num_samples <- ncol(methylation_levels_subset)
perplexity_value <- min(10, floor((num_samples - 1) / 3)) 

tsne_res <- Rtsne(t(methylation_levels_subset), dims = 2, perplexity = perplexity_value)
tsne_df <- data.frame(tsne_res$Y, Sample = colnames(methylation_levels_subset))
tsne_file <- file.path(base_dir, paste(type_prefix, "_samples_tsne.svg", sep = ""))
svg(tsne_file, width = 8, height = 8)
ggplot(tsne_df, aes(X1, X2, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel() +
    theme_minimal() +
    labs(title = sprintf("t-SNE of %s Samples", type_prefix))
dev.off()

print("t-SNE saved")

num_clusters <- clusters
set.seed(123)
pca_res <- prcomp(t(methylation_levels_subset), scale. = FALSE)
pca_df <- data.frame(pca_res$x, Sample = colnames(methylation_levels_subset))
kmeans_res <- kmeans(pca_res$x[, 1:2], centers = num_clusters)
pca_df$Cluster <- as.factor(kmeans_res$cluster)

pca_clustered_file <- file.path(base_dir, paste(type_prefix, "_samples_pca_clusters.svg", sep = ""))
svg(pca_clustered_file, width = 8, height = 8)
ggplot(pca_df, aes(PC1, PC2, color = Cluster, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel() +
    theme_minimal() +
    labs(title = sprintf("PCA of %s Samples with Clusters", type_prefix))
dev.off()
print("PCA with clusters saved")

cluster_assignments_file <- file.path(base_dir, paste(type_prefix, "_cluster_assignments.csv", sep = ""))
write.csv(pca_df, cluster_assignments_file, row.names = FALSE)
print("Cluster assignments saved")
