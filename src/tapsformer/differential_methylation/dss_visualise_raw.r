start_time <- Sys.time()

install_and_load <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Function for Bioconductor packages
bioc_install_and_load <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, update = FALSE)
  sapply(pkg, require, character.only = TRUE)
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Define packages
bioc_packages <- c("DSS", "GenomicRanges", "bsseq", "org.Hs.eg.db",
                   "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub")
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "ggplot2", "svglite","pheatmap")

# Install and load packages
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)

# Print loaded package versions
sessionInfo()

base_dir <- "/users/zetzioni/sharedscratch/tapsformer/data/methylation"
plot_dir <- file.path(base_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
flog.appender(appender.file("/users/zetzioni/sharedscratch/logs/dss_visualise_raw.log"))

flog.info("Starting DSS analysis script")

# Load the data
tumour_rds_path <- file.path(base_dir, "tumour_data.rds")
control_rds_path <- file.path(base_dir, "control_data.rds")
tumour_data <- readRDS(tumour_rds_path)
control_data <- readRDS(control_rds_path)

# Function to safely create a plot
safe_plot <- function(filename, plot_func) {
  tryCatch({
    svglite(filename, width = 10, height = 8)
    plot_func()
    dev.off()
    flog.info(paste("Plot saved as", filename))
  }, error = function(e) {
    flog.error(paste("Error creating plot:", filename, "-", conditionMessage(e)))
  })
}

# Combine and clean data for histogram plot
histogram_data <- rbind(
  data.table(sample_type = "tumour", beta = tumour_data$beta),
  data.table(sample_type = "Control", beta = control_data$beta)
)
histogram_data <- histogram_data[!is.na(beta) & !is.nan(beta)]

# Log the number of removed values
flog.info(paste("Number of NA/NaN beta values removed:", nrow(histogram_data) - nrow(histogram_data[!is.na(beta) & !is.nan(beta)])))

flog.info("preparing beta histogram")

# Create the histogram plot
beta_histogram <- ggplot(histogram_data, aes(x = beta, color = sample_type, fill = sample_type)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = c("Control" = "blue", "tumour" = "red")) +
  scale_fill_manual(values = c("Control" = "blue", "tumour" = "red")) +
  labs(title = "Distribution of Unmethylation Levels",
       x = "Unmethylation Level (Beta Value)",
       y = "Density",
       color = "Sample Type",
       fill = "Sample Type") +
  theme_minimal()

safe_plot(
  file.path(plot_dir, "beta_value_distribution.svg"),
  function() {
    print(beta_histogram)
  }
)

summary_stats <- histogram_data[, .(
  min = min(beta, na.rm = TRUE),
  q1 = quantile(beta, 0.25, na.rm = TRUE),
  median = median(beta, na.rm = TRUE),
  mean = mean(beta, na.rm = TRUE),
  q3 = quantile(beta, 0.75, na.rm = TRUE),
  max = max(beta, na.rm = TRUE),
  na_count = sum(is.na(beta)),
  total_count = .N
), by = sample_type]

print(summary_stats)
flog.info("Summary statistics of beta values:")
flog.info(capture.output(print(summary_stats)))
flog.info("saved beta histogram")

# Function to prepare data
prepare_data <- function(tumour_data, control_data, window_size = 1e4) {
  # Ensure data is sorted
  setkey(tumour_data, chr, pos)
  setkey(control_data, chr, pos)
  
  # Create windows
  tumour_data[, window := floor(pos / window_size) * window_size, by = chr]
  control_data[, window := floor(pos / window_size) * window_size, by = chr]
  
  # Calculate average methylation for each window
  tumour_avg <- tumour_data[, .(avg_beta = mean(beta, na.rm = TRUE)), by = .(chr, window)]
  control_avg <- control_data[, .(avg_beta = mean(beta, na.rm = TRUE)), by = .(chr, window)]
  
  # Merge and calculate difference (tumour - control)
  merged <- merge(tumour_avg, control_avg, by = c("chr", "window"), suffixes = c("_tumour", "_control"))
  merged[, diff := avg_beta_tumour - avg_beta_control]
  
  # Add cumulative position for plotting
  merged[, chr_num := as.numeric(sub("chr", "", chr))]
  merged <- merged[order(chr_num, window)]
  merged[, cum_pos := cumsum(c(0, diff(window))) + 1e8 * (chr_num - 1)]
  
  return(merged)
}

# Prepare data
window_size <- 1e4
methylation_data <- prepare_data(tumour_data, control_data, window_size)

# Remove NA values
methylation_data <- methylation_data[!is.na(diff)]

# Calculate chromosome midpoints for x-axis labels
chr_midpoints <- methylation_data[, .(mid_pos = median(cum_pos)), by = chr]

# Calculate color scale limits based on data
color_limits <- c(min(methylation_data$diff), max(methylation_data$diff))

# Create the heatmap plot
flog.info("preparing heatmap plot")

p <- ggplot(methylation_data, aes(x = cum_pos, y = 1, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                     midpoint = 0, 
                     limits = c(-max(abs(methylation_data$diff)), max(abs(methylation_data$diff))),
                     name = "Unmethylation\nDifference") + 
  scale_x_continuous(
    breaks = chr_midpoints$mid_pos,
    labels = chr_midpoints$chr,
    expand = c(0, 0)
  ) +
  labs(title = "Unmethylation Difference Across the Genome (Tumour - Control)",
       subtitle = paste("Range:", round(min(methylation_data$diff), 4), "to", round(max(methylation_data$diff), 4)),
       x = "Chromosome", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the heatmap plot
safe_plot(
  file.path(plot_dir, "unmethylation_difference_heatmap.svg"),
  function() {
    print(p)
  }
)
flog.info("saved heatmap plot")
print(summary(methylation_data$diff))

# Aggregated Heatmap
flog.info("preparing aggregated heatmap plot")
aggregate_data <- function(data, window_size = 1e4) {
  data[, window := floor(pos / window_size) * window_size, by = chr]
  data_avg <- data[, .(avg_beta = mean(beta, na.rm = TRUE)), by = .(chr, window)]
  return(data_avg)
}

tumour_avg <- aggregate_data(tumour_data)
control_avg <- aggregate_data(control_data)

# Merge the aggregated data
heatmap_data <- merge(tumour_avg, control_avg, by = c("chr", "window"), suffixes = c("_tumour", "_control"))
heatmap_data[, diff := avg_beta_tumour - avg_beta_control]

# Check for NA values and remove them
heatmap_data <- heatmap_data[!is.na(diff)]

# Subsample if necessary to reduce the matrix size
max_size <- 10000  # Adjust based on your system's capability
if (nrow(heatmap_data) > max_size) {
  set.seed(123)
  heatmap_data <- heatmap_data[sample(.N, max_size)]
}

# Create the matrix for the heatmap
heatmap_matrix <- dcast(heatmap_data, chr + window ~ ., value.var = "diff")
heatmap_matrix <- as.matrix(heatmap_matrix[, -c("chr", "window"), with = FALSE])

# Ensure no NA values in the matrix
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Check if matrix has at least 2 rows/columns
if (nrow(heatmap_matrix) >= 2 && ncol(heatmap_matrix) >= 2) {
  # Plot the heatmap
  safe_plot(
    file.path(plot_dir, "aggregated_unmethylation_heatmap.svg"),
    function() {
      pheatmap(heatmap_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Aggregated Unmethylation Difference Heatmap (Tumour - Control)")
    }
  )
  flog.info("saved aggregated heatmap plot")
} else {
  flog.error("Aggregated heatmap data has less than 2 rows/columns; skipping plot.")
}

# Boxplot of Beta Values by Chromosome
flog.info("preparing boxplot of beta values by chromosome plot")

# Function to subsample data
subsample_data <- function(data, max_points_per_group) {
  sampled_data <- data[, .SD[sample(.N, min(.N, max_points_per_group))], by = chr]
  return(sampled_data)
}

# Subsample data to reduce file size
max_points_per_group <- 100000  # Adjust as needed
subsampled_tumour_data <- subsample_data(tumour_data, max_points_per_group)
subsampled_control_data <- subsample_data(control_data, max_points_per_group)

# Combine subsampled data
boxplot_data <- rbind(
  data.table(sample_type = "tumour", chr = subsampled_tumour_data$chr, beta = subsampled_tumour_data$beta),
  data.table(sample_type = "Control", chr = subsampled_control_data$chr, beta = subsampled_control_data$beta)
)
boxplot_data <- boxplot_data[!is.na(beta) & !is.nan(beta)]

# Create the boxplot
boxplot <- ggplot(boxplot_data, aes(x = chr, y = beta, fill = sample_type)) +
  geom_boxplot(outlier.size = 0.5, lwd = 0.25) +
  scale_fill_manual(values = c("Control" = "blue", "tumour" = "red")) +
  labs(title = "Boxplot of Unmethylation Levels by Chromosome",
       x = "Chromosome",
       y = "Unmethylation Level (Beta Value)",
       fill = "Sample Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the boxplot
safe_plot(
  file.path(plot_dir, "beta_value_boxplot_by_chr.svg"),
  function() {
    print(boxplot)
  }
)

flog.info("saved boxplot of beta values by chromosome plot")

# PCA Plot
set.seed(123)
sampled_indices <- sample(seq_len(nrow(tumour_data)), size = 10000)

flog.info("preparing pca plot")

# Combine and clean data for PCA
pca_data <- rbind(
  data.table(sample = "tumour", beta = tumour_data$beta[sampled_indices]),
  data.table(sample = "Control", beta = control_data$beta[sampled_indices])
)

# Ensure no NA or infinite values
pca_data <- pca_data[!is.na(beta) & !is.infinite(beta)]

# Check variance in the data
flog.info(paste("Variance in beta values:", var(pca_data$beta)))

# Perform PCA
pca <- prcomp(pca_data[, .(beta)], center = TRUE, scale. = TRUE)

# Inspect PCA object
flog.info("PCA summary:")
flog.info(capture.output(summary(pca)))

# Extract PC1 and PC2 for plotting
pca_results <- as.data.table(pca$x)
pca_results$sample <- pca_data$sample

# Check if PC2 exists in pca_results
if (!"PC2" %in% names(pca_results)) {
  flog.error("PC2 not found in PCA results; cannot create PCA plot.")
} else {
  # Create PCA plot
  pca_plot <- ggplot(pca_results, aes(x = PC1, y = PC2, color = sample)) +
    geom_point() +
    labs(title = "PCA of Unmethylation Data",
         x = "Principal Component 1",
         y = "Principal Component 2",
         color = "Sample Type") +
    theme_minimal()

  safe_plot(
    file.path(plot_dir, "pca_plot.svg"),
    function() {
      print(pca_plot)
    }
  )
  flog.info("saved pca plot")
}

flog.info("Visualisation script completed successfully")

end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))