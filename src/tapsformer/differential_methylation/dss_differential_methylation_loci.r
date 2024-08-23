# this script aims to identify differential methylation loci between tumour tissue and healthy control cfDNA.
# the main output is a bed file of positions of tumour hypomethylated loci.
# in addition it generated a number of visualisations for analysis.

start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript dml_analysis.R <delta> <p.threshold> <fdr.threshold> <smoothing> <suffix>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])

smoothing_arg <- tolower(args[4])
smoothing <- if (smoothing_arg == "true") TRUE else if (smoothing_arg == "false") FALSE else NA
if (is.na(smoothing)) {
  stop("Invalid smoothing argument. Please use 'TRUE' or 'FALSE'.")
}
suffix <- args[5]

# setup parallel processing
library(parallel)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(DSS)
  library(data.table)
  library(GenomicRanges)
  library(bsseq)
})

# initialise dirs
base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg", suffix)
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_%s_dml_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", suffix, delta, p.threshold, fdr.threshold))

# set up logging
flog.logger("dss_logger", appender.file(log_dir))
flog.threshold(INFO, name = "dss_logger")
flog.threshold(WARN, name = "dss_logger")
flog.threshold(ERROR, name = "dss_logger")
flog.info("Starting DSS DML analysis", name = "dss_logger")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

# generate a plot of top DMLs and their methylation levels in each sample
plot_top_DMLs <- function(top_hypo_dmls, combined_bsseq, output_dir) {
  # Get the raw methylation data
  methylation_data <- getMeth(combined_bsseq, type = "raw")
  sample_names <- colnames(methylation_data) # Get the actual sample names
  plot_data <- data.table(chr = top_hypo_dmls$chr, pos = top_hypo_dmls$pos)

  for (i in 1:nrow(top_hypo_dmls)) {
    dml <- top_hypo_dmls[i, ]

    # Find the corresponding methylation levels
    matching_indices <- which(seqnames(combined_bsseq) == dml$chr & start(combined_bsseq) == dml$pos)

    # Extract the methylation levels for this DML
    meth_levels <- methylation_data[matching_indices, , drop = FALSE]

    # Check if meth_levels has valid data
    if (is.null(meth_levels) || nrow(meth_levels) == 0 || ncol(meth_levels) == 0) {
      flog.warn(sprintf("No valid methylation data found for DML at %s:%d", dml$chr, dml$pos), name = "dss_logger")
      next # Skip to the next DML if no valid data is found
    }

    # Add methylation levels to plot_data with sample names as columns
    plot_data[i, (sample_names) := as.list(meth_levels)]
  }

  # Melt the data for plotting
  plot_data <- melt(plot_data,
    id.vars = c("chr", "pos"),
    variable.name = "Sample", value.name = "MethylationLevel"
  )

  # Add a new column to identify tumour and control samples
  plot_data[, SampleType := ifelse(grepl("^tumour", Sample), "Tumour", "Control")]

  # Check if plot_data is empty
  if (nrow(plot_data) == 0) {
    flog.error("Plot data is empty after processing. No plot will be generated.", name = "dss_logger")
    return(NULL)
  }

  # Generate and save the plot directly
  output_filename <- file.path(output_dir, "dml_methylation_plot.svg")

  tryCatch(
    {
      # Open a new svg device with larger size
      svglite::svglite(output_filename, width = 15, height = 10) # Adjust width and height for larger rectangles

      # Generate the plot
      plot <- ggplot(plot_data, aes(x = Sample, y = MethylationLevel, shape = SampleType, color = SampleType)) +
        geom_point(size = 2) +
        facet_wrap(~ chr + pos, scales = "free_y", ncol = 10) + # Set grid layout to 10 columns
        theme_bw() +
        labs(
          title = "Methylation Levels Across Samples for Each DML",
          x = "Sample", y = "Methylation Level"
        ) +
        scale_shape_manual(values = c(Tumour = 16, Control = 1)) + # Full circle for Tumour, open circle for Control
        scale_color_manual(values = c(Tumour = "blue", Control = "red")) + # Blue for Tumour, Red for Control
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 12), # Increase strip text size for readability
          plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust margins if needed
        )

      print(plot) # Ensure the plot is rendered

      # Close the svg device
      dev.off()
    },
    error = function(e) {
      # Handle any errors that occur during plotting
      flog.error(sprintf("Error generating plot: %s", conditionMessage(e)), name = "dss_logger")
      print(sprintf("Error generating plot: %s", conditionMessage(e)))
      dev.off() # Ensure device is closed even if there's an error
    }
  )

  return(output_filename)
}

sliding_window_filter <- function(dmls, window_size, min_cpgs = 3, consistency_threshold = 0.8) {
  dmls[, window := cut(pos, breaks = seq(min(pos), max(pos) + window_size, by = window_size))]
  windowed_dmls <- dmls[, .(
    mean_diff = mean(diff),
    mean_pval = mean(pval),
    n_dmls = .N,
    consistency = mean(sign(diff) == sign(mean(diff)))
  ), by = .(chr, window)]

  significant_windows <- windowed_dmls[
    abs(mean_diff) >= delta &
      mean_pval < p.threshold &
      n_dmls >= min_cpgs &
      consistency >= consistency_threshold
  ]

  dmls[chr %in% significant_windows$chr & window %in% significant_windows$window]
}

cross_validate_dmls <- function(bsseq_data, group1, group2, n_iterations = 10, subsample_fraction = 0.8,
                                delta, p.threshold, fdr.threshold, smoothing) {
  all_dmls <- list()

  for (i in 1:n_iterations) {
    # Subsample from each group separately
    subsample1 <- sample(group1, length(group1) * subsample_fraction)
    subsample2 <- sample(group2, length(group2) * subsample_fraction)
    subsample <- c(subsample1, subsample2)

    sub_bsseq <- bsseq_data[, subsample]

    dml_test <- DMLtest(sub_bsseq, group1 = subsample1, group2 = subsample2, smoothing = smoothing)
    dmls <- callDML(dml_test, delta = delta, p.threshold = p.threshold)

    dml_dt <- as.data.table(dmls)
    dml_dt[, significant_after_fdr := fdr < fdr.threshold]

    all_dmls[[i]] <- dml_dt[significant_after_fdr == TRUE, .(chr, pos)]
  }

  # Keep DMLs that appear in at least half of the iterations
  dml_counts <- table(unlist(lapply(all_dmls, function(x) paste(x$chr, x$pos))))
  consistent_dmls <- names(dml_counts[dml_counts >= n_iterations / 2])

  return(consistent_dmls)
}

adaptive_smooth <- function(bsseq_data, min_span = 200, max_span = 1000, min_cpgs = 20) {
  gr <- granges(bsseq_data)
  cpg_density <- width(disjoin(resize(gr, 1000, fix = "center")))
  span <- pmax(min_span, pmin(max_span, 1000 / cpg_density * min_cpgs))

  BSmooth(bsseq_data, h = span, ns = min_cpgs)
}

# this is the core function here, doing the DML analysis choosing hypomethylated DMLs and
# saving the output and visualisations.
perform_dml_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold,
                                 min_coverage, window_size, smoothing, cl, n_iterations = 10,
                                 subsample_fraction = 0.8, min_span = 200, max_span = 1000, min_cpgs = 20) {
  smoothing_string <- ifelse(smoothing, "smooth", "unsmooth")
  output_dir <- file.path(base_dir, sprintf(
    "dml_delta_%.2f_p_%.4f_fdr_%.2f_%s",
    delta, p.threshold, fdr.threshold, smoothing_string
  ))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (smoothing) {
    flog.info("Performing adaptive smoothing", name = "dss_logger")
    combined_bsseq <- tryCatch(
      {
        adaptive_smooth(combined_bsseq, min_span = min_span, max_span = max_span, min_cpgs = min_cpgs)
      },
      error = function(e) {
        flog.error(paste("Error in adaptive smoothing:", e$message), name = "dss_logger")
        stop("Adaptive smoothing failed")
      }
    )
    smoothing_for_dmltest <- FALSE
  } else {
    smoothing_for_dmltest <- FALSE
  }

  flog.info("Performing DML test", name = "dss_logger")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)

  # Cross-validation step
  flog.info("Performing cross-validation", name = "dss_logger")
  consistent_dmls <- cross_validate_dmls(combined_bsseq, group1, group2,
    n_iterations = n_iterations,
    subsample_fraction = subsample_fraction,
    delta = delta, p.threshold = p.threshold,
    fdr.threshold = fdr.threshold, smoothing = smoothing_for_dmltest
  )
  saveRDS(consistent_dmls, file.path(output_dir, "consistent_dmls.rds"))

  # Run the dml test with smoothing (TRUE/FALSE).
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = smoothing_for_dmltest)

  flog.info("Calling DMLs", name = "dss_logger")

  # Identify DML given the delta and pvalue threshold.
  dmls <- callDML(dml_test, delta = delta, p.threshold = p.threshold)

  dml_dt <- as.data.table(dmls)

  z_score <- qnorm(0.975) # Two-tailed 95% CI

  # Identify hypomethylated loci in the tumour and check significance
  dml_dt[, `:=`(
    hypo_in_tumour = diff < 0,
    significant_after_fdr = fdr < fdr.threshold,
    mean_methylation_diff = abs(diff),
    lower_ci = diff - (z_score * diff.se),
    upper_ci = diff + (z_score * diff.se),
    ci_excludes_zero = sign(diff - (z_score * diff.se)) == sign(diff + (z_score * diff.se)),
    consistent = paste(chr, pos) %in% consistent_dmls # Add this line to mark consistent DMLs
  )]

  # Select top hypomethylated DMLs
  top_hypo_dmls <- dml_dt[
    hypo_in_tumour == TRUE &
      significant_after_fdr == TRUE &
      mean_methylation_diff >= delta &
      ci_excludes_zero == TRUE &
      consistent == TRUE # Use the consistent flag here
  ]

  top_hypo_dmls <- sliding_window_filter(top_hypo_dmls, window_size)
  top_hypo_dmls[, composite_score := (abs(stat) * abs(diff)) / (diff.se * sqrt(fdr))]

  setorder(top_hypo_dmls, -composite_score)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmls, "stat", output_dir)

  # Tag differential methylation strength by quantile of areaStat
  top_hypo_dmls[, hypomethylation_strength := case_when(
    stat < thresholds$very_strong ~ "Very Strong",
    stat < thresholds$strong ~ "Strong",
    stat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # Save output
  fwrite(top_hypo_dmls[, .(chr, pos, pos + 2, stat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmls.bed"),
    sep = "\t"
  )

  flog.info("Creating visualizations", name = "dss_logger")
  strongest_dmls <- tail(top_hypo_dmls[order(top_hypo_dmls$stat), ], 50)
  str(strongest_dmls)
  plot_top_DMLs(strongest_dmls, combined_bsseq, output_dir)
  create_volcano_plot(dml_dt, diff_col = "diff", pval_col = "pval", output_dir)
  create_methylation_diff_plot(dml_dt, diff_col = "diff", output_dir)
  create_chromosome_coverage_plot(dml_dt, diff_col = "diff", output_dir)
  create_manhattan_plot(dml_dt, output_dir)
  create_qq_plot(dml_dt, output_dir)
  create_genomic_context_visualization(dml_dt, diff_col = "diff", output_dir)
  flog.info("Analysis complete", name = "dss_logger")
  return(dml_dt)
}

# perform analysis with provided parameters
flog.info(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f",
  delta, p.threshold, fdr.threshold
), name = "dss_logger")
result <- perform_dml_analysis(
  combined_bsseq,
  base_dir,
  delta = delta,
  smoothing = smoothing,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  min_coverage = 10,
  window_size = 500,
  cl = cl
)

stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:", name = "dss_logger")
flog.info(sprintf("Total DMLs found: %d", nrow(result)), name = "dss_logger")
flog.info(sprintf("Hypomethylated DMLs in tumor: %d", sum(result$hypo_in_tumour)), name = "dss_logger")
flog.info(sprintf("Hypermethylated DMLs in tumor: %d", sum(!result$hypo_in_tumour)), name = "dss_logger")
mean_diff <- mean(result$diff.meth)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff), name = "dss_logger")
median_diff <- median(result$diff.meth)
flog.info(sprintf("Median methylation difference: %.4f", median_diff), name = "dss_logger")
end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"), name = "dss_logger")
options(warn = 0)
