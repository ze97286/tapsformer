# git pull;clear;Rscript dss_differential_methylation.r 0.4 0.05 0.01 raw
# git pull;clear;Rscript dss_differential_methylation.r 0.4 0.05 0.01 raw_with_liver

start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript dml_analysis.R <delta> <p.threshold> <fdr.threshold> <suffix>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
suffix <- args[4]

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
flog.info("This is an INFO message.", name = "myLogger")
flog.warn("This is a WARN message.", name = "myLogger")
flog.error("This is an ERROR message.", name = "myLogger")

flog.info("Starting DSS DML analysis", name = "dss_logger")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

# generate a plot of top DMLs and their methylation levels in each sample
plot_top_DMLs <- function(top_hypo_dmls, combined_bsseq, output_dir) {
  methylation_data <- getMeth(combined_bsseq, type = "raw")
  plot_data <- data.table(chr = top_hypo_dmls$chr, pos = top_hypo_dmls$pos)
  for (i in 1:nrow(top_hypo_dmls)) {
    dml <- top_hypo_dmls[i, ]
    meth_levels <- methylation_data[which(seqnames(combined_bsseq) == dml$chr &
      start(combined_bsseq) == dml$pos), ]
    plot_data[i, (paste0("Sample_", 1:ncol(meth_levels))) := as.list(meth_levels)]
  }
  plot_data <- melt(plot_data,
    id.vars = c("chr", "pos"),
    variable.name = "Sample", value.name = "MethylationLevel"
  )
  plot_func <- function() {
    ggplot(plot_data, aes(x = Sample, y = MethylationLevel)) +
      geom_point() +
      facet_wrap(~ chr + pos, scales = "free_y") +
      theme_bw() +
      labs(
        title = "Methylation Levels Across Samples for Each DML",
        x = "Sample", y = "Methylation Level"
      ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  output_filename <- file.path(output_dir, "dml_methylation_plot.svg")
  safe_plot(output_filename, plot_func, width = 10, height = 8)
  return(output_filename)
}

# this is the core function here, doing the DML analysis + FDR correction, choosing hypomethylated DMLs and
# saving the output and visualisations.
perform_dml_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, cl) {
  output_dir <- file.path(base_dir, sprintf(
    "dml_delta_%.2f_p_%.4f_fdr_%.2f",
    delta, p.threshold, fdr.threshold
  ))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  flog.info("Performing DML test", name = "dss_logger")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = TRUE)

  flog.info("Calling DMLs", name = "dss_logger")
  dmls <- callDML(dml_test,
    delta = delta,
    p.threshold = p.threshold
  )

  dml_dt <- as.data.table(dmls)

  # identify hypomethylated regions in the tumour and check significance
  dml_dt[, `:=`(
    hypo_in_tumour = diff < 0,
    significant_after_fdr = fdr < fdr.threshold
  )]

  # Select top hypomethylated DMLs
  top_hypo_dmls <- dml_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmls, -stat)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmls, "stat", output_dir)

  # Tag differential methylation strength by quantile of areaStat
  top_hypo_dmls[, hypomethylation_strength := case_when(
    stat < thresholds$very_strong ~ "Very Strong",
    stat < thresholds$strong ~ "Strong",
    stat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # Save output
  fwrite(top_hypo_dmls[, .(chr, start, end, stat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmls.bed"),
    sep = "\t"
  )

  flog.info("Creating visualizations", name = "dss_logger")
  strongest_dmls <- tail(top_hypo_dmls[order(top_hypo_dmls$areaStat), ], 12)
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
result <- perform_dml_analysis(combined_bsseq, base_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
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
