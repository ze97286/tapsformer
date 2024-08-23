# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw
# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw_with_liver

start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: Rscript dmr_analysis.R <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
min.CpG <- as.numeric(args[4])
min.len <- as.numeric(args[5])
dis.merge <- as.numeric(args[6])
smoothing_arg <- tolower(args[7])  # Convert to lowercase for consistency
smoothing <- if (smoothing_arg == "true") TRUE else if (smoothing_arg == "false") FALSE else NA
if (is.na(smoothing)) {
  stop("Invalid smoothing argument. Please use 'TRUE' or 'FALSE'.")
}
suffix <- args[8]

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
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_%s_dmr_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", suffix, delta, p.threshold, fdr.threshold))

# Set up logging
flog.logger("dss_logger", appender.file(log_dir))
flog.threshold(INFO, name = "dss_logger")
flog.threshold(WARN, name = "dss_logger")
flog.threshold(ERROR, name = "dss_logger")
flog.info("Starting DSS DMR analysis", name = "dss_logger")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

# generate a plot of top DMRs and their methylation levels in each sample
plot_top_DMRs <- function(top_hypo_dmrs, combined_bsseq, output_dir, n = 50, ext = 0) {
  dmr_plot_dir <- file.path(output_dir, "strongest_hypomethylated_dmr_plots")
  dir.create(dmr_plot_dir, showWarnings = FALSE, recursive = TRUE)

  strongest_dmrs <- tail(top_hypo_dmrs[order(top_hypo_dmrs$areaStat), ], n)

  for (i in 1:nrow(strongest_dmrs)) {
    dmr <- strongest_dmrs[i, ]

    flog.info(sprintf("Processing DMR %d: chr%s:%d-%d", i, dmr$chr, dmr$start, dmr$end), name = "dss_logger")

    filename <- file.path(dmr_plot_dir, sprintf(
      "DMR_%d_chr%s_%d-%d.svg",
      i, dmr$chr, dmr$start, dmr$end
    ))

    plot_single_dmr(filename, dmr, combined_bsseq, i, ext)
  }

  flog.info(sprintf("Completed plotting %d strongest hypomethylated DMRs", n), name = "dss_logger")
}

# this is the core function here, doing the DMR analysis + FDR correction, choosing hypomethylated DMRs and
# saving the output and visualisations.
perform_dmr_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, smoothing, cl) {
  smoothing_string <- ifelse(smoothing, "smooth", "unsmooth")
  output_dir <- file.path(base_dir, sprintf(
    "dmr_delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d_%s",
    delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, smoothing_string
  ))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  flog.info("Performing DMR test", name = "dss_logger")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = smoothing, test = "BB")

  flog.info("Calling DMRs", name = "dss_logger")
  dmrs <- callDMR(
    dml_test,
    delta = delta,
    p.threshold = p.threshold,
    minlen = min.len,
    minCG = min.CpG,
    dis.merge = dis.merge
  )

  dmr_dt <- as.data.table(dmrs)
  dml_results <- as.data.table(dml_test)
  setkey(dml_results, chr, pos)

  get_min_pval <- function(dmr) {
    pvals <- dml_results[.(dmr$chr, dmr$start:dmr$end), on = .(chr, pos), nomatch = 0]$pval
    if (length(pvals) == 0) {
      return(NA)
    }
    min(pvals, na.rm = TRUE)
  }
  dmr_dt[, pval := parSapply(cl, split(dmr_dt, 1:nrow(dmr_dt)), get_min_pval)]
  dmr_dt[, fdr := p.adjust(pval, method = "BH")]
  dmr_dt[, `:=`(
    hypo_in_tumour = diff.Methy < 0,
    significant_after_fdr = fdr < fdr.threshold
  )]

  # select top hypomethylated DMLs
  top_hypo_dmrs <- dmr_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmrs, -areaStat)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmrs, "areaStat", output_dir)

  # tag differential methylation strength by quantile of areastat
  top_hypo_dmrs[, hypomethylation_strength := case_when(
    areaStat < thresholds$very_strong ~ "Very Strong",
    areaStat < thresholds$strong ~ "Strong",
    areaStat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # save output
  fwrite(top_hypo_dmrs[, .(chr, start, end, areaStat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmrs.bed"),
    sep = "\t"
  )

  # Visualizations
  flog.info("Creating visualizations", name = "dss_logger")
  plot_top_DMRs(top_hypo_dmrs, combined_bsseq, output_dir, n = 20)
  create_volcano_plot(dmr_dt, diff_col = "diff.Methy", pval_col = "pval", output_dir)
  create_methylation_diff_plot(dmr_dt, diff_col = "diff.Methy", output_dir)
  create_chromosome_coverage_plot(dmr_dt, diff_col = "diff.Methy", output_dir)
  create_dmr_length_plot(dmr_dt, output_dir)
  create_manhattan_plot(dmr_dt, output_dir)
  create_qq_plot(dmr_dt, output_dir)
  create_genomic_context_visualization(dmr_dt, diff_col = "diff.Methy", output_dir)
  flog.info("Analysis complete", name = "dss_logger")
  return(dmr_dt)
}

# perform analysis with provided parameters
flog.info(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
), name = "dss_logger")
result <- perform_dmr_analysis(combined_bsseq, base_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  min.CpG = min.CpG,
  min.len = min.len,
  dis.merge = dis.merge,
  smoothing = smoothing,
  cl = cl
)

stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:", name = "dss_logger")
flog.info(sprintf("Total DMRs found: %d", nrow(result)), name = "dss_logger")
flog.info(sprintf("Hypomethylated DMRs in tumor: %d", sum(result$hypo_in_tumour)), name = "dss_logger")
flog.info(sprintf("Hypermethylated DMRs in tumor: %d", sum(!result$hypo_in_tumour)), name = "dss_logger")
mean_diff <- mean(result$diff.Methy)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff), name = "dss_logger")
median_diff <- median(result$diff.Methy)
flog.info(sprintf("Median methylation difference: %.4f", median_diff), name = "dss_logger")
end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"), name = "dss_logger")
options(warn = 0)
