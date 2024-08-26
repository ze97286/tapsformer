start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript dss_dmr.r <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge> <suffix>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
min.CpG <- as.numeric(args[4])
min.len <- as.numeric(args[5])
dis.merge <- as.numeric(args[6])
suffix <- args[7]

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

print("Starting dss dmr analysis")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

# generate a plot of top DMRs and their methylation levels in each sample
plot_top_DMRs <- function(top_hypo_dmrs, combined_bsseq, output_dir, n = 20, ext = 0, prefix = "dss") {
  dmr_plot_dir <- file.path(output_dir, prefix, "strongest_hypomethylated_dmr_plots")
  dir.create(dmr_plot_dir, showWarnings = FALSE, recursive = TRUE)

  strongest_dmrs <- head(top_hypo_dmrs[order(-top_hypo_dmrs$areaStat), ], n)

  for (i in 1:nrow(strongest_dmrs)) {
    dmr <- strongest_dmrs[i, ]
    print(sprintf("Processing DMR %d: chr%s:%d-%d", i, dmr$chr, dmr$start, dmr$end))
    filename <- file.path(dmr_plot_dir, sprintf(
      "DMR_%d_chr%s_%d-%d.svg",
      i, dmr$chr, dmr$start, dmr$end
    ))
    plot_single_dmr(filename, dmr, combined_bsseq, i, ext)
  }
  print(sprintf("Completed plotting %d strongest hypomethylated DMRs", n))
}

# Generate plots from output
create_visualisations <- function(top_hypo_dmrs, combined_bsseq, output_dir, prefix, n) {
  plot_top_DMRs(top_hypo_dmrs, combined_bsseq, output_dir, n = 100, prefix = prefix)
  create_volcano_plot(top_hypo_dmrs, diff_col = "diff.Methy", pval_col = "pval", output_dir, prefix = prefix)
  create_methylation_diff_plot(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
  create_chromosome_coverage_plot(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
  create_dmr_length_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_manhattan_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_qq_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_genomic_context_visualization(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
}

# Perform DMR analysis + FDR correction, choosing hypomethylated DMRs and saving the output and visualisations.
perform_dmr_analysis <- function(
    combined_bsseq, base_dir, output_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, cl, areaStat_percentile = 0.75,
    n_iterations = 5,
    subsample_fraction = 0.8, min_span = 200, max_span = 1000, min_cpgs = 20) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print("Performing DMR analysis with DML pre-filtering")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)

  if (length(group1) < 2 || length(group2) < 2) {
    stop("Each group must have at least 2 samples. Found ", length(group1), " tumour samples and ", length(group2), " control samples.")
  }
  print(sprintf("Analysis will be performed on %d tumour samples and %d control samples", length(group1), length(group2)))

  print("Performing quality control checks")

  # Check distribution of methylation values
  meth_rates <- bsseq::getMeth(combined_bsseq, type = "raw")
  meth_summary <- summary(as.vector(meth_rates))
  print("Summary of methylation rates across all samples:")
  print(meth_summary)

  # Check coverage across samples
  coverage <- bsseq::getCoverage(combined_bsseq)
  cov_summary <- summary(as.vector(coverage))
  print("Summary of coverage across all samples:")
  print(cov_summary)

  # run cross validation dml
  consistent_dmls <- cross_validate_dmls(combined_bsseq, group1, group2,
    n_iterations = n_iterations,
    subsample_fraction = subsample_fraction,
    delta = delta, p.threshold = p.threshold,
    fdr.threshold = fdr.threshold, smoothing = FALSE
  )
  saveRDS(consistent_dmls, file.path(output_dir, "consistent_dmls.rds"))

  print("Performing DMLtest")
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = TRUE)
  dml_dt <- as.data.table(dml_test)

  dml_dt[, `:=`(
    hypo_in_tumour = diff < 0,
    significant_after_fdr = p.adjust(pval, method = "BH") < fdr.threshold,
    consistent = paste(chr, pos) %in% consistent_dmls,
    mean_methylation_diff = abs(diff),
    ci_excludes_zero = sign(diff - (z_score * diff.se)) == sign(diff + (z_score * diff.se))
  )]

  # Filter DMLs
  # Stage 1: Hypomethylation
  hypo_dmls <- dml_dt[hypo_in_tumour == TRUE]
  print(sprintf("DMLs after hypomethylation filter: %d", nrow(hypo_dmls)))

  # Stage 2: FDR significance
  sig_hypo_dmls <- hypo_dmls[significant_after_fdr == TRUE]
  print(sprintf("DMLs after FDR significance filter: %d", nrow(sig_hypo_dmls)))

  diff_sig_hypo_dmls <- sig_hypo_dmls[mean_methylation_diff >= delta]
  print(sprintf("DMLs after methylation difference filter: %d", nrow(diff_sig_hypo_dmls)))

  # Stage 4: Confidence interval
  ci_diff_sig_hypo_dmls <- diff_sig_hypo_dmls[ci_excludes_zero == TRUE]
  print(sprintf("DMLs after confidence interval filter: %d", nrow(ci_diff_sig_hypo_dmls)))

  # Stage 5: Consistency (from cross-validation)
  consistent_ci_diff_sig_hypo_dmls <- ci_diff_sig_hypo_dmls[consistent == TRUE]
  print(sprintf("DMLs after consistency filter: %d", nrow(consistent_ci_diff_sig_hypo_dmls)))

  # Stage 6: Sliding window filter
  filtered_dmls <- sliding_window_filter(consistent_ci_diff_sig_hypo_dmls, window_size)
  print(sprintf("DMLs after sliding window filter: %d", nrow(top_hypo_dmls)))

  # Call DMRs using filtered DMLs
  print("Calling DMRs on filtered DMLs")
  dmrs <- callDMR(filtered_dmls, delta = delta, p.threshold = p.threshold, minlen = min.len, minCG = min.CpG, dis.merge = dis.merge)

  print(sprintf("Initial DMR analysis found %d DMRs", nrow(dmrs)))

  # Convert to data.table for further processing
  dmr_dt <- as.data.table(dmrs)

  # Calculate more stringent p-values for DMRs
  setkey(filtered_dmls, chr, pos)

  get_min_pval <- function(dmr_row) {
    pvals <- filtered_dmls[.(dmr_row[["chr"]], dmr_row[["start"]]:dmr_row[["end"]]), on = .(chr, pos), nomatch = 0]$pval
    if (length(pvals) == 0) {
      return(NA)
    }
    min(pvals, na.rm = TRUE)
  }

  dmr_dt[, pval := apply(.SD, 1, get_min_pval)]
  dmr_dt[, fdr := p.adjust(pval, method = "BH")]

  # Categorize DMR strength based on areaStat thresholds before filtering
  thresholds <- analyze_areastat_thresholds(dmr_dt, "areaStat", output_dir)
  dmr_dt[, hypomethylation_strength := case_when(
    areaStat >= thresholds$very_strong ~ "Very Strong",
    areaStat >= thresholds$strong ~ "Strong",
    areaStat >= thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # Filter based on hypomethylation, significant FDR, and high areaStat
  dmr_dt[, `:=`(
    hypo_in_tumour = diff.Methy < 0,
    significant_after_fdr = fdr < fdr.threshold,
    high_areaStat = areaStat > quantile(areaStat, areaStat_percentile)
  )]

  top_hypo_dmrs <- dmr_dt[hypo_in_tumour == TRUE &
    significant_after_fdr == TRUE &
    high_areaStat == TRUE]

  setorder(top_hypo_dmrs, -areaStat)

  print(sprintf("Identified %d high-confidence hypomethylated DMRs", nrow(top_hypo_dmrs)))

  if (nrow(top_hypo_dmrs) == 0) {
    print("No high-confidence hypomethylated DMRs found. Stopping analysis.")
    return(NULL)
  }

  # Save output
  fwrite(top_hypo_dmrs[, .(chr, start, end, areaStat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "dss_high_confidence_hypomethylated_dmrs.bed"),
    sep = "\t"
  )

  print(table(top_hypo_dmrs$hypomethylation_strength))

  saveRDS(top_hypo_dmrs, file.path(output_dir, "dss_top_hypo_dmrs.rds"))
  saveRDS(combined_bsseq, file.path(output_dir, "combined_bsseq.rds"))

  # Visualizations
  print("Creating visualizations")
  create_visualisations(top_hypo_dmrs, combined_bsseq, output_dir, "dss", n = 100)
  print("DSS dmr analysis complete")
  return(dmr_dt)
}

# perform analysis with provided parameters
print(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
))
output_dir <- file.path(base_dir, sprintf(
  "dmr_delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
))

print("running dss analysis")
result <- perform_dmr_analysis(combined_bsseq, base_dir, output_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  min.CpG = min.CpG,
  min.len = min.len,
  dis.merge = dis.merge,
  cl = cl
)

stopCluster(cl)

end_time <- Sys.time()
print(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
