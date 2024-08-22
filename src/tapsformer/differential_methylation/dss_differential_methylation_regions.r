# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw
# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw_with_liver

start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript dmr_analysis.R <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
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
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_%s_dmr_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", suffix, delta, p.threshold, fdr.threshold))

# Set up logging
flog.appender(appender.file(log_dir))
flog.info("Starting DSS DMR analysis")

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

    flog.info(sprintf("Processing DMR %d: chr%s:%d-%d", i, dmr$chr, dmr$start, dmr$end))

    filename <- file.path(dmr_plot_dir, sprintf(
      "DMR_%d_chr%s_%d-%d.svg",
      i, dmr$chr, dmr$start, dmr$end
    ))

    safe_plot(filename, function() {
      tryCatch(
        {
          par(mar = c(5, 4, 4, 2) + 0.1)
          showOneDMR(dmr, combined_bsseq, ext = ext)
          title(
            main = sprintf(
              "DMR %d: %s:%d-%d\nStrength: %s, areaStat: %.2f",
              i, dmr$chr, dmr$start, dmr$end, dmr$hypomethylation_strength, dmr$areaStat
            ),
            cex.main = 0.8
          )
        },
        error = function(e) {
          flog.error(sprintf("Error plotting DMR %d: %s", i, conditionMessage(e)))
          plot(1, type = "n", xlab = "", ylab = "", main = sprintf("Error plotting DMR %d", i))
          text(1, 1, labels = conditionMessage(e), cex = 0.8, col = "red")
        }
      )
    }, flog, width = 14, height = 12)
  }

  flog.info(sprintf("Completed plotting %d strongest hypomethylated DMRs", n))
}

# this is the core function here, doing the DMR analysis + FDR correction, choosing hypomethylated DMRs and
# saving the output and visualisations.
perform_dmr_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, cl) {
  output_dir <- file.path(base_dir, sprintf(
    "dmr_delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d",
    delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
  ))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  flog.info("Performing DML test")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = TRUE)

  flog.info("Calling DMRs")
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
  thresholds <- analyze_areastat_thresholds(top_hypo_dmrs, output_dir)

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
  flog.info("Creating visualizations")
  plot_top_DMRs(top_hypo_dmrs, combined_bsseq, output_dir, n = 50)
  create_volcano_plot(dmr_dt, output_dir, flog)
  create_methylation_diff_plot(dmr_dt, output_dir, flog)
  create_dmr_length_plot(dmr_dt, output_dir, flog)
  create_manhattan_plot(dmr_dt, output_dir, flog)
  create_qq_plot(dmr_dt, output_dir, flog)
  create_genomic_context_visualization(dmr_dt, output_dir, flog)

  flog.info("Analysis complete")
  return(dmr_dt)
}

# perform analysis with provided parameters
flog.info(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
))
result <- perform_dmr_analysis(combined_bsseq, base_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  min.CpG = min.CpG,
  min.len = min.len,
  dis.merge = dis.merge,
  cl = cl
)

stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:")
flog.info(sprintf("Total DMRs found: %d", nrow(result)))
flog.info(sprintf("Hypomethylated DMRs in tumor: %d", sum(result$hypo_in_tumour)))
flog.info(sprintf("Hypermethylated DMRs in tumor: %d", sum(!result$hypo_in_tumour)))
mean_diff <- mean(result$diff.Methy)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff))
median_diff <- median(result$diff.Methy)
flog.info(sprintf("Median methylation difference: %.4f", median_diff))
end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
