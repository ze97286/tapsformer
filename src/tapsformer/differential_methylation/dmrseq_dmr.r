start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: Rscript dmrseq_dmr.r <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
min.CpG <- as.numeric(args[4])
min.len <- as.numeric(args[5])
dis.merge <- as.numeric(args[6])
smoothing_arg <- tolower(args[7]) # Convert to lowercase for consistency
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
print("Starting dmrseq analysis")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

# generate a plot of top DMRs and their methylation levels in each sample
plot_top_DMRs <- function(top_hypo_dmrs, combined_bsseq, output_dir, n = 20, ext = 0, prefix = "dss_") {
  dmr_plot_dir <- file.path(output_dir, prefix, "strongest_hypomethylated_dmr_plots")
  dir.create(dmr_plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Sort in descending order of areaStat and take top n
  strongest_dmrs <- head(top_hypo_dmrs[order(-areaStat)], n)

  for (i in 1:nrow(strongest_dmrs)) {
    dmr <- strongest_dmrs[i, ]
    print(sprintf("Processing DMR %d: chr%s:%d-%d", i, dmr$chr, dmr$start, dmr$end))
    filename <- file.path(dmr_plot_dir, sprintf(
      "DMR_%d_chr%s_%d-%d.svg",
      i, dmr$chr, dmr$start, dmr$end
    ))
    tryCatch(
      {
        plot_single_dmr(filename, dmr, combined_bsseq, i, ext)
      },
      error = function(e) {
        print(sprintf("Error plotting DMR %d: %s", i, conditionMessage(e)))
      }
    )
  }
  print(sprintf("Completed plotting %d strongest hypomethylated DMRs", n))
}

create_visualisations <- function(top_hypo_dmrs, combined_bsseq, output_dir, prefix, n) {
  plot_top_DMRs(top_hypo_dmrs, combined_bsseq, output_dir, n = 10, prefix = prefix)
  create_volcano_plot(top_hypo_dmrs, diff_col = "diff.Methy", pval_col = "pval", output_dir, prefix = prefix)
  create_methylation_diff_plot(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
  create_chromosome_coverage_plot(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
  create_dmr_length_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_manhattan_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_qq_plot(top_hypo_dmrs, output_dir, prefix = prefix)
  create_genomic_context_visualization(top_hypo_dmrs, diff_col = "diff.Methy", output_dir, prefix = prefix)
}


perform_dmrseq_analysis <- function(combined_bsseq, output_dir,
                                    testCovariate = "condition",
                                    cutoff = 0.05,
                                    beta_threshold = -0.4, # For hypomethylation
                                    minNumRegion = 5,
                                    maxGap = 1000
                                    ) {
  print("starting dmrseq analysis")
  sampleNames(combined_bsseq) <- gsub("tumour_", "tumour_", sampleNames(combined_bsseq))

  # Create condition vector
  condition <- factor(ifelse(grepl("tumour_", sampleNames(combined_bsseq)), "tumour", "control"))

  # Run DMRseq
  dmrs <- dmrseq(
    bs = combined_bsseq,
    cutoff = cutoff,
    testCovariate = testCovariate,
    minNumRegion = minNumRegion,
    maxGap = maxGap,
  )

  print("finished initial dmrseq analysis")

  # Convert to data.table and filter for significant hypomethylated DMRs
  dmr_dt <- as.data.table(dmrs)
  dmr_dt <- dmr_dt[pval <= p_threshold & beta < beta_threshold, .(
    chr = seqnames,
    start = start,
    end = end,
    diff.Methy = beta,
    areaStat = stat * width, # Using stat * width as a proxy for area
    pval = pval,
    fdr = qval,
    length = width,
    numCG = numCpGs
  )]

  # Add hypomethylation flag (all should be TRUE at this point)
  dmr_dt[, hypo_in_tumour := TRUE]

  print(sprintf("dmrseq Identified %d high-confidence hypomethylated DMRs", nrow(dmr_dt)))

  # Save as BED file
  bed_file <- file.path(output_dir, "dmrseq_high_confidence_hypomethylated_dmrs.bed")
  fwrite(dmr_dt[, .(chr, start, end, diff.Methy, areaStat, pval, fdr)],
    file = bed_file, sep = "\t", col.names = FALSE
  )
  saveRDS(dmr_dt, file.path(output_dir, "dmrseq_hypomethylated_dmrs.rds"))

  print(sprintf("DMRSeq analysis identified %d significant hypomethylated DMRs", nrow(dmr_dt)))
  print("dmrseq visualisation")
  create_visualisations(top_hypo_dmrs, combined_bsseq, output_dir, "dmrseq", n = 10)
  print("finished dmrseq")
}


# perform analysis with provided parameters
print(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
))
smoothing_string <- ifelse(smoothing, "smooth", "unsmooth")
output_dir <- file.path(base_dir, sprintf(
  "dmr_delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d_%s",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, smoothing_string
))

print("running dmrseq analysis")
dmrseq_results <- perform_dmrseq_analysis(
  combined_bsseq,
  output_dir,
  beta_threshold = -delta,
  minNumRegion = min.CpG,
  maxGap = dis.merge,
)

stopCluster(cl)

end_time <- Sys.time()
print(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
