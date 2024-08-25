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
flog.logger("dss_logger", appender.file(log_dir))
flog.threshold(INFO)
flog.threshold(WARN)
flog.threshold(ERROR)
print("Starting DSS DMR analysis")

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

# this is the core function here, doing the DMR analysis + FDR correction, choosing hypomethylated DMRs and
# saving the output and visualisations.
perform_dmr_analysis <- function(combined_bsseq, base_dir, output_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, smoothing, cl) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print("Performing DMR analysis with DML pre-filtering")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)

  # Perform DMLtest
  print("Performing DMLtest")
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = smoothing)

  # Convert DMLtest results to data.table for filtering
  dml_dt <- as.data.table(dml_test)

  # Apply filters similar to DML analysis
  z_score <- qnorm(0.975) # Two-tailed 95% CI
  dml_dt[, `:=`(
    hypo_in_tumour = diff < 0,
    significant_after_fdr = p.adjust(pval, method = "BH") < fdr.threshold,
    mean_methylation_diff = abs(diff),
    ci_excludes_zero = sign(diff - (z_score * diff.se)) == sign(diff + (z_score * diff.se))
  )]

  # Filter DMLs
  filtered_dmls <- dml_dt[
    hypo_in_tumour == TRUE &
      significant_after_fdr == TRUE &
      mean_methylation_diff >= delta &
      ci_excludes_zero == TRUE
  ]

  print(sprintf("Filtered DMLs: %d out of %d", nrow(filtered_dmls), nrow(dml_dt)))

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

  # Apply stringent filters to DMRs
  dmr_dt[, `:=`(
    hypo_in_tumour = diff.Methy < 0,
    significant_after_fdr = fdr < fdr.threshold,
    high_areaStat = areaStat > quantile(areaStat, 0.75) # Top 25% by areaStat
  )]

  # Select top hypomethylated DMRs
  top_hypo_dmrs <- dmr_dt[hypo_in_tumour == TRUE &
    significant_after_fdr == TRUE &
    high_areaStat == TRUE]

  setorder(top_hypo_dmrs, -areaStat)

  print(sprintf("Identified %d high-confidence hypomethylated DMRs", nrow(top_hypo_dmrs)))

  if (nrow(top_hypo_dmrs) == 0) {
    print("No high-confidence hypomethylated DMRs found. Stopping analysis.")
    return(NULL)
  }

  # Categorize DMR strength
  thresholds <- analyze_areastat_thresholds(top_hypo_dmrs, "areaStat", output_dir)
  top_hypo_dmrs[, hypomethylation_strength := case_when(
    areaStat >= thresholds$very_strong ~ "Very Strong",
    areaStat >= thresholds$strong ~ "Strong",
    areaStat >= thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

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
  create_visualisations(top_hypo_dmrs, combined_bsseq, output_dir, "dss_", n = 10)
  print("DSS Analysis complete")
  return(dmr_dt)
}


perform_dmrseq_analysis <- function(combined_bsseq, output_dir,
                                    testCovariate = "condition",
                                    p_threshold = 0.05,
                                    beta_threshold = -0.4, # For hypomethylation
                                    minNumRegion = 5,
                                    minCpGs = 3,
                                    maxGap = 1000,
                                    maxPerms = 10) {
  print("starting dmrseq analysis")
  sampleNames(combined_bsseq) <- gsub("tumour_", "tumour_", sampleNames(combined_bsseq))

  # Create condition vector
  condition <- factor(ifelse(grepl("tumour_", sampleNames(combined_bsseq)), "tumour", "control"))

  # Run DMRseq
  dmrs <- dmrseq(
    bs = combined_bsseq,
    cutoff = p_threshold,
    testCovariate = testCovariate,
    minNumRegion = minNumRegion,
    minCpGs = minCpGs,
    maxGap = maxGap,
    maxPerms = maxPerms
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

perform_dmrcate_analysis <- function(combined_bsseq, output_dir, delta, p.threshold, fdr.threshold, lambda, C = 2) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print("Performing DMR analysis using DMRcate")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)

  # Design matrix for DMRcate analysis
  design <- model.matrix(~ factor(c(rep("tumour", length(group1)), rep("control", length(group2)))))
  colnames(design) <- c("Intercept", "TumourVsControl")

  # Annotate CpGs with the design matrix
  print("Annotating CpGs for DMRcate")
  myAnnotation <- cpg.annotate(
    object = combined_bsseq,
    datatype = "BS",
    what = "perRegion", # "perBase" or "perRegion" depending on your analysis
    analysis.type = "differential",
    design = design,
    contrasts = TRUE,
    cont.matrix = "TumourVsControl",
    coef = 2
  )

  # Perform DMR analysis using DMRcate
  print("Calling DMRs with DMRcate")
  dmrcoutput <- dmrcate(
    object = myAnnotation,
    lambda = lambda,
    C = C,
    p.adjust.method = "BH",
    p.threshold = p.threshold
  )

  # Extract and filter the DMRs
  print("Extracting and filtering DMRs")
  dmrs <- extractRanges(dmrcoutput, delta = delta, minlen = 3) # Adjust `minlen` as necessary

  # Filter for hypomethylated regions
  print("Filtering for hypomethylated DMRs")
  dmrs_hypo <- dmrs[dmrs$stat < 0] # Retain only hypomethylated regions

  # Save the hypomethylated DMRs to an RDS file
  saveRDS(dmrs_hypo, file.path(output_dir, "dmrcate_hypomethylated_dmrs.rds"))

  # Create BED format dataframe for hypomethylated DMRs
  dmrs_hypo_bed <- data.frame(
    chr = as.character(seqnames(dmrs_hypo)),
    start = start(dmrs_hypo),
    end = end(dmrs_hypo),
    name = paste("HypoDMR", seq_along(dmrs_hypo), sep = "_"),
    score = round(dmrs_hypo$stat, 2), # Score indicates level of hypomethylation
    strand = rep(".", length(dmrs_hypo))
  )

  # Save BED file
  bed_file <- file.path(output_dir, "dmrcate_hypomethylated_dmrs.bed")
  write.table(dmrs_hypo_bed, file = bed_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  print(sprintf("Saved %d hypomethylated DMRs to %s", nrow(dmrs_hypo_bed), bed_file))

  # Visualization of hypomethylated DMRs
  print("Creating visualizations for hypomethylated DMRs")
  pdf(file.path(output_dir, "dmrcate_top10_hypomethylated_dmrs.pdf"))
  DMR.plot(ranges = dmrs_hypo, dmr_output = dmrcoutput, bsseq = combined_bsseq, numRegions = 10, genome = "hg38")
  dev.off()

  print("DMRcate Analysis complete")
  return(dmrs_hypo)
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

print("running bump hunter analysis")
dmrcate_results <- perform_dmrcate_analysis(
  combined_bsseq,
  output_dir,
  delta = delta,
  fdr.threshold = fdr.threshold,
  p.threshold = p.threshold,
  C = 2,
  lambda = 500,
)

print("running dmrseq analysis")
dmrseq_results <- perform_dmrseq_analysis(
  combined_bsseq,
  output_dir,
  p_threshold = p.threshold,
  beta_threshold = -delta,
  minCpGs = min.CpG,
)

print("running dss analysis")
result <- perform_dmr_analysis(combined_bsseq, base_dir, output_dir,
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

end_time <- Sys.time()
print(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
