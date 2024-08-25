start_time <- Sys.time()

source("dss_common.r")

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: Rscript dmrcate_dmr.R <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
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

print("Starting dmrcate analysis")

# data loading
combined_bsseq <- load_and_combine_bsseq(base_dir, "tumour", "control")
gc()

perform_dmrcate_analysis <- function(combined_bsseq, output_dir, delta, lambda, C = 2, fdr_threshold = 0.05) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print("Performing DMR analysis using DMRcate")

  # Directly use bsseq::getMeth
  meth_data <- bsseq::getMeth(combined_bsseq, type = "raw")
  cov_data <- bsseq::getCoverage(combined_bsseq)

  # Ensure that the sample names in combined_bsseq are correctly labeled
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)

  # Design matrix for sequencing.annotate
  design <- model.matrix(~ 0 + factor(c(rep("tumour", length(group1)), rep("control", length(group2)))))
  colnames(design) <- c("tumour", "control")

  # Create the contrast matrix
  cont.matrix <- makeContrasts(TumourVsControl = tumour - control, levels = design)

  # Create GenomicRatioSet
  gr <- GenomicRanges::granges(combined_bsseq)
  grs <- GenomicRatioSet(
    gr = gr, beta = meth_data, M = cov_data,
    colData = SummarizedExperiment::colData(combined_bsseq)
  )

  # Run sequencing.annotate
  myAnnotation <- sequencing.annotate(
    grs,
    methdesign = design,
    contrasts = TRUE,
    cont.matrix = cont.matrix,
    coef = "TumourVsControl",
    fdr = fdr_threshold
  )

  # Perform DMR analysis using DMRcate
  print("Calling DMRs with DMRcate")
  dmrcoutput <- dmrcate(
    object = myAnnotation,
    lambda = lambda,
    C = C,
    min.cpgs = 3
  )

  # Extract and filter the DMRs
  print("Extracting and filtering DMRs")
  dmrs <- extractRanges(dmrcoutput, delta = delta)
  # Filter for hypomethylated regions
  print("Filtering for hypomethylated DMRs")
  dmrs_hypo <- dmrs[dmrs$stat < 0]

  # Save the hypomethylated DMRs to an RDS file
  saveRDS(dmrs_hypo, file.path(output_dir, "dmrcate_hypomethylated_dmrs.rds"))

  # Create BED format dataframe for hypomethylated DMRs
  dmrs_hypo_bed <- data.frame(
    chr = as.character(seqnames(dmrs_hypo)),
    start = start(dmrs_hypo),
    end = end(dmrs_hypo),
    name = paste("HypoDMR", seq_along(dmrs_hypo), sep = "_"),
    score = round(dmrs_hypo$stat, 2),
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

print("running dmrcate analysis")
perform_dmrcate_analysis(
  combined_bsseq, # Your BSseq object
  output_dir = output_dir,
  delta = delta,
  lambda = 500, # Smoothing parameter
  fdr_threshold = fdr.threshold,
  C = 2 # Clustering parameter
)

stopCluster(cl)

end_time <- Sys.time()
print(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
