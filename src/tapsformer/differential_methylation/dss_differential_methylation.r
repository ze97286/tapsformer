args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: Rscript dmr_analysis.R <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
}

delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
min.CpG <- as.numeric(args[4])
min.len <- as.numeric(args[5])
dis.merge <- as.numeric(args[6])

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
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite")

# Install and load packages
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)

# Print loaded package versions
sessionInfo()

base_dir <- "/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg"
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_dmr_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", delta, p.threshold, fdr.threshold))

# Set up logging
flog.appender(appender.file(log_dir))
flog.info("Starting DSS differential methylation analysis script")

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

# Load the data
tumour_rds_path <- file.path(base_dir, "tumour_data.rds")
control_rds_path <- file.path(base_dir, "control_data.rds")
tumour_data <- readRDS(tumour_rds_path)
control_data <- readRDS(control_rds_path)

if (nrow(tumour_data) == 0 || nrow(control_data) == 0) {
  flog.error("No data available after preprocessing. Check your input files and filtering criteria.")
  stop("No data available after preprocessing")
}

# Create BSseq objects
flog.info("Creating BSseq objects")
tryCatch({
  tumour_bsseq <- BSseq(chr = tumour_data$chr,  
                       pos = tumour_data$pos, 
                       M = as.matrix(tumour_data$X), 
                       Cov = as.matrix(tumour_data$N), 
                       sampleNames = "tumour")
  
  control_bsseq <- BSseq(chr = control_data$chr, 
                         pos = control_data$pos, 
                         M = as.matrix(control_data$X), 
                         Cov = as.matrix(control_data$N), 
                         sampleNames = "Control")
}, error = function(e) {
  flog.error("Error in creating BSseq objects:")
  flog.error(conditionMessage(e))
  flog.error("tumour data structure:")
  flog.error(capture.output(str(tumour_data)))
  flog.error("Control data structure:")
  flog.error(capture.output(str(control_data)))
  stop("Error in creating BSseq objects")
})

# Combine datasets
flog.info("Combining datasets")
combined_bsseq <- bsseq::combine(tumour_bsseq, control_bsseq)
saveRDS(combined_bsseq, file.path(base_dir, "combined_bsseq.rds"))
gc()

perform_dmr_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge) {
  # Define output directory based on parameters
  output_dir <- file.path(base_dir, sprintf("delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d", 
                                            delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Perform DML test
  flog.info("Performing DML test")
  flog.info(sprintf("Starting analysis with delta = %.2f, p.threshold = %.4f", delta, p.threshold))

  dml_test <- DMLtest(combined_bsseq, group1 = c("tumour"), group2 = c("Control"), smoothing = TRUE)

  # Call DMRs with new parameters
  flog.info(paste("Calling DMRs with delta =", delta, ", p.threshold =", p.threshold, 
                  ", min.CpG =", min.CpG, ", min.len =", min.len, ", dis.merge =", dis.merge))
  dmrs <- callDMR(dml_test, delta = delta, p.threshold = p.threshold, 
                  minlen = min.len, minCG = min.CpG, dis.merge = dis.merge)

  # Convert DMRs to data.table
  dmr_dt <- as.data.table(dmrs)

  # Log the columns in dmr_dt
  flog.info(paste("Columns in dmr_dt after callDMR:", paste(names(dmr_dt), collapse = ", ")))

  # Extract p-values from DML test results for the CpGs in our DMRs
  dml_results <- as.data.table(dml_test)
  setkey(dml_results, chr, pos)

  # Function to get minimum p-value for each DMR
  get_min_pval <- function(chr, start, end) {
    pvals <- dml_results[.(chr, start:end), on = .(chr, pos), nomatch = 0]$pval
    if (length(pvals) == 0) return(NA)
    min(pvals, na.rm = TRUE)
  }

  # Add minimum p-value for each DMR
  dmr_dt[, pval := mapply(get_min_pval, chr, start, end)]

  # Log the columns in dmr_dt again
  flog.info(paste("Columns in dmr_dt after adding pval:", paste(names(dmr_dt), collapse = ", ")))

  # Check if 'pval' column exists and log its summary
  if ("pval" %in% names(dmr_dt)) {
    flog.info(paste("Summary of pval column:", 
                    paste(capture.output(summary(dmr_dt$pval)), collapse = "\n")))
  } else {
    flog.error("'pval' column not found in dmr_dt after attempting to add it")
  }

  # Calculate FDR
  tryCatch({
    dmr_dt[, fdr := p.adjust(pval, method = "BH")]
    flog.info("FDR calculation successful")
  }, error = function(e) {
    flog.error(paste("Error in FDR calculation:", conditionMessage(e)))
  })
  
  saveRDS(dmr_dt, file = file.path(output_dir, "dmr_results.rds"))
  flog.info("DMR results saved")

  # Save DML test results to disk
  saveRDS(dml_test, file = file.path(output_dir, "dml_results.rds"))
  flog.info("DML test results saved")

  flog.info(paste("Columns in dmr_dt:", paste(names(dmr_dt), collapse = ", ")))

  # Flag hypo-methylated regions in tumour and apply FDR threshold
  dmr_dt[, `:=`(
    hypo_in_tumour = diff.Methy < 0,  # Lower methylation in tumor
    significant_after_fdr = fdr < fdr.threshold
  )]

  # Select top 100 hypomethylated DMRs after FDR correction
  top_hypo_dmrs <- dmr_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmrs, -areaStat)
  
  # Log the number of significant hypomethylated DMRs
  flog.info(paste("Number of significant hypomethylated DMRs:", nrow(top_hypo_dmrs)))

  # Write results to BED file
  bed_file <- file.path(output_dir, "hypomethylated_dmrs.bed")
  flog.info("Writing results to top 100 BED file")
  fwrite(top_hypo_dmrs[, .(chr, start, end, areaStat, pval, fdr)], bed_file, sep = "\t")

  # Visualizations
  flog.info("Creating visualizations")

  # 1. Volcano plot
  if ("pval" %in% names(dmr_dt)) {
    flog.info("Creating volcano plot")
    safe_plot(
      file.path(output_dir, "volcano_plot.svg"),
      function() {
        plot <- ggplot(dmr_dt, aes(x = diff.Methy, y = -log10(pval))) +
                geom_point(aes(color = hypo_in_tumour)) +
                scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), 
                                  labels = c("TRUE" = "Hypermethylated in Tumor", "FALSE" = "Hypomethylated in Tumor")) +
                labs(title = "Volcano plot of DMRs",
                    x = "Methylation Difference (Positive = Hypomethylation in Tumor)",
                    y = "-log10(p-value)",
                    color = "Methylation State") +
                theme_minimal()
        print(plot)
      }
    )
  } else {
    flog.warn("Skipping volcano plot due to missing 'pval' column")
  }

  # 2. Methylation difference distribution
  flog.info("Creating methylation difference distribution plot")
  safe_plot(
    file.path(output_dir, "methylation_difference_distribution.svg"),
    function() {
      plot <- ggplot(dmr_dt, aes(x = diff.Methy)) +
        geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
        labs(title = "Distribution of Methylation Differences",
             x = "Methylation Difference",
             y = "Count") +
        theme_minimal()
      print(plot)
    }
  )

  # 3. DMR length distribution
  flog.info("Creating DMR length distribution plot")
  safe_plot(
    file.path(output_dir, "dmr_length_distribution.svg"),
    function() {
      plot <- ggplot(dmr_dt, aes(x = length)) +
        geom_histogram(binwidth = 50, fill = "lightgreen", color = "black") +
        labs(title = "Distribution of DMR Lengths",
             x = "DMR Length (bp)",
             y = "Count") +
        theme_minimal()
      print(plot)
    }
  )

  # Create a data frame with chromosome sizes
  chr_sizes <- data.frame(
    chr = paste0("chr", c(1:22, "X", "Y")),
    size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
             145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
             101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
             156040895, 57227415)
  )

  # Add cumulative position
  chr_sizes$cumpos <- cumsum(as.numeric(chr_sizes$size))
  chr_sizes$pos <- chr_sizes$cumpos - chr_sizes$size/2

  # Add chromosome and cumulative position to DMR data
  dmr_dt$chr_num <- as.numeric(sub("chr", "", dmr_dt$chr))
  dmr_dt$cumpos <- dmr_dt$start + chr_sizes$cumpos[match(dmr_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dmr_dt$chr, chr_sizes$chr)]

  # Create the chromosome coverage plot
  p <- ggplot(dmr_dt, aes(x = cumpos, y = diff.Methy, color = diff.Methy > 0)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("blue", "red"), name = "Hypomethylated") +
    scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
    labs(title = "DMRs across chromosomes",
      x = "Chromosome", 
      y = "Methylation Difference (Positive = Hypomethylation in Tumor)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_plot(file.path(output_dir, "chromosome_coverage_plot.svg"), function() { print(p) })

  # Genomic Context Visualization
  dmr_gr <- GRanges(seqnames = dmr_dt$chr,
                    ranges = IRanges(start = dmr_dt$start, end = dmr_dt$end),
                    diff.Methy = dmr_dt$diff.Methy)

  # Get genomic features
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promoters <- promoters(txdb)
  genes <- genes(txdb)
  exons <- exons(txdb)
  introns <- gaps(exons)

  # Annotate DMRs
  dmr_annotation <- data.frame(
    DMR = seq_along(dmr_gr),
    Promoter = overlapsAny(dmr_gr, promoters),
    Gene = overlapsAny(dmr_gr, genes),
    Exon = overlapsAny(dmr_gr, exons),
    Intron = overlapsAny(dmr_gr, introns)
  )

  tryCatch({
    ah <- AnnotationHub()
    enhancers <- ah[["AH46978"]]  # GeneHancer enhancers for hg38
    dmr_annotation$Enhancer <- overlapsAny(dmr_gr, enhancers)
  }, error = function(e) {
    flog.warn("Failed to fetch enhancer data. Skipping enhancer annotation.")
    dmr_annotation$Enhancer <- FALSE
  })

  # Add annotation to original data
  dmr_dt$genomic_context <- apply(dmr_annotation[,-1], 1, function(x) {
    paste(names(x)[x], collapse = ";")
  })

  # Visualize genomic context distribution
  genomic_context_summary <- dmr_dt %>%
    tidyr::separate_rows(genomic_context, sep = ";") %>%
    dplyr::group_by(genomic_context) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(percentage = count / sum(count) * 100)

  p <- ggplot(genomic_context_summary, aes(x = "", y = percentage, fill = genomic_context)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = "Distribution of DMRs across Genomic Features")

  safe_plot(file.path(output_dir, "genomic_context_distribution.svg"), function() { print(p) })
  flog.info("Analysis complete")  
  return(dmr_dt)
}

# Perform analysis with provided parameters
flog.info(sprintf("Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d", 
                  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge))
result <- perform_dmr_analysis(combined_bsseq, base_dir, 
                               delta = delta, 
                               p.threshold = p.threshold, 
                               fdr.threshold = fdr.threshold,
                               min.CpG = min.CpG,
                               min.len = min.len,
                               dis.merge = dis.merge)
gc()

# Save the result
saveRDS(result, file.path(base_dir, sprintf("delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d", 
                                            delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge), "dmr_results.rds"))
flog.info("Analysis completed successfully")

end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))