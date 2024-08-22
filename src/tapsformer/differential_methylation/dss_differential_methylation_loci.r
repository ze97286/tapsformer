# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw
# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw_with_liver

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop("Usage: Rscript dml_analysis.R <delta> <p.threshold> <fdr.threshold> <min.CpG> <min.len> <dis.merge>")
}

delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
min.CpG <- as.numeric(args[4])
min.len <- as.numeric(args[5])
dis.merge <- as.numeric(args[6])
suffix <- args[7]

start_time <- Sys.time()

install_and_load <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Function for Bioconductor packages
bioc_install_and_load <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, update = FALSE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Define packages
bioc_packages <- c(
  "DSS", "GenomicRanges", "bsseq", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub"
)
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite", "pheatmap", "reshape")

# Install and load packages
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)

# Print loaded package versions
sessionInfo()

# Set up parallel processing
library(parallel)
num_cores <- detectCores() - 1 # Use all but one core
cl <- makeCluster(num_cores)

# Load required packages on all cores
clusterEvalQ(cl, {
  library(DSS)
  library(data.table)
  library(GenomicRanges)
  library(bsseq)
})

options(warn = -1)
suppressMessages(library(DSS))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(bsseq))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(AnnotationHub))

base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg", suffix)
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_dml_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", delta, p.threshold, fdr.threshold))

# Set up logging
flog.appender(appender.file(log_dir))
flog.info("Starting optimized DSS differential methylation analysis script")

safe_plot <- function(filename, plot_func, width = 10, height = 8) {
  tryCatch(
    {
      svglite::svglite(filename, width = width, height = height)
      plot_func()
      dev.off()
      flog.info(paste("Plot saved as", filename))
    },
    error = function(e) {
      if (dev.cur() > 1) dev.off() # Close device if open
      flog.error(paste("Error creating plot:", filename, "-", conditionMessage(e)))
    }
  )
}

load_and_create_bsseq <- function(base_dir, prefix) {
  sample_files <- list.files(path = base_dir, pattern = paste0("^", prefix, "_.*\\.rds$"), full.names = TRUE)
  if (length(sample_files) == 0) {
    stop("Error: No RDS files found with the given prefix.")
  }

  # Create the BSseq objects for all samples
  bsseq_list <- lapply(sample_files, function(file_path) {
    sample_data <- readRDS(file_path)
    sample_name <- gsub("\\.rds$", "", basename(file_path))  # Keep the full name including the prefix
    BSseq(
      chr = sample_data$chr,
      pos = sample_data$pos,
      M = as.matrix(sample_data$X),
      Cov = as.matrix(sample_data$N),
      sampleNames = sample_name
    )
  })

  # Combine the list of BSseq objects into one BSseq object
  combined_bsseq <- do.call(combineList, bsseq_list)

  return(combined_bsseq)
}

# Load the data
tumour_bsseq <- load_and_create_bsseq(base_dir, "tumour")
control_bsseq <- load_and_create_bsseq(base_dir, "control")

# Combine tumour and control BSseq objects
combined_bsseq <- bsseq::combine(tumour_bsseq, control_bsseq)

saveRDS(combined_bsseq, file.path(base_dir, "combined_bsseq.rds"))
gc()

analyze_areastat_thresholds <- function(top_hypo_dmls, output_dir) {
  flog.info("Analyzing areaStat distribution and thresholds")

  # Calculate quantiles
  quantiles <- quantile(top_hypo_dmls$areaStat, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))

  flog.info("areaStat quantiles:")
  print(quantiles)

  # Suggest thresholds
  moderate_threshold <- quantiles["75%"]
  strong_threshold <- quantiles["90%"]
  very_strong_threshold <- quantiles["99%"]

  flog.info(sprintf("Suggested thresholds for hypomethylation strength:"))
  flog.info(sprintf("Moderate: < %.2f", moderate_threshold))
  flog.info(sprintf("Strong: < %.2f", strong_threshold))
  flog.info(sprintf("Very strong: < %.2f", very_strong_threshold))

  # Create histogram
  safe_plot(
    file.path(output_dir, "areastat_distribution.svg"),
    function() {
      p <- ggplot(top_hypo_dmls, aes(x = areaStat)) +
        geom_histogram(bins = 50, fill = "skyblue", color = "black") +
        geom_vline(
          xintercept = c(moderate_threshold, strong_threshold, very_strong_threshold),
          color = c("green", "orange", "red"), linetype = "dashed"
        ) +
        annotate("text",
          x = c(moderate_threshold, strong_threshold, very_strong_threshold),
          y = Inf, label = c("Moderate", "Strong", "Very Strong"),
          color = c("green", "orange", "red"), hjust = 1, vjust = 2, angle = 90
        ) +
        labs(
          title = "Distribution of areaStat Values",
          x = "areaStat",
          y = "Count"
        ) +
        theme_minimal()
      print(p)
    }
  )

  # Return thresholds
  return(list(
    moderate = moderate_threshold,
    strong = strong_threshold,
    very_strong = very_strong_threshold
  ))
}

perform_dml_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, cl) {
  output_dir <- file.path(base_dir, sprintf(
    "delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d",
    delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
  ))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Perform DML test
  flog.info("Performing DML test")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = TRUE)

  # Call DMLs
  flog.info("Calling DMLs")
  dmls <- callDML(dml_test,
    delta = delta,
    p.threshold = p.threshold,
    smoothing = TRUE,
    minCoverage = 5,
    stat = "fdr",
  )

  # Convert DMLs to data.table
  dml_dt <- as.data.table(dmls)

  # Extract p-values from DML test results for the CpGs in our DMLs
  dml_results <- as.data.table(dml_test)
  setkey(dml_results, chr, pos)

  # Since p-values are already corrected, you don't need to calculate the minimum p-value again.
  # The callDML already applies FDR, so we skip this step.

  # Flag hypo-methylated regions in the tumour and check significance
  dml_dt[, `:=`(
    hypo_in_tumour = diff.meth < 0,
    significant_after_fdr = fdr < fdr.threshold
  )]

  # Select top hypomethylated DMLs after FDR correction
  top_hypo_dmls <- dml_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmls, -areaStat)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmls, output_dir)

  # You can then use these thresholds in your analysis, e.g.:
  top_hypo_dmls[, hypomethylation_strength := case_when(
    areaStat < thresholds$very_strong ~ "Very Strong",
    areaStat < thresholds$strong ~ "Strong",
    areaStat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # Write results to BED file
  fwrite(top_hypo_dmls[, .(chr, start, end, areaStat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmls.bed"),
    sep = "\t"
  )

  # Visualizations
  flog.info("Creating visualizations")

  plot_top_DMLs <- function(top_hypo_dmls, combined_bsseq, output_dir) {
    # Extract methylation data for the DMLs
    methylation_data <- getMeth(combined_bsseq, type = "raw")

    # Prepare the data table
    plot_data <- data.table(chr = top_hypo_dmls$chr, pos = top_hypo_dmls$pos)

    # Add methylation levels from all samples
    for (i in 1:nrow(top_hypo_dmls)) {
      dml <- top_hypo_dmls[i, ]
      meth_levels <- methylation_data[which(seqnames(combined_bsseq) == dml$chr &
        start(combined_bsseq) == dml$pos), ]
      plot_data[i, (paste0("Sample_", 1:ncol(meth_levels))) := as.list(meth_levels)]
    }

    # Reshape the data for ggplot
    plot_data <- melt(plot_data,
      id.vars = c("chr", "pos"),
      variable.name = "Sample", value.name = "MethylationLevel"
    )

    # Define the plot function
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

    # Define the output filename
    output_filename <- file.path(output_dir, "dml_methylation_plot.svg")

    # Use safe_plot to save the plot to disk
    safe_plot(output_filename, plot_func, width = 10, height = 8)

    # Return the path to the saved plot
    return(output_filename)
  }
  strongest_dmls <- tail(top_hypo_dmls[order(top_hypo_dmls$areaStat), ], 12)
  plot_top_DMLs(strongest_dmls, combined_bsseq, output_dir)

  # 1. Volcano plot
  if ("pval" %in% names(dml_dt)) {
    flog.info("Creating volcano plot")
    safe_plot(
      file.path(output_dir, "volcano_plot.svg"),
      function() {
        plot <- ggplot(dml_dt, aes(x = diff.meth, y = -log10(pval))) +
          geom_point(aes(color = hypo_in_tumour)) +
          scale_color_manual(
            values = c("TRUE" = "blue", "FALSE" = "red"),
            labels = c("TRUE" = "Hypomethylated in Tumor", "FALSE" = "Hypermethylated in Tumor")
          ) +
          labs(
            title = "Volcano plot of DMLs",
            x = "Methylation Difference (Negative = Hypomethylation in Tumor)",
            y = "-log10(p-value)",
            color = "Methylation State"
          ) +
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
      plot <- ggplot(dml_dt, aes(x = diff.meth)) +
        geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
        labs(
          title = "Distribution of Methylation Differences",
          x = "Methylation Difference",
          y = "Count"
        ) +
        theme_minimal()
      print(plot)
    }
  )

  # 3. DML length distribution
  flog.info("Creating DML length distribution plot")
  safe_plot(
    file.path(output_dir, "dml_length_distribution.svg"),
    function() {
      plot <- ggplot(dml_dt, aes(x = length)) +
        geom_histogram(binwidth = 50, fill = "lightgreen", color = "black") +
        labs(
          title = "Distribution of DML Lengths",
          x = "DML Length (bp)",
          y = "Count"
        ) +
        theme_minimal()
      print(plot)
    }
  )

  # Create a data frame with chromosome sizes
  chr_sizes <- data.frame(
    chr = paste0("chr", c(1:22, "X", "Y")),
    size = c(
      248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
      145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
      101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
      156040895, 57227415
    )
  )

  # Add cumulative position
  chr_sizes$cumpos <- cumsum(as.numeric(chr_sizes$size))
  chr_sizes$pos <- chr_sizes$cumpos - chr_sizes$size / 2

  # Add chromosome and cumulative position to DML data
  dml_dt$chr_num <- as.numeric(sub("chr", "", dml_dt$chr))
  dml_dt$cumpos <- dml_dt$start + chr_sizes$cumpos[match(dml_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dml_dt$chr, chr_sizes$chr)]

  # Create the chromosome coverage plot
  p <- ggplot(dml_dt, aes(x = cumpos, y = diff.meth, color = diff.meth < 0)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(
      values = c("FALSE" = "red", "TRUE" = "blue"),
      name = "Methylation State",
      labels = c("FALSE" = "Hypermethylated", "TRUE" = "Hypomethylated")
    ) +
    scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
    labs(
      title = paste(
        "DMLs across chromosomes\n",
        "Hypomethylated:", sum(dml_dt$diff.meth < 0),
        "Hypermethylated:", sum(dml_dt$diff.meth >= 0)
      ),
      x = "Chromosome",
      y = "Methylation Difference (Negative = Hypomethylation in Tumor)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_plot(file.path(output_dir, "chromosome_coverage_plot.svg"), function() {
    print(p)
  })

  # Manhattan Plot
  create_manhattan_plot <- function(dml_dt, output_dir) {
    flog.info("Creating Manhattan plot")
    safe_plot(
      file.path(output_dir, "manhattan_plot.svg"),
      function() {
        plot <- ggplot(dml_dt, aes(x = cumpos, y = -log10(pval), color = factor(chr))) +
          geom_point(alpha = 0.8, size = 1) +
          scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
          scale_color_manual(values = rep(c("#1B9E77", "#D95F02"), 12)) +
          labs(
            title = "Manhattan Plot of DMLs",
            x = "Chromosome",
            y = "-log10(p-value)"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
          )
        print(plot)
      }
    )
  }

  # Q-Q Plot
  create_qq_plot <- function(dml_dt, output_dir) {
    flog.info("Creating Q-Q plot")
    safe_plot(
      file.path(output_dir, "qq_plot.svg"),
      function() {
        observed <- sort(-log10(dml_dt$pval))
        expected <- -log10(ppoints(length(observed)))
        qq_df <- data.frame(observed = observed, expected = expected)

        plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
          geom_point() +
          geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
          labs(
            title = "Q-Q Plot of Observed vs Expected p-values",
            x = "Expected -log10(p-value)",
            y = "Observed -log10(p-value)"
          ) +
          theme_minimal()
        print(plot)
      }
    )
  }

  # Heatmap of Top DMLs
  create_top_dml_heatmap <- function(combined_bsseq, dml_dt, output_dir, top_n = 100) {
    flog.info("Creating heatmap of top DMLs")
    safe_plot(
      file.path(output_dir, "top_dml_heatmap.png"), # Changed to PNG for testing
      function() {
        # Check if dml_dt is not empty
        if (nrow(dml_dt) == 0) {
          flog.error("dml_dt is empty. Cannot create heatmap.")
          return(NULL)
        }

        top_dmls <- head(dml_dt[order(-abs(areaStat))], n = top_n)

        # Check if top_dmls is not empty
        if (nrow(top_dmls) == 0) {
          flog.error("No top DMLs selected. Cannot create heatmap.")
          return(NULL)
        }

        dml_ranges <- GRanges(
          seqnames = top_dmls$chr,
          ranges = IRanges(start = top_dmls$start, end = top_dmls$end)
        )

        meth_mat <- getMeth(combined_bsseq, regions = dml_ranges, type = "smooth", what = "perRegion")

        # Check if meth_mat is not empty
        if (is.null(meth_mat) || all(is.na(meth_mat))) {
          flog.error("Methylation matrix is empty or all NA. Cannot create heatmap.")
          return(NULL)
        }

        # Print dimensions of meth_mat for debugging
        flog.info(sprintf("Methylation matrix dimensions: %d rows, %d columns", nrow(meth_mat), ncol(meth_mat)))

        pheatmap(meth_mat,
          scale = "none",
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          show_rownames = FALSE,
          main = paste("Methylation Levels of Top", top_n, "DMLs"),
          filename = file.path(output_dir, "top_dml_heatmap.png")
        ) # Changed to PNG
      }
    )
  }

  create_manhattan_plot(dml_dt, output_dir)
  create_qq_plot(dml_dt, output_dir)
  create_top_dml_heatmap(combined_bsseq, dml_dt, output_dir)

  # Genomic Context Visualization
  flog.info("Annotating DMLs with genomic context")
  dml_gr <- GRanges(
    seqnames = dml_dt$chr,
    ranges = IRanges(start = dml_dt$start, end = dml_dt$end),
    diff.meth = dml_dt$diff.meth
  )

  # Get genomic features
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  # Filter out non-standard chromosomes
  standard_chromosomes <- paste0("chr", c(1:22))
  dml_gr <- dml_gr[seqnames(dml_gr) %in% standard_chromosomes]

  promoters <- promoters(txdb)
  genes <- suppressWarnings(genes(txdb, single.strand.genes.only = TRUE))
  exons <- exons(txdb)
  introns <- gaps(exons)

  # Ensure all genomic features are on standard chromosomes
  promoters <- promoters[seqnames(promoters) %in% standard_chromosomes]
  genes <- genes[seqnames(genes) %in% standard_chromosomes]
  exons <- exons[seqnames(exons) %in% standard_chromosomes]
  introns <- introns[seqnames(introns) %in% standard_chromosomes]

  # Annotate DMLs
  dml_annotation <- data.frame(
    DML = seq_along(dml_gr),
    Promoter = overlapsAny(dml_gr, promoters),
    Gene = overlapsAny(dml_gr, genes),
    Exon = overlapsAny(dml_gr, exons),
    Intron = overlapsAny(dml_gr, introns)
  )

  tryCatch(
    {
      ah <- AnnotationHub()
      enhancers <- ah[["AH46978"]] # GeneHancer enhancers for hg38
      enhancers <- enhancers[seqnames(enhancers) %in% standard_chromosomes]
      dml_annotation$Enhancer <- overlapsAny(dml_gr, enhancers)
    },
    error = function(e) {
      flog.warn("Failed to fetch enhancer data. Skipping enhancer annotation.")
      dml_annotation$Enhancer <- FALSE
    }
  )

  # Add annotation to original data
  dml_dt$genomic_context <- apply(dml_annotation[, -1], 1, function(x) {
    paste(names(x)[x], collapse = ";")
  })

  # Visualize genomic context distribution
  genomic_context_summary <- dml_dt %>%
    tidyr::separate_rows(genomic_context, sep = ";") %>%
    dplyr::group_by(genomic_context) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(percentage = count / sum(count) * 100)

  p <- ggplot(genomic_context_summary, aes(x = "", y = percentage, fill = genomic_context)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = "Distribution of DMLs across Genomic Features")

  safe_plot(file.path(output_dir, "genomic_context_distribution.svg"), function() {
    print(p)
  })
  flog.info("Analysis complete")
  return(dml_dt)
}

# Perform analysis with provided parameters
flog.info(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f, min.CpG = %d, min.len = %d, dis.merge = %d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
))
result <- perform_dml_analysis(combined_bsseq, base_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  min.CpG = min.CpG,
  min.len = min.len,
  dis.merge = dis.merge,
  cl = cl
)

# Save the result
saveRDS(result, file.path(base_dir, sprintf(
  "delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
), "dml_results.rds"))

# Stop the cluster
stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:")
flog.info(sprintf("Total DMLs found: %d", nrow(result)))
flog.info(sprintf("Hypomethylated DMLs in tumor: %d", sum(result$hypo_in_tumour)))
flog.info(sprintf("Hypermethylated DMLs in tumor: %d", sum(!result$hypo_in_tumour)))

# Calculate mean methylation difference
mean_diff <- mean(result$diff.meth)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff))

# Calculate median methylation difference
median_diff <- median(result$diff.meth)
flog.info(sprintf("Median methylation difference: %.4f", median_diff))

end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))

# Restore warning level at the end of the script
options(warn = 0)
