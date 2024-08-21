# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw
# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.2 0.05 0.01 4 50 50 raw_with_liver

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
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite", "pheatmap")

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
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_dmr_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", delta, p.threshold, fdr.threshold))

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
  sample_list <- lapply(sample_files, readRDS)

  # Create the BSseq objects for all samples
  bsseq_list <- lapply(seq_along(sample_list), function(i) {
    sample_data <- sample_list[[i]]
    BSseq(
      chr = sample_data$chr,
      pos = sample_data$pos,
      M = as.matrix(sample_data$X),
      Cov = as.matrix(sample_data$N),
      sampleNames = paste0(prefix, "_", i)
    )
  })

  # Combine all BSseq objects into one
  combined_bsseq <- do.call(bsseq::combine, bsseq_list)
  return(combined_bsseq)
}


# Load the data
tumour_bsseq <- load_and_create_bsseq(base_dir, "tumour")
control_bsseq <- load_and_create_bsseq(base_dir, "control")

# Combine tumour and control BSseq objects
combined_bsseq <- bsseq::combine(tumour_bsseq, control_bsseq)

saveRDS(combined_bsseq, file.path(base_dir, "combined_bsseq.rds"))
gc()

analyze_areastat_thresholds <- function(top_hypo_dmrs, output_dir) {
  flog.info("Analyzing areaStat distribution and thresholds")

  # Calculate quantiles
  quantiles <- quantile(top_hypo_dmrs$areaStat, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))

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
      p <- ggplot(top_hypo_dmrs, aes(x = areaStat)) +
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

perform_dmr_analysis <- function(combined_bsseq, base_dir, delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge, cl) {
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

  # Call DMRs
  flog.info("Calling DMRs")
  dmrs <- callDMR(dml_test,
    delta = delta, p.threshold = p.threshold,
    minlen = min.len, minCG = min.CpG, dis.merge = dis.merge
  )

  # Convert DMRs to data.table
  dmr_dt <- as.data.table(dmrs)

  # Extract p-values from DML test results for the CpGs in our DMRs
  dml_results <- as.data.table(dml_test)
  setkey(dml_results, chr, pos)

  # Function to get minimum p-value for each DMR (to be used in parallel)
  get_min_pval <- function(dmr) {
    pvals <- dml_results[.(dmr$chr, dmr$start:dmr$end), on = .(chr, pos), nomatch = 0]$pval
    if (length(pvals) == 0) {
      return(NA)
    }
    min(pvals, na.rm = TRUE)
  }

  # Add minimum p-value for each DMR (in parallel)
  dmr_dt[, pval := parSapply(cl, split(dmr_dt, 1:nrow(dmr_dt)), get_min_pval)]

  # Calculate FDR
  dmr_dt[, fdr := p.adjust(pval, method = "BH")]

  # Flag hypo-methylated regions in tumour and apply FDR threshold
  dmr_dt[, `:=`(
    hypo_in_tumour = diff.Methy < 0,
    significant_after_fdr = fdr < fdr.threshold
  )]

  # Select top hypomethylated DMRs after FDR correction
  top_hypo_dmrs <- dmr_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmrs, -areaStat)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmrs, output_dir)

  # You can then use these thresholds in your analysis, e.g.:
  top_hypo_dmrs[, hypomethylation_strength := case_when(
    areaStat < thresholds$very_strong ~ "Very Strong",
    areaStat < thresholds$strong ~ "Strong",
    areaStat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # Write results to BED file
  fwrite(top_hypo_dmrs[, .(chr, start, end, areaStat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmrs.bed"),
    sep = "\t"
  )

  # Visualizations
  flog.info("Creating visualizations")

  plot_top_DMRs <- function(top_hypo_dmrs, tumour_bsseq, control_bsseq, output_dir, n = 50, ext = 0) {
    dmr_plot_dir <- file.path(output_dir, "strongest_hypomethylated_dmr_plots")
    dir.create(dmr_plot_dir, showWarnings = FALSE, recursive = TRUE)

    strongest_dmrs <- tail(top_hypo_dmrs[order(top_hypo_dmrs$areaStat), ], n)
    n_samples <- ncol(tumour_bsseq) + ncol(control_bsseq)

    # Get the sample names from the BSseq objects
    sample_names <- c(sampleNames(tumour_bsseq), sampleNames(control_bsseq))

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
            par(mfrow = c(n_samples, 1), mar = c(3, 3, 2, 1))

            # Plot all samples in the BSseq object
            for (j in 1:length(sample_names)) {
              showOneDMR(dmr, bsseq = bsseq::subsetBySample(tumour_bsseq, sample_names[j]), ext = ext)
              title(
                main = sprintf(
                  "%s - DMR %d: %s:%d-%d\nStrength: %s, areaStat: %.2f",
                  sample_names[j], i, dmr$chr, dmr$start, dmr$end,
                  dmr$hypomethylation_strength, dmr$areaStat
                ),
                cex.main = 0.9
              )
            }
          },
          error = function(e) {
            flog.error(sprintf("Error plotting DMR %d: %s", i, conditionMessage(e)))
            # Create a simple error plot
            plot(1, type = "n", xlab = "", ylab = "", main = sprintf("Error plotting DMR %d", i))
            text(1, 1, labels = conditionMessage(e), cex = 0.8, col = "red")
          }
        )
      }, width = 12, height = 10)
    }

    flog.info(sprintf("Completed plotting %d strongest hypomethylated DMRs", n))
  }


  # Check BSseq object structure
  flog.info(sprintf(
    "Tumour BSseq object: %d samples, %d features",
    ncol(tumour_bsseq), nrow(tumour_bsseq)
  ))
  flog.info(sprintf(
    "Control BSseq object: %d samples, %d features",
    ncol(control_bsseq), nrow(control_bsseq)
  ))

  # Then call the function
  plot_top_DMRs(top_hypo_dmrs, tumour_bsseq, control_bsseq, output_dir, n = 50)

  # 1. Volcano plot
  if ("pval" %in% names(dmr_dt)) {
    flog.info("Creating volcano plot")
    safe_plot(
      file.path(output_dir, "volcano_plot.svg"),
      function() {
        plot <- ggplot(dmr_dt, aes(x = diff.Methy, y = -log10(pval))) +
          geom_point(aes(color = hypo_in_tumour)) +
          scale_color_manual(
            values = c("TRUE" = "blue", "FALSE" = "red"),
            labels = c("TRUE" = "Hypomethylated in Tumor", "FALSE" = "Hypermethylated in Tumor")
          ) +
          labs(
            title = "Volcano plot of DMRs",
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
      plot <- ggplot(dmr_dt, aes(x = diff.Methy)) +
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

  # 3. DMR length distribution
  flog.info("Creating DMR length distribution plot")
  safe_plot(
    file.path(output_dir, "dmr_length_distribution.svg"),
    function() {
      plot <- ggplot(dmr_dt, aes(x = length)) +
        geom_histogram(binwidth = 50, fill = "lightgreen", color = "black") +
        labs(
          title = "Distribution of DMR Lengths",
          x = "DMR Length (bp)",
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

  # Add chromosome and cumulative position to DMR data
  dmr_dt$chr_num <- as.numeric(sub("chr", "", dmr_dt$chr))
  dmr_dt$cumpos <- dmr_dt$start + chr_sizes$cumpos[match(dmr_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dmr_dt$chr, chr_sizes$chr)]

  # Create the chromosome coverage plot
  p <- ggplot(dmr_dt, aes(x = cumpos, y = diff.Methy, color = diff.Methy < 0)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(
      values = c("FALSE" = "red", "TRUE" = "blue"),
      name = "Methylation State",
      labels = c("FALSE" = "Hypermethylated", "TRUE" = "Hypomethylated")
    ) +
    scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
    labs(
      title = paste(
        "DMRs across chromosomes\n",
        "Hypomethylated:", sum(dmr_dt$diff.Methy < 0),
        "Hypermethylated:", sum(dmr_dt$diff.Methy >= 0)
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
  create_manhattan_plot <- function(dmr_dt, output_dir) {
    flog.info("Creating Manhattan plot")
    safe_plot(
      file.path(output_dir, "manhattan_plot.svg"),
      function() {
        plot <- ggplot(dmr_dt, aes(x = cumpos, y = -log10(pval), color = factor(chr))) +
          geom_point(alpha = 0.8, size = 1) +
          scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
          scale_color_manual(values = rep(c("#1B9E77", "#D95F02"), 12)) +
          labs(
            title = "Manhattan Plot of DMRs",
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
  create_qq_plot <- function(dmr_dt, output_dir) {
    flog.info("Creating Q-Q plot")
    safe_plot(
      file.path(output_dir, "qq_plot.svg"),
      function() {
        observed <- sort(-log10(dmr_dt$pval))
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

  # Heatmap of Top DMRs
  create_top_dmr_heatmap <- function(combined_bsseq, dmr_dt, output_dir, top_n = 100) {
    flog.info("Creating heatmap of top DMRs")
    safe_plot(
      file.path(output_dir, "top_dmr_heatmap.png"), # Changed to PNG for testing
      function() {
        # Check if dmr_dt is not empty
        if (nrow(dmr_dt) == 0) {
          flog.error("dmr_dt is empty. Cannot create heatmap.")
          return(NULL)
        }

        top_dmrs <- head(dmr_dt[order(-abs(areaStat))], n = top_n)

        # Check if top_dmrs is not empty
        if (nrow(top_dmrs) == 0) {
          flog.error("No top DMRs selected. Cannot create heatmap.")
          return(NULL)
        }

        dmr_ranges <- GRanges(
          seqnames = top_dmrs$chr,
          ranges = IRanges(start = top_dmrs$start, end = top_dmrs$end)
        )

        meth_mat <- getMeth(combined_bsseq, regions = dmr_ranges, type = "smooth", what = "perRegion")

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
          main = paste("Methylation Levels of Top", top_n, "DMRs"),
          filename = file.path(output_dir, "top_dmr_heatmap.png")
        ) # Changed to PNG
      }
    )
  }

  create_manhattan_plot(dmr_dt, output_dir)
  create_qq_plot(dmr_dt, output_dir)
  create_top_dmr_heatmap(combined_bsseq, dmr_dt, output_dir)

  # Genomic Context Visualization
  flog.info("Annotating DMRs with genomic context")
  dmr_gr <- GRanges(
    seqnames = dmr_dt$chr,
    ranges = IRanges(start = dmr_dt$start, end = dmr_dt$end),
    diff.Methy = dmr_dt$diff.Methy
  )

  # Get genomic features
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  # Filter out non-standard chromosomes
  standard_chromosomes <- paste0("chr", c(1:22))
  dmr_gr <- dmr_gr[seqnames(dmr_gr) %in% standard_chromosomes]

  promoters <- promoters(txdb)
  genes <- suppressWarnings(genes(txdb, single.strand.genes.only = TRUE))
  exons <- exons(txdb)
  introns <- gaps(exons)

  # Ensure all genomic features are on standard chromosomes
  promoters <- promoters[seqnames(promoters) %in% standard_chromosomes]
  genes <- genes[seqnames(genes) %in% standard_chromosomes]
  exons <- exons[seqnames(exons) %in% standard_chromosomes]
  introns <- introns[seqnames(introns) %in% standard_chromosomes]

  # Annotate DMRs
  dmr_annotation <- data.frame(
    DMR = seq_along(dmr_gr),
    Promoter = overlapsAny(dmr_gr, promoters),
    Gene = overlapsAny(dmr_gr, genes),
    Exon = overlapsAny(dmr_gr, exons),
    Intron = overlapsAny(dmr_gr, introns)
  )

  tryCatch(
    {
      ah <- AnnotationHub()
      enhancers <- ah[["AH46978"]] # GeneHancer enhancers for hg38
      enhancers <- enhancers[seqnames(enhancers) %in% standard_chromosomes]
      dmr_annotation$Enhancer <- overlapsAny(dmr_gr, enhancers)
    },
    error = function(e) {
      flog.warn("Failed to fetch enhancer data. Skipping enhancer annotation.")
      dmr_annotation$Enhancer <- FALSE
    }
  )

  # Add annotation to original data
  dmr_dt$genomic_context <- apply(dmr_annotation[, -1], 1, function(x) {
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

  safe_plot(file.path(output_dir, "genomic_context_distribution.svg"), function() {
    print(p)
  })
  flog.info("Analysis complete")
  return(dmr_dt)
}

# Perform analysis with provided parameters
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

# Save the result
saveRDS(result, file.path(base_dir, sprintf(
  "delta_%.2f_p_%.4f_fdr_%.2f_minCpG_%d_minLen_%d_disMerge_%d",
  delta, p.threshold, fdr.threshold, min.CpG, min.len, dis.merge
), "dmr_results.rds"))

# Stop the cluster
stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:")
flog.info(sprintf("Total DMRs found: %d", nrow(result)))
flog.info(sprintf("Hypomethylated DMRs in tumor: %d", sum(result$hypo_in_tumour)))
flog.info(sprintf("Hypermethylated DMRs in tumor: %d", sum(!result$hypo_in_tumour)))

# Calculate mean methylation difference
mean_diff <- mean(result$diff.Methy)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff))

# Calculate median methylation difference
median_diff <- median(result$diff.Methy)
flog.info(sprintf("Median methylation difference: %.4f", median_diff))

end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))

# Restore warning level at the end of the script
options(warn = 0)
