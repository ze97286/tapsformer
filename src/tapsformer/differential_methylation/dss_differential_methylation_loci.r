# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.4 0.05 0.01 raw
# git pull;clear;Rscript src/tapsformer/differential_methylation/dss_differential_methylation.r 0.4 0.05 0.01 raw_with_liver

start_time <- Sys.time()

# initialise command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript dml_analysis.R <delta> <p.threshold> <fdr.threshold> <suffix>")
}
delta <- as.numeric(args[1])
p.threshold <- as.numeric(args[2])
fdr.threshold <- as.numeric(args[3])
suffix <- args[4]

# load libraries
install_and_load <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

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

options(repos = c(CRAN = "https://cloud.r-project.org"))
bioc_packages <- c(
  "DSS", "GenomicRanges", "bsseq", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub"
)
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite", "pheatmap", "reshape")
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)
sessionInfo()

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

options(warn = -1)
suppressMessages(library(DSS))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(bsseq))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(AnnotationHub))

# initialise dirs
base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg", suffix)
log_dir <- file.path("/users/zetzioni/sharedscratch/logs", sprintf("dss_%s_dml_analysis_delta_%.2f_p_%.4f_fdr_%.2f.log", suffix, delta, p.threshold, fdr.threshold))

# set up logging
flog.appender(appender.file(log_dir))
flog.info("Starting optimized DSS DML analysis")

# setup function for safe saving of svg plots to the output dir
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

# data loading
# load all rds files with methylation counters per cpg for tumour and control (=prefix)
load_and_create_bsseq <- function(base_dir, prefix) {
  sample_files <- list.files(path = base_dir, pattern = paste0("^", prefix, "_.*\\.rds$"), full.names = TRUE)
  if (length(sample_files) == 0) {
    stop("Error: No RDS files found with the given prefix.")
  }
  bsseq_list <- lapply(sample_files, function(file_path) {
    sample_data <- readRDS(file_path)
    sample_name <- gsub("\\.rds$", "", basename(file_path))
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

tumour_bsseq <- load_and_create_bsseq(base_dir, "tumour")
control_bsseq <- load_and_create_bsseq(base_dir, "control")
combined_bsseq <- bsseq::combine(tumour_bsseq, control_bsseq)
gc()

# a function for tagging dmls based on area stat quantiles. Saves a histogram of DMLs and their strength tag.
analyze_areastat_thresholds <- function(top_hypo_dmls, output_dir) {
  flog.info("Analyzing areaStat distribution and thresholds")
  quantiles <- quantile(top_hypo_dmls$areaStat, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
  flog.info("areaStat quantiles:")
  moderate_threshold <- quantiles["75%"]
  strong_threshold <- quantiles["90%"]
  very_strong_threshold <- quantiles["99%"]
  flog.info(sprintf("Suggested thresholds for hypomethylation strength:"))
  flog.info(sprintf("Moderate: < %.2f", moderate_threshold))
  flog.info(sprintf("Strong: < %.2f", strong_threshold))
  flog.info(sprintf("Very strong: < %.2f", very_strong_threshold))
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
  return(list(
    moderate = moderate_threshold,
    strong = strong_threshold,
    very_strong = very_strong_threshold
  ))
}

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

  flog.info("Performing DML test")
  group1 <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
  group2 <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
  dml_test <- DMLtest(combined_bsseq, group1 = group1, group2 = group2, smoothing = TRUE)

  flog.info("Calling DMLs")
  dmls <- callDML(dml_test,
    delta = delta,
    p.threshold = p.threshold
  )

  # FDR correction manually to the p-values
  dmls$fdr <- p.adjust(dmls$pval, method = "BH")
  dml_dt <- as.data.table(dmls)

  # identify hypomethylated regions in the tumour and check significance
  dml_dt[, `:=`(
    hypo_in_tumour = diff.meth < 0,
    significant_after_fdr = fdr < fdr.threshold
  )]

  # select top hypomethylated DMLs
  top_hypo_dmls <- dml_dt[hypo_in_tumour == TRUE & significant_after_fdr == TRUE]
  setorder(top_hypo_dmls, -areaStat)
  thresholds <- analyze_areastat_thresholds(top_hypo_dmls, output_dir)

  # tag differential methylation strength by quantile of areastat
  top_hypo_dmls[, hypomethylation_strength := case_when(
    areaStat < thresholds$very_strong ~ "Very Strong",
    areaStat < thresholds$strong ~ "Strong",
    areaStat < thresholds$moderate ~ "Moderate",
    TRUE ~ "Weak"
  )]

  # save output
  fwrite(top_hypo_dmls[, .(chr, start, end, areaStat, pval, fdr, hypomethylation_strength)],
    file.path(output_dir, "hypomethylated_dmls.bed"),
    sep = "\t"
  )

  flog.info("Creating visualizations")
  strongest_dmls <- tail(top_hypo_dmls[order(top_hypo_dmls$areaStat), ], 12)
  plot_top_DMLs(strongest_dmls, combined_bsseq, output_dir)

  # volcano plot
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

  # methylation difference distribution
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

  # chromosome coverage plot
  chr_sizes <- data.frame(
    chr = paste0("chr", c(1:22, "X", "Y")),
    size = c(
      248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
      145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
      101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
      156040895, 57227415
    )
  )
  chr_sizes$cumpos <- cumsum(as.numeric(chr_sizes$size))
  chr_sizes$pos <- chr_sizes$cumpos - chr_sizes$size / 2
  dml_dt$chr_num <- as.numeric(sub("chr", "", dml_dt$chr))
  dml_dt$cumpos <- dml_dt$start + chr_sizes$cumpos[match(dml_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dml_dt$chr, chr_sizes$chr)]
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

  # manhattan plot
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

  create_manhattan_plot(dml_dt, output_dir)
  create_qq_plot(dml_dt, output_dir)

  # genomic context visualization
  flog.info("Annotating DMLs with genomic context")
  dml_gr <- GRanges(
    seqnames = dml_dt$chr,
    ranges = IRanges(start = dml_dt$start, end = dml_dt$end),
    diff.meth = dml_dt$diff.meth
  )
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  standard_chromosomes <- paste0("chr", c(1:22))
  dml_gr <- dml_gr[seqnames(dml_gr) %in% standard_chromosomes]
  promoters <- promoters(txdb)
  genes <- suppressWarnings(genes(txdb, single.strand.genes.only = TRUE))
  exons <- exons(txdb)
  introns <- gaps(exons)
  promoters <- promoters[seqnames(promoters) %in% standard_chromosomes]
  genes <- genes[seqnames(genes) %in% standard_chromosomes]
  exons <- exons[seqnames(exons) %in% standard_chromosomes]
  introns <- introns[seqnames(introns) %in% standard_chromosomes]
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
  dml_dt$genomic_context <- apply(dml_annotation[, -1], 1, function(x) {
    paste(names(x)[x], collapse = ";")
  })
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

# perform analysis with provided parameters
flog.info(sprintf(
  "Starting analysis with delta = %.2f, p.threshold = %.4f, fdr.threshold = %.2f",
  delta, p.threshold, fdr.threshold
))
result <- perform_dml_analysis(combined_bsseq, base_dir,
  delta = delta,
  p.threshold = p.threshold,
  fdr.threshold = fdr.threshold,
  cl = cl
)

stopCluster(cl)

# Result Summary
flog.info("Analysis Results Summary:")
flog.info(sprintf("Total DMLs found: %d", nrow(result)))
flog.info(sprintf("Hypomethylated DMLs in tumor: %d", sum(result$hypo_in_tumour)))
flog.info(sprintf("Hypermethylated DMLs in tumor: %d", sum(!result$hypo_in_tumour)))
mean_diff <- mean(result$diff.meth)
flog.info(sprintf("Mean methylation difference: %.4f", mean_diff))
median_diff <- median(result$diff.meth)
flog.info(sprintf("Median methylation difference: %.4f", median_diff))
end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
options(warn = 0)
