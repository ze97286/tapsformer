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
    "DSS", "GenomicRanges", "bsseq", "org.Hs.eg.db", "bumphunter", "dmrseq","ComplexHeatmap",
    "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub", "BiocParallel", "Gviz", "edgeR", "limma", "DMRcate", "SummarizedExperiment"
)
cran_packages <- c("cluster", "Rtsne", "jpeg", "minfi", "Biobase", "data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite", "pheatmap", "grid", "gridExtra")
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)
sessionInfo()

options(warn = -1)
suppressMessages(library(DSS))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(bsseq))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(AnnotationHub))

test_logging <- function() {
    print("This is a test log message.")
}

# Setup function for safe saving of svg plots to the output dir
safe_plot <- function(filename, plot_func, width = 10, height = 8) {
    tryCatch(
        {
            svglite::svglite(filename, width = width, height = height)
            plot_func()
            dev.off()
            print(paste("Plot saved as", filename))
        },
        error = function(e) {
            if (dev.cur() > 1) dev.off() # Close device if open
            print(paste("Error creating plot:", filename, "-", conditionMessage(e)))
        }
    )
}

# Load tumour/control samples data
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

    combined_bsseq <- do.call(combineList, bsseq_list)
    return(combined_bsseq)
}

min_coverage_threshold <- function(group_samples) {
    return(max(floor(length(group_samples) / 2), 2))
}

# Load the combined bsseq data set with tumour and control samples.
load_and_combine_bsseq <- function(base_dir, tumour_prefix, control_prefix) {
    tumour_bsseq <- load_and_create_bsseq(base_dir, tumour_prefix)
    control_bsseq <- load_and_create_bsseq(base_dir, control_prefix)
    combined_bsseq <- bsseq::combine(tumour_bsseq, control_bsseq)
    coverage_matrix <- getCoverage(combined_bsseq)
    control_samples <- grep("control_", sampleNames(combined_bsseq), value = TRUE)
    tumour_samples <- grep("tumour_", sampleNames(combined_bsseq), value = TRUE)
    threshold <- min_coverage_threshold(control_samples)
    loci_to_keep <- rowSums(coverage_matrix[, tumour_samples] >= 1) >= threshold &
        rowSums(coverage_matrix[, control_samples] >= 1) >= threshold
    filtered_bsseq <- combined_bsseq[loci_to_keep, ]
    return(filtered_bsseq)
}

# a function for tagging dmls based on area stat quantiles. Saves a histogram of DMLs and their strength tag.
analyze_areastat_thresholds <- function(top_hypo_dmxs, column_name, output_dir) {
    quantiles <- quantile(top_hypo_dmxs[[column_name]], probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
    moderate_threshold <- quantiles["75%"]
    strong_threshold <- quantiles["90%"]
    very_strong_threshold <- quantiles["99%"]
    safe_plot(
        file.path(output_dir, paste0(column_name, "_distribution.svg")),
        function() {
            p <- ggplot(top_hypo_dmxs, aes_string(x = column_name)) +
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
                    title = paste("Distribution of", column_name, "Values"),
                    x = column_name,
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

# Volcano plot
create_volcano_plot <- function(dmx_dt, diff_col, pval_col, output_dir, prefix = "dss") {
    if (pval_col %in% names(dmx_dt)) {
        print("Creating volcano plot")
        safe_plot(
            file.path(output_dir, prefix, "volcano_plot.svg"),
            function() {
                plot <- ggplot(dmx_dt, aes_string(x = diff_col, y = paste0("-log10(", pval_col, ")"))) +
                    geom_point(aes(color = hypo_in_tumour)) +
                    scale_color_manual(
                        values = c("TRUE" = "blue", "FALSE" = "red"),
                        labels = c("TRUE" = "Hypomethylated in Tumor", "FALSE" = "Hypermethylated in Tumor")
                    ) +
                    labs(
                        title = "Volcano plot",
                        x = "Methylation Difference (Negative = Hypomethylation in Tumor)",
                        y = "-log10(p-value)",
                        color = "Methylation State"
                    ) +
                    theme_minimal()
                print(plot)
            }
        )
        return(TRUE)
    } else {
        print("Skipping volcano plot due to missing p-value column")
        return(FALSE)
    }
}

# Methylation difference distribution plot
create_methylation_diff_plot <- function(dmx_dt, diff_col, output_dir, prefix = "dss") {
    print("Creating methylation difference distribution plot")
    safe_plot(
        file.path(output_dir, prefix, "methylation_difference_distribution.svg"),
        function() {
            plot <- ggplot(dmx_dt, aes_string(x = diff_col)) +
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
    return(TRUE)
}

create_dmr_length_plot <- function(dmx_dt, output_dir, prefix = "dss") {
    print("Creating DMR length distribution plot")
    safe_plot(
        file.path(output_dir, "dmr_length_distribution.svg"),
        function() {
            plot <- ggplot(dmx_dt, aes(x = length)) +
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
    return(TRUE) # Indicate success
}

# dmx by chromosome plot
create_chromosome_coverage_plot <- function(dmx_dt, diff_col, output_dir, prefix) {
    print("Creating chromosome coverage plot")

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

    if (!("start" %in% colnames(dmx_dt)) || !("end" %in% colnames(dmx_dt))) {
        dmx_dt[, start := pos]
        dmx_dt[, end := pos + 2]
    }

    # Add cumulative position
    chr_sizes$cumpos <- cumsum(as.numeric(chr_sizes$size))
    chr_sizes$pos <- chr_sizes$cumpos - chr_sizes$size / 2

    # Add chromosome and cumulative position to DMR/DML data
    dmx_dt$chr_num <- as.numeric(sub("chr", "", dmx_dt$chr))
    dmx_dt$cumpos <- dmx_dt$start + chr_sizes$cumpos[match(dmx_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dmx_dt$chr, chr_sizes$chr)]

    # Create the chromosome coverage plot
    p <- ggplot(dmx_dt, aes_string(x = "cumpos", y = diff_col, color = paste0(diff_col, " < 0"))) +
        geom_point(alpha = 0.5) +
        scale_color_manual(
            values = c("FALSE" = "red", "TRUE" = "blue"),
            name = "Methylation State",
            labels = c("FALSE" = "Hypermethylated", "TRUE" = "Hypomethylated")
        ) +
        scale_x_continuous(label = chr_sizes$chr, breaks = chr_sizes$pos) +
        labs(
            title = paste(
                "Methylation differences across chromosomes\n",
                "Hypomethylated:", sum(dmx_dt[[diff_col]] < 0),
                "Hypermethylated:", sum(dmx_dt[[diff_col]] >= 0)
            ),
            x = "Chromosome",
            y = "Methylation Difference (Negative = Hypomethylation in Tumor)"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    safe_plot(file.path(output_dir, prefix, "chromosome_coverage_plot.svg"), function() {
        print(p)
    })

    return(TRUE) # Indicate success
}

# Manhattan plot
create_manhattan_plot <- function(dmx_dt, output_dir, prefix = "dss") {
    print("Creating Manhattan plot")
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

    if (!("start" %in% colnames(dmx_dt)) || !("end" %in% colnames(dmx_dt))) {
        dmx_dt[, start := pos]
        dmx_dt[, end := pos + 2]
    }

    # Add cumulative position
    chr_sizes$cumpos <- cumsum(as.numeric(chr_sizes$size))
    chr_sizes$pos <- chr_sizes$cumpos - chr_sizes$size / 2
    dmx_dt$chr_num <- as.numeric(sub("chr", "", dmx_dt$chr))
    dmx_dt$cumpos <- dmx_dt$start + chr_sizes$cumpos[match(dmx_dt$chr, chr_sizes$chr)] - chr_sizes$size[match(dmx_dt$chr, chr_sizes$chr)]

    safe_plot(
        file.path(output_dir, prefix, "manhattan_plot.svg"),
        function() {
            plot <- ggplot(dmx_dt, aes(x = cumpos, y = -log10(pval), color = factor(chr))) +
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
create_qq_plot <- function(dmx_dt, output_dir, prefix = "dss") {
    print("Creating Q-Q plot")
    safe_plot(
        file.path(output_dir, prefix, "qq_plot.svg"),
        function() {
            observed <- sort(-log10(dmx_dt$pval))
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

# Pie chart of genomic context of differentially methylated regions/loci
create_genomic_context_visualization <- function(dmx_dt, diff_col, output_dir, prefix = "dss") {
    print("Annotating regions with genomic context")
    if (!("start" %in% colnames(dmx_dt)) || !("end" %in% colnames(dmx_dt))) {
        dmx_dt[, start := pos]
        dmx_dt[, end := pos + 2]
    }

    # Convert DMR/DML data to GRanges object
    dmx_gr <- GRanges(
        seqnames = dmx_dt$chr,
        ranges = IRanges(start = dmx_dt$start, end = dmx_dt$end),
        diff = dmx_dt[[diff_col]]
    )

    # Get genomic features
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    # Filter out non-standard chromosomes
    standard_chromosomes <- paste0("chr", c(1:22))
    dmx_gr <- dmx_gr[seqnames(dmx_gr) %in% standard_chromosomes]

    promoters <- promoters(txdb)
    genes <- suppressWarnings(genes(txdb, single.strand.genes.only = TRUE))
    exons <- exons(txdb)
    introns <- gaps(exons)

    # Ensure all genomic features are on standard chromosomes
    promoters <- promoters[seqnames(promoters) %in% standard_chromosomes]
    genes <- genes[seqnames(genes) %in% standard_chromosomes]
    exons <- exons[seqnames(exons) %in% standard_chromosomes]
    introns <- introns[seqnames(introns) %in% standard_chromosomes]

    # Annotate regions (DMRs/DMLs)
    dmx_annotation <- data.frame(
        Region = seq_along(dmx_gr),
        Promoter = overlapsAny(dmx_gr, promoters),
        Gene = overlapsAny(dmx_gr, genes),
        Exon = overlapsAny(dmx_gr, exons),
        Intron = overlapsAny(dmx_gr, introns)
    )

    tryCatch(
        {
            ah <- AnnotationHub()
            enhancers <- ah[["AH46978"]] # GeneHancer enhancers for hg38
            enhancers <- enhancers[seqnames(enhancers) %in% standard_chromosomes]
            dmx_annotation$Enhancer <- overlapsAny(dmx_gr, enhancers)
        },
        error = function(e) {
            print("Failed to fetch enhancer data. Skipping enhancer annotation.")
            dmx_annotation$Enhancer <- FALSE
        }
    )

    # Add annotation to original data
    dmx_dt$genomic_context <- apply(dmx_annotation[, -1], 1, function(x) {
        paste(names(x)[x], collapse = ";")
    })

    genomic_context_summary <- dmx_dt %>%
        tidyr::separate_rows(genomic_context, sep = ";") %>%
        dplyr::group_by(genomic_context) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::mutate(percentage = count / sum(count) * 100)

    p <- ggplot(genomic_context_summary, aes(x = "", y = percentage, fill = genomic_context)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void() +
        labs(title = "Distribution of Regions across Genomic Features")

    safe_plot(file.path(output_dir, prefix, "genomic_context_distribution.svg"), function() {
        print(p)
    })

    return(list(dmx_dt = dmx_dt, plot = p))
}


showOneDMRTwoPrefixes <- function(OneDMR, BSobj, prefix1, prefix2, ext = 500, ylim = c(0, 1)) {
    ## get chr, position, and counts
    allchr <- as.character(seqnames(BSobj))
    allpos <- start(BSobj)
    X <- getBSseq(BSobj, "M")
    N <- getBSseq(BSobj, "Cov")

    ## locate the data for plotting
    chr <- as.character(OneDMR$chr)
    ix.chr <- which(allchr == chr)
    thispos <- allpos[ix.chr]
    thisN <- N[ix.chr, ]
    thisX <- X[ix.chr, ]
    xlim <- c(OneDMR$start - ext, OneDMR$end + ext)
    ix1 <- which(thispos <= xlim[2] & thispos >= xlim[1])

    ## separate samples based on prefixes
    sNames <- sampleNames(BSobj)
    prefix1_samples <- grep(paste0("^", prefix1), sNames, value = TRUE)
    prefix2_samples <- grep(paste0("^", prefix2), sNames, value = TRUE)

    nSample1 <- length(prefix1_samples)
    nSample2 <- length(prefix2_samples)
    nSample <- max(nSample1, nSample2)

    # Set layout to 2 columns, with rows equal to the max number of samples in either group
    layout(matrix(c(1:nSample, (nSample + 1):(2 * nSample)), nrow = nSample, ncol = 2))

    # Adjust margins based on the number of samples
    mar_vertical <- max(0.5, 2.5 / sqrt(nSample))
    mar_horizontal <- max(0.5, 2.5 / sqrt(2)) # Always 2 columns
    par(
        mar = c(mar_vertical, mar_horizontal, mar_vertical, mar_horizontal),
        mgp = c(1.5, 0.5, 0),
        oma = c(2, 2, 2, 2)
    ) # Add outer margins for labels

    thisP <- thisX / thisN

    plotSample <- function(sample_name) {
        i <- which(sNames == sample_name)
        plot(thispos[ix1], thisP[ix1, i],
            type = "h", col = "blue", axes = FALSE, lwd = 1.5,
            xlab = "", ylab = "", ylim = ylim, xlim = xlim,
            main = ""
        )
        box(col = "black")
        axis(1, cex.axis = 0.8)
        axis(2, col = "blue", col.axis = "blue", cex.axis = 0.8)
        title(main = sample_name, cex.main = 0.9)

        thisN.norm <- thisN[ix1, i] / max(thisN[ix1, ]) * ylim[2]
        lines(thispos[ix1], thisN.norm, type = "l", col = "gray", lwd = 1.5)
        axis(
            side = 4, at = seq(0, ylim[2], length.out = 5),
            labels = round(seq(0, max(thisN[ix1, ]), length.out = 5)),
            cex.axis = 0.8
        )

        rect(OneDMR$start, ylim[1], OneDMR$end, ylim[2], col = "#FF00001A", border = NA)
    }

    # Plot all samples in their respective columns
    for (i in 1:nSample) {
        if (i <= nSample1) {
            plotSample(prefix1_samples[i])
        } else {
            plot.new()
        }
    }
    for (i in 1:nSample) {
        if (i <= nSample2) {
            plotSample(prefix2_samples[i])
        } else {
            plot.new()
        }
    }

    # Add overall title and axis labels
    mtext(chr, side = 1, outer = TRUE, line = 0.5, cex = 0.8)
    mtext("methyl%", side = 2, outer = TRUE, line = 0.5, col = "blue", cex = 0.8)
    mtext("read depth", side = 4, outer = TRUE, line = 0.5, cex = 0.8)
    mtext(sprintf("DMR: %s:%d-%d", OneDMR$chr, OneDMR$start, OneDMR$end),
        side = 3, outer = TRUE, line = 0.5, cex = 1
    )
}

# The plot_single_dmr and plot_top_DMRs functions remain the same
plot_single_dmr <- function(filename, dmr, combined_bsseq, i, ext = 0) {
    print(sprintf("Processing DMR %d: %s:%d-%d with ext = %d", i, dmr$chr, dmr$start, dmr$end, ext))

    sNames <- sampleNames(combined_bsseq)
    nSample1 <- length(grep("^tumour", sNames))
    nSample2 <- length(grep("^control", sNames))
    nSample <- max(nSample1, nSample2)

    # Dynamically set plot dimensions
    plot_width <- max(10, 5 * sqrt(2 * nSample)) # Adjusted for 2 columns
    plot_height <- max(8, 4 * sqrt(nSample))

    tryCatch(
        {
            svglite::svglite(filename, width = plot_width, height = plot_height)
            showOneDMRTwoPrefixes(dmr, combined_bsseq, "tumour", "control", ext = ext)
            dev.off()
            print(sprintf("Plot saved to %s", filename))
        },
        error = function(e) {
            print(sprintf("Error encountered while processing DMR %d: %s", i, conditionMessage(e)))
        }
    )
}

cross_validate_dmls <- function(bsseq_data, group1, group2, n_iterations = 5, subsample_fraction = 0.8,
                                delta, p.threshold, fdr.threshold, smoothing) {
    all_dmls <- list()

    for (i in 1:n_iterations) {
        # Subsample from each group separately
        subsample1 <- sample(group1, length(group1) * subsample_fraction)
        subsample2 <- sample(group2, length(group2) * subsample_fraction)
        subsample <- c(subsample1, subsample2)

        sub_bsseq <- bsseq_data[, subsample]

        dml_test <- DMLtest(sub_bsseq, group1 = subsample1, group2 = subsample2, smoothing = smoothing)
        dmls <- callDML(dml_test, delta = delta, p.threshold = p.threshold)

        dml_dt <- as.data.table(dmls)
        dml_dt[, significant_after_fdr := fdr < fdr.threshold]

        all_dmls[[i]] <- dml_dt[significant_after_fdr == TRUE, .(chr, pos)]
    }

    # Keep DMLs that appear in at least half of the iterations
    dml_counts <- table(unlist(lapply(all_dmls, function(x) paste(x$chr, x$pos))))
    consistent_dmls <- names(dml_counts[dml_counts >= n_iterations / 2])

    return(consistent_dmls)
}

sliding_window_filter <- function(dmls, window_size, min_cpgs = 3, consistency_threshold = 0.8) {
  dmls[, window := cut(pos, breaks = seq(min(pos), max(pos) + window_size, by = window_size))]
  windowed_dmls <- dmls[, .(
    mean_diff = mean(diff),
    mean_pval = mean(pval),
    n_dmls = .N,
    consistency = mean(sign(diff) == sign(mean(diff)))
  ), by = .(chr, window)]

  significant_windows <- windowed_dmls[
    abs(mean_diff) >= delta &
      mean_pval < p.threshold &
      n_dmls >= min_cpgs &
      consistency >= consistency_threshold
  ]

  dmls[chr %in% significant_windows$chr & window %in% significant_windows$window]
}

# generate a plot of top DMLs and their methylation levels in each sample
plot_top_DMLs <- function(top_hypo_dmls, combined_bsseq, output_dir) {
  # Get the raw methylation data
  methylation_data <- bsseq::getMeth(combined_bsseq, type = "raw")
  sample_names <- colnames(methylation_data)
  plot_data <- data.table(chr = top_hypo_dmls$chr, pos = top_hypo_dmls$pos)

  for (i in 1:nrow(top_hypo_dmls)) {
    dml <- top_hypo_dmls[i, ]

    # Find the corresponding methylation levels
    matching_indices <- which(seqnames(combined_bsseq) == dml$chr & start(combined_bsseq) == dml$pos)

    # Extract the methylation levels for this DML
    meth_levels <- methylation_data[matching_indices, , drop = FALSE]

    # Check if meth_levels has valid data
    if (is.null(meth_levels) || nrow(meth_levels) == 0 || ncol(meth_levels) == 0) {
      print(sprintf("No valid methylation data found for DML at %s:%d", dml$chr, dml$pos))
      next # Skip to the next DML if no valid data is found
    }

    # Add methylation levels to plot_data with sample names as columns
    plot_data[i, (sample_names) := as.list(meth_levels)]
  }

  # Melt the data for plotting
  plot_data <- melt(plot_data,
    id.vars = c("chr", "pos"),
    variable.name = "Sample", value.name = "MethylationLevel"
  )

  # Add a new column to identify tumour and control samples
  plot_data[, SampleType := ifelse(grepl("^tumour", Sample), "Tumour", "Control")]

  # Check if plot_data is empty
  if (nrow(plot_data) == 0) {
    print("Plot data is empty after processing. No plot will be generated.")
    return(NULL)
  }

  # Generate and save the plot directly
  output_filename <- file.path(output_dir, "dml_methylation_plot.svg")

  tryCatch(
    {
      # Calculate optimal number of columns based on the number of loci
      num_loci <- nrow(top_hypo_dmls)
      num_samples <- length(sample_names)
      ncol <- min(5, num_loci) 

      # Calculate the number of rows required
      nrow <- ceiling(num_loci / ncol)

      # Adjust plot dimensions based on the number of rows and columns
      plot_width <- max(18, ncol * 2) # Increase width to scale with more columns
      plot_height <- max(10, nrow * 1.5) # Increase height to scale with more rows

      # Generate the plot
      svglite::svglite(output_filename, width = plot_width, height = plot_height)

      plot <- ggplot(plot_data, aes(x = Sample, y = MethylationLevel, shape = SampleType, color = SampleType)) +
        geom_point(size = 2) +
        facet_wrap(~ chr + pos, scales = "free_y", ncol = ncol) + # Use the updated ncol
        theme_bw() +
        labs(
          title = "Methylation Levels Across Samples for Each DML",
          x = "Sample", y = "Methylation Level"
        ) +
        scale_shape_manual(values = c(Tumour = 16, Control = 1)) +
        scale_color_manual(values = c(Tumour = "blue", Control = "red")) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), # Adjust size for better readability
          strip.text = element_text(size = 12), # Adjust if needed
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          panel.spacing = unit(1.5, "lines") # Increase spacing between panels
        )

      print(plot)
      dev.off()
    },
    error = function(e) {
      # Handle any errors that occur during plotting
      print(sprintf("Error generating plot: %s", conditionMessage(e)))
      dev.off() # Ensure device is closed even if there's an error
    }
  )
  return(output_filename)
}