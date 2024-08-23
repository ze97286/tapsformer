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
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "tidyr", "ggplot2", "svglite", "pheatmap")
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
    flog.info("This is a test log message.", name = "dss_logger")
}

# Setup function for safe saving of svg plots to the output dir
safe_plot <- function(filename, plot_func, width = 10, height = 8) {
    tryCatch(
        {
            svglite::svglite(filename, width = width, height = height)
            plot_func()
            dev.off()
            flog.info(paste("Plot saved as", filename), name = "dss_logger")
        },
        error = function(e) {
            if (dev.cur() > 1) dev.off() # Close device if open
            flog.error(paste("Error creating plot:", filename, "-", conditionMessage(e)), name = "dss_logger")
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
    threshold <- min_coverage_threshold(control_samples)
    loci_to_keep <- rowSums(coverage_matrix >= 1) >= threshold
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
create_volcano_plot <- function(dmx_dt, diff_col, pval_col, output_dir) {
    if (pval_col %in% names(dmx_dt)) {
        flog.info("Creating volcano plot", name = "dss_logger")
        safe_plot(
            file.path(output_dir, "volcano_plot.svg"),
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
        flog.warn("Skipping volcano plot due to missing p-value column", name = "dss_logger")
        return(FALSE)
    }
}

# Methylation difference distribution plot
create_methylation_diff_plot <- function(dmx_dt, diff_col, output_dir) {
    flog.info("Creating methylation difference distribution plot", name = "dss_logger")
    safe_plot(
        file.path(output_dir, "methylation_difference_distribution.svg"),
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

create_dmr_length_plot <- function(dmx_dt, output_dir) {
    flog.info("Creating DMR length distribution plot", name = "dss_logger")
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
create_chromosome_coverage_plot <- function(dmx_dt, diff_col, output_dir) {
    flog.info("Creating chromosome coverage plot", name = "dss_logger")

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

    safe_plot(file.path(output_dir, "chromosome_coverage_plot.svg"), function() {
        print(p)
    })

    return(TRUE) # Indicate success
}

# Manhattan plot
create_manhattan_plot <- function(dmx_dt, output_dir) {
    flog.info("Creating Manhattan plot", name = "dss_logger")
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
        file.path(output_dir, "manhattan_plot.svg"),
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
create_qq_plot <- function(dmx_dt, output_dir) {
    flog.info("Creating Q-Q plot", name = "dss_logger")
    safe_plot(
        file.path(output_dir, "qq_plot.svg"),
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
create_genomic_context_visualization <- function(dmx_dt, diff_col, output_dir) {
    flog.info("Annotating regions with genomic context", name = "dss_logger")
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
            flog.warn("Failed to fetch enhancer data. Skipping enhancer annotation.", name = "dss_logger")
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

    safe_plot(file.path(output_dir, "genomic_context_distribution.svg"), function() {
        print(p)
    })

    return(list(dmx_dt = dmx_dt, plot = p)) # Return updated dmx_dt and the plot
}

plot_single_dmr <- function(filename, dmr, combined_bsseq, i, ext) {
    # Determine the number of samples in the combined_bsseq object
    num_samples <- length(sampleNames(combined_bsseq))

    # Adjust margins and text size based on the number of samples
    if (num_samples > 20) {
        mar <- c(4, 4, 2, 2) + 0.1 # Smaller margins for larger sample sizes
        cex_main <- 0.6 # Smaller text for large sample sizes
        plot_width <- 20
        plot_height <- 14
    } else if (num_samples > 10) {
        mar <- c(5, 4, 3, 2) + 0.1
        cex_main <- 0.7
        plot_width <- 18
        plot_height <- 12
    } else {
        mar <- c(7, 4, 4, 2) + 0.1 # Larger margins for fewer samples
        cex_main <- 0.8
        plot_width <- 16
        plot_height <- 12
    }

    safe_plot(filename, function() {
        tryCatch(
            {
                par(mar = mar)
                showOneDMR(dmr, combined_bsseq, ext = ext)
                title(
                    main = sprintf(
                        "DMR %d: %s:%d-%d\nStrength: %s, areaStat: %.2f",
                        i, dmr$chr, dmr$start, dmr$end, dmr$hypomethylation_strength, dmr$areaStat
                    ),
                    cex.main = cex_main,
                    line = 2
                )
            },
            error = function(e) {
                flog.error(sprintf("Error plotting DMR %d: %s", i, conditionMessage(e)), name = "dss_logger")
                plot(1, type = "n", xlab = "", ylab = "", main = sprintf("Error plotting DMR %d", i))
                text(1, 1, labels = conditionMessage(e), cex = 0.8, col = "red")
            }
        )
    }, width = plot_width, height = plot_height) # Use dynamically adjusted dimensions
}
