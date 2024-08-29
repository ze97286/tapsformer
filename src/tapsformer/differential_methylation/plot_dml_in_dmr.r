start_time <- Sys.time()

source("dss_common.r")

main <- function() {
    # Set paths
    base_dir <- "/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg/subset3/"
    output_dir <- file.path(base_dir, "dmr_delta_0.50_p_0.0100_fdr_0.01_minCpG_4_minLen_50_disMerge_50_smooth_cv_window","dml_plots")
    combined_bsseq_file <- file.path(base_dir, "combined_bsseq.rds")
    dml_bed_file <- file.path(base_dir, "top_dmls/dmls_in_high_confidence_hypomethylated_dmrs.bed")

    # Load combined BSseq object
    combined_bsseq <- readRDS(combined_bsseq_file)

    # Load DML BED file
    dml_data <- fread(dml_bed_file)
    setnames(dml_data, c("chr", "start", "end", "name", "score", "strand", "pval", "stat", "areaStat", "dmr_id"))

    # Convert to GRanges for easier manipulation
    dml_gr <- GRanges(
        seqnames = dml_data$chr,
        ranges = IRanges(start = dml_data$start + 1, end = dml_data$end),
        mcols = dml_data[, .(name, score, pval, stat, areaStat, dmr_id)]
    )

    # Convert back to data.table with correct format for plot_dmls_per_dmr function
    filtered_dmls_in_dmrs <- as.data.table(dml_gr)
    setnames(filtered_dmls_in_dmrs, "start", "pos")

    # Call the plotting function
    plot_dmls_per_dmr(filtered_dmls_in_dmrs, combined_bsseq, output_dir, top_n_dmrs = 20)
}

# Run the main function
main()
